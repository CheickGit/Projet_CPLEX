#include "modele1_lineaire.h"
#include <set> // Ajout de l'inclusion pour std::set
ILOSTLBEGIN
#define CPLEX_EPS 1e-6
const double FRACTIONALITY_TOLERANCE_PC = 1e-6; // Nom différent pour éviter conflit si CPLEX_EPS est global

ILOBRANCHCALLBACK1(MyBranch, IloNumVarArray, vars_y) {
    if ( getBranchType() != BranchOnVariable )
       return;
       
    IloNumArray y;
    IloNumArray obj;
    IntegerFeasibilityArray feas;
 
    try {
       y    = IloNumArray(getEnv());
       obj  = IloNumArray(getEnv());
       feas = IntegerFeasibilityArray(getEnv());
       getValues(y, vars_y);
       getObjCoefs(obj, vars_y);
       getFeasibilities(feas, vars_y);
 
       IloInt bestj  = -1;
       IloNum maxinf = 0.0;
       IloNum maxobj = 0.0;
       IloInt cols = vars_y.getSize();
       for (IloInt j = 0; j < cols; j++) {
          if ( feas[j] == Infeasible ) {
             IloNum yj_inf = y[j] - IloFloor (y[j]);
             if ( yj_inf > 0.5 )
                yj_inf = 1.0 - yj_inf;
             if ( yj_inf >= maxinf                              &&
                  (yj_inf > maxinf || IloAbs (obj[j]) >= maxobj)  ) {
                bestj  = j;
                maxinf = yj_inf;
                maxobj = IloAbs (obj[j]);
             }
          }
       }
 
       if ( bestj >= 0 ) {
          makeBranch(vars_y[bestj], y[bestj], IloCplex::BranchUp,   getObjValue());
          makeBranch(vars_y[bestj], y[bestj], IloCplex::BranchDown, getObjValue());
       }
    }
    catch (...) {
       y.end();
       obj.end();
       feas.end();
       throw;
    }
    y.end();
    obj.end();
    feas.end();
}

bool isInteger(IloNum value) {
    const double tolerance = 1e-5;
    return std::fabs(value - std::round(value)) <= tolerance;
}

ILOBRANCHCALLBACK4(MyHybridBranch, IloNumVarArray, vars_y, std::vector<Pseudocost>*, custom_pseudocosts, IloCplex, cplex_solver, int, max_depth){
    if (getBranchType() != BranchOnVariable) {
        return;
    }

    IloEnv cb_env = getEnv();
    IloNumArray y_values(cb_env);
    IntegerFeasibilityArray feas(cb_env);

    try {
        getValues(y_values, vars_y);
        getFeasibilities(feas, vars_y);

        IloInt best_j_idx = -1;
        double max_branching_score = -std::numeric_limits<double>::infinity(); // On cherche le score le plus élevé

       
        int current_depth = getCurrentNodeDepth(); // Pour logs ou logique future
        bool use_strong = (current_depth <= max_depth);

        if (use_strong) {
            int cpx_call_status;
            CPXENVptr c_api_env = CPXopenCPLEX(&cpx_call_status);
            if (!c_api_env) {
                cb_env.warning() << "MyHybridBranch SB: Could not open CPLEX C API environment." << std::endl;
                y_values.end(); feas.end(); return;
            }

            CPXLPptr c_api_lp = CPXcloneprob(c_api_env, const_cast<CPXLPptr>(cplex_solver.getImpl()->getCplexLp()), &cpx_call_status);
            if (!c_api_lp) {
                cb_env.warning() << "MyHybridBranch SB: Could not clone CPLEX problem." << std::endl;
                CPXcloseCPLEX(&c_api_env);
                y_values.end(); feas.end(); return;
            }

            cpx_call_status = CPXchgprobtype(c_api_env, c_api_lp, CPXPROB_LP);
            if (cpx_call_status) {
                cb_env.warning() << "MyHybridBranch SB: Could not change problem type to LP." << std::endl;
                CPXfreeprob(c_api_env, &c_api_lp);
                CPXcloseCPLEX(&c_api_env);
                y_values.end(); feas.end(); return;
            }

            std::vector<int> candidate_c_api_indices;
            std::vector<IloInt> candidate_vars_y_indices;

            for (IloInt j = 0; j < vars_y.getSize(); ++j) {
                if (feas[j] == Infeasible && !isInteger(y_values[j])) {
                    const char* var_name = vars_y[j].getName();
                    if (!var_name || strlen(var_name) == 0) {
                        char namebuf[64];
                        sprintf(namebuf, "y_sb_var_idx_%d", (int)j);
                        vars_y[j].setName(namebuf);
                        var_name = vars_y[j].getName();
                    }
                    int c_api_idx;
                    if (CPXgetcolindex(c_api_env, c_api_lp, var_name, &c_api_idx) == 0) {
                        candidate_c_api_indices.push_back(c_api_idx);
                        candidate_vars_y_indices.push_back(j);
                    }
                }
            }
            if (candidate_c_api_indices.empty()) {
                CPXfreeprob(c_api_env, &c_api_lp);
                CPXcloseCPLEX(&c_api_env);
                y_values.end(); feas.end(); return;
            }
            std::vector<double> down_obj_values_sb(candidate_c_api_indices.size());
            std::vector<double> up_obj_values_sb(candidate_c_api_indices.size());

            CPXsetintparam(c_api_env, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
            CPXsetintparam(c_api_env, CPXPARAM_ScreenOutput, CPX_OFF);
            int strong_branch_iterations = 10; // Limite d'itérations Simplex

            cpx_call_status = CPXstrongbranch(c_api_env, c_api_lp,
                                      candidate_c_api_indices.data(),
                                      candidate_c_api_indices.size(),
                                      down_obj_values_sb.data(), up_obj_values_sb.data(),
                                      strong_branch_iterations);

            if (cpx_call_status == 0) { // Succès de CPXstrongbranch
                double current_obj_at_node = getObjValue();
                for (size_t i = 0; i < candidate_c_api_indices.size(); ++i) {
                    IloInt original_j_idx = candidate_vars_y_indices[i];
                    double var_val = y_values[original_j_idx];

                    double down_obj_sb = (down_obj_values_sb[i] < 1e+19) ? down_obj_values_sb[i] : current_obj_at_node + 1e7;
                    double up_obj_sb   = (up_obj_values_sb[i] < 1e+19)   ? up_obj_values_sb[i]   : current_obj_at_node + 1e7;

                    double down_gain = down_obj_sb - current_obj_at_node; // Doit être >= 0 pour une minimisation
                    double up_gain   = up_obj_sb   - current_obj_at_node;   // Doit être >= 0

                    // Score de Strong Branching (formule de l'article [cite: 38])
                    double mu_sb = 1.0/6.0; // Comme dans l'article [cite: 39]
                    double sb_score = (1.0 - mu_sb) * std::min(down_gain, up_gain) +
                                       mu_sb  * std::max(down_gain, up_gain);

                    // Mise à jour des pseudocoûts personnalisés
                    if (original_j_idx < (IloInt)(*custom_pseudocosts).size()) {
                        double val_floor = floor(var_val);
                        double val_ceil = ceil(var_val);

                        // Fractionnalité pour la branche BAS (x_i - floor(x_i))
                        double frac_down = var_val - val_floor;
                        // Fractionnalité pour la branche HAUT (ceil(x_i) - x_i)
                        double frac_up = val_ceil - var_val;

                        if (frac_down > FRACTIONALITY_TOLERANCE_PC) {
                            double contrib_down = down_gain / frac_down; // Psi_i^-
                            int& count_d = (*custom_pseudocosts)[original_j_idx].down_count;
                            double& pc_down_score_ref = (*custom_pseudocosts)[original_j_idx].down_score;
                            pc_down_score_ref = (count_d * pc_down_score_ref + contrib_down) / (count_d + 1.0);
                            count_d++;
                        }
                        if (frac_up > FRACTIONALITY_TOLERANCE_PC) {
                            double contrib_up = up_gain / frac_up; // Psi_i^+
                            int& count_u = (*custom_pseudocosts)[original_j_idx].up_count;
                            double& pc_up_score_ref = (*custom_pseudocosts)[original_j_idx].up_score;
                            pc_up_score_ref = (count_u * pc_up_score_ref + contrib_up) / (count_u + 1.0);
                            count_u++;
                        }
                    }
                    if (sb_score > max_branching_score) {
                        max_branching_score = sb_score;
                        best_j_idx = original_j_idx;
                    }
                }
            }
            CPXfreeprob(c_api_env, &c_api_lp);
            CPXcloseCPLEX(&c_api_env);

        } else { 
            double avg_up_pc_val = 0.0, avg_down_pc_val = 0.0;
            int total_up_init_counts = 0, total_down_init_counts = 0;

            // Calculer les moyennes globales des pseudocoûts initialisés [cite: 52, 53]
            for (IloInt j = 0; j < vars_y.getSize(); ++j) {
                if (j < (IloInt)(*custom_pseudocosts).size()){
                    if ((*custom_pseudocosts)[j].up_count > 0) {
                        avg_up_pc_val += (*custom_pseudocosts)[j].up_score; // Somme des Psi_i^+ initialisés
                        total_up_init_counts++;
                    }
                    if ((*custom_pseudocosts)[j].down_count > 0) {
                        avg_down_pc_val += (*custom_pseudocosts)[j].down_score; // Somme des Psi_i^- initialisés
                        total_down_init_counts++;
                    }
                }
            }
            avg_up_pc_val = (total_up_init_counts > 0) ? avg_up_pc_val / total_up_init_counts : 1.0; // Psi_avg^+ [cite: 53]
            avg_down_pc_val = (total_down_init_counts > 0) ? avg_down_pc_val / total_down_init_counts : 1.0; // Psi_avg^-

            for (IloInt j = 0; j < vars_y.getSize(); ++j) {
                if (feas[j] == Infeasible && !isInteger(y_values[j]) && j < (IloInt)(*custom_pseudocosts).size()) {
                    double var_val = y_values[j];
                    double val_floor = floor(var_val);
                    double val_ceil = ceil(var_val);
                     double f_j_down = var_val - val_floor;
                    // f_i^+ = ceil(x_i) - x_i (pour variable binaire ou entière entre 0 et 1, cela devient 1 - (x_i - floor(x_i)))
                    double f_j_up = val_ceil - var_val;
                    // Obtenir Psi_j^- et Psi_j^+ (pseudocoûts par unité) [cite: 52]
                    double psi_j_down = ((*custom_pseudocosts)[j].down_count > 0) ?
                                        (*custom_pseudocosts)[j].down_score : avg_down_pc_val;
                    double psi_j_up = ((*custom_pseudocosts)[j].up_count > 0) ?
                                      (*custom_pseudocosts)[j].up_score : avg_up_pc_val;
                    // Estimations de dégradation d'objectif (q_j^-, q_j^+) [cite: 51] (implicitement, texte avant eq 5)
                    double q_j_down = 0.0;
                    if (f_j_down > FRACTIONALITY_TOLERANCE_PC) q_j_down = f_j_down * psi_j_down;
                    else q_j_down = std::numeric_limits<double>::max(); // Branche infaisable ou négligeable

                    double q_j_up = 0.0;
                    if (f_j_up > FRACTIONALITY_TOLERANCE_PC) q_j_up = f_j_up * psi_j_up;
                    else q_j_up = std::numeric_limits<double>::max(); // Branche infaisable ou négligeable
                    if (q_j_down == std::numeric_limits<double>::max() && q_j_up == std::numeric_limits<double>::max()) {
                        continue; // Ne peut pas brancher sur cette variable
                    }


                    // Calcul du Score Combiné (Équation 3 de l'article [cite: 38])
                    double mu_pc = 1.0/6.0; // Comme utilisé dans l'article [cite: 39]
                    double pc_score = (1.0 - mu_pc) * std::min(q_j_down, q_j_up) +
                                       mu_pc  * std::max(q_j_down, q_j_up);
                    if (pc_score > max_branching_score) {
                        max_branching_score = pc_score;
                        best_j_idx = j;
                    }
                }
            }
        }

        if (best_j_idx >= 0) {
            makeBranch(vars_y[best_j_idx], y_values[best_j_idx], IloCplex::BranchUp, getObjValue());
            makeBranch(vars_y[best_j_idx], y_values[best_j_idx], IloCplex::BranchDown, getObjValue());
        }

    } catch (const IloException& e) {
        std::cerr << "Exception Concert/CPLEX dans MyHybridBranch: " << e << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception std C++ dans MyHybridBranch: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Exception inconnue dans MyHybridBranch." << std::endl;
    }

    y_values.end();
    feas.end();
}
//----------------------------------------------------------------------------------

MCND::MCND(const std::string & instance) : cplex(env), model(env) {
    read_data(instance);
}

//----------------------------------------------------------------------------------

void MCND::set_parameters() {
    cplex.setParam(IloCplex::Param::ClockType, 1);
    //cplex.setParam(IloCplex::Param::Threads, 1);
    //cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Param::TimeLimit, 3600.0); // Time limit in seconds
    
}

//----------------------------------------------------------------------------------

void MCND::create_model() {
    x = IloNumVarArray(env, ndemands * narcs, 0, 1, ILOINT); // Variables de flux
    y = IloNumVarArray(env, narcs, 0, 1, ILOINT); // Variables continues forcées
    IloExpr obj(env);
    for (int a = 0; a < narcs; ++a) {
        obj += arcs[a].f * y[a];
        for (int k = 0; k < ndemands; ++k) {
            obj += arcs[a].c[k] * d_k[k].quantity*x[k * narcs + a];
        }
    }
    model.add(IloMinimize(env, obj)); // Ajout de la fonction objectif
    obj.end();
    
    // Contraintes de conservation des flux
    for (int k = 0; k < ndemands; ++k) {
        for (int i = 1; i <= nnodes; i++) {
            IloExpr constraint(env);
            for (int e = 0; e < narcs; ++e) {
                if (i == arcs[e].i)
                    constraint += x[k * narcs + e]; // Flux sortant
                else if (i == arcs[e].j)
                    constraint -= x[k * narcs + e]; // Flux entrant
            }
            if (i == d_k[k].D)
                constraint += 1; // Demande au nœud destination
            else if (i == d_k[k].O)
                constraint -= 1; // Offre au nœud origine
            model.add(constraint == 0); // Conservation des flux
            constraint.end();
        }
    }

    // Contraintes de capacité globale
    for (int e = 0; e < narcs; ++e) {
        IloExpr constraint(env);
        for (int k = 0; k < ndemands; ++k) {
            constraint += d_k[k].quantity*x[k * narcs + e]; // Somme des flux sur l'arc
        }
        constraint -= arcs[e].capa * y[e]; // Capacité de l'arc activé
        model.add(constraint <= 0);
        constraint.end();
    }

    // Contraintes de capacité individuelle
    for (int k = 0; k < ndemands; ++k) {
        for (int e = 0; e < narcs; ++e) {
            IloExpr constraint(env);
            constraint += d_k[k].quantity*x[k * narcs + e]; // Flux de la demande k sur l'arc e
            constraint -= arcs[e].b[k] * y[e]; // Borne supérieure pour cette demande
            model.add(constraint <= 0);
            constraint.end();
        }
    }

    cplex.extract(model); // Extraire le modèle dans CPLEX
}

//----------------------------------------------------------------------------------

bool MCND::solve(const std::string & instance){
    set_parameters();
    create_model();
    

   // Initialisation
   pseudocosts.resize(narcs);

   cplex.use(MyHybridBranch(env, y, &pseudocosts, cplex, 10));  
   cplex.solve();

    if (cplex.getStatus() == IloAlgorithm::Infeasible) {
        std::cout << "Le problème est infaisable : aucune solution trouvée." << std::endl;
        return false;
    } else if (cplex.getStatus() == IloAlgorithm::Unbounded) {
        std::cout << "Le problème est non borné : aucune solution trouvée." << std::endl;
        return false;
    } else if (cplex.getStatus() == IloAlgorithm::Feasible ||
               cplex.getStatus() == IloAlgorithm::Optimal) {
        std::cout << "Une solution a été trouvée." << std::endl;

        IloNumArray y_(env);
        IloNumArray x_(env);
        cplex.getValues(y_, y); // Récupérer les valeurs des variables y
        cplex.getValues(x_, x); // Récupérer les valeurs des variables x
        std::string inst(instance);
        while(inst.find("/")!=std::string::npos){
            inst = inst.substr(inst.find("/")+1);
        }
        std::ofstream file("fileout", std::ios::app);
        file<<std::setprecision(10)<<instance<<" lb: "<<cplex.getBestObjValue()<<" ub: "<<cplex.getObjValue()<<" gap: "<<cplex.getMIPRelativeGap()*100<<" nodes: "<<cplex.getNnodes()<<" t: "<<cplex.getCplexTime()<<std::endl;
        file.close();

        int nbarcsActifs=0;
        std::vector<int> taille_chemin;
        taille_chemin.resize(ndemands, 0);
        // Afficher les arcs activés
        std::cout << "Arcs activés:\n";
        for (int e = 0; e < narcs; ++e) {
            if (y_[e] > 0.5) {
                ++nbarcsActifs;
                std::cout << "Arc activé : " << arcs[e].i << " -> " << arcs[e].j << std::endl;
            }
        }

        // Afficher les flots pour chaque demande
        std::cout << "\nDemandes transportées:\n";
        for (int k = 0; k < ndemands; ++k) {
            std::cout << "Demande " << k + 1 << " :\n";
            for (int e = 0; e < narcs; ++e) {
                double flow = x_[k * narcs + e];
                if (flow > 1e-6) {
                    taille_chemin[k]++;
                    std::cout << "  Arc " << arcs[e].i << " -> " << arcs[e].j
                              << " : flot = " << flow*d_k[k].quantity << std::endl;
                }
            }
        }

        double val=0;
        for (int a = 0; a < narcs; ++a) {
            val += arcs[a].f * y_[a];
            for (int k = 0; k < ndemands; ++k) {
                val += arcs[a].c[k] *d_k[k].quantity* x_[k * narcs + a];
            }
        }
        std::cout<<std::setprecision(8)<<val<<std::endl;
        ecrire_solution(instance, x_, y_, nbarcsActifs, val, taille_chemin);
        y_.end();
        x_.end();
    }
    return true;
    
}

void MCND::ecrire_solution(const std::string & instance, const IloNumArray& x_, const IloNumArray& y_, int nbarcsActifs, double coutTotal, std::vector<int> taille_chemin) {
    std::string inst(instance);
    size_t lastSlash = inst.find_last_of("/\\");
    size_t lastDot = inst.find_last_of(".");
    
    std::string filename = inst.substr(
        lastSlash == std::string::npos ? 0 : lastSlash + 1,
        lastDot == std::string::npos ? inst.length() : lastDot - (lastSlash == std::string::npos ? 0 : lastSlash + 1)
    );
    
    // Création du dossier RESULTAT s'il n'existe pas
    if (system("mkdir -p RESULTAT") != 0) {
        std::cout << "Erreur lors de la création du dossier RESULTAT" << std::endl;
        return;
    }
    
    std::string fullPath = "RESULTAT/" + filename + ".txt";
    std::ofstream file(fullPath);
    
    if (!file.is_open()) {
        std::cout << "Impossible d'ouvrir le fichier: " << fullPath << std::endl;
        return;
    }

    file << "Nombre d'arcs actifs: " << nbarcsActifs << "\n";
    for (int e = 0; e < narcs; ++e) {
        if (y_[e] > 0.5) {
            file << e+1 << std::endl;
        }
    }
    // Afficher les flots pour chaque demande
    for (int k = 0; k < ndemands; ++k) {
        file << "Demande " << k + 1 << " (Taille: " << taille_chemin[k] << ")\n";
        for (int e = 0; e < narcs; ++e) {
            double flow = x_[k * narcs + e];
            if (flow > 1e-6) {
                file << e+1 << std::endl;
            }
        }
    }
    
    file << "Valeur totale: " << setprecision(8)<<coutTotal<< "\n";
    file.close();
}

//----------------------------------------------------------------------------------

void
MCND::read_data(const std::string & instance) {
    std::ifstream file;
    file.open(instance.c_str());
    if (!file.is_open()) {
        std::cerr << "Failure to open datafile: " << instance << std::endl;
        exit(1); // Utilisation de exit() au lieu de abort() pour une sortie plus propre
    }
    
    
    
    std::string s;
    std::istringstream ss;
    int style=0;
    getline(file,s);
    if(s == "MULTIGEN.DAT:" ||s == " MULTIGEN.DAT:"){
        style=1;
        getline(file,s);
    }
    ss.str(s);
    
    // read number of locations and number of customers
    
    ss>>nnodes;
    ss>>narcs;
    ss>>ndemands;
    
    ss.clear();
    
    d_k.resize(ndemands);
    arcs.resize(narcs);
    if(style == 1){
        double cost=0;
        for(int i=0; i<narcs; ++i){
            getline(file,s);
            ss.str(s);
            ss>>arcs[i].i;
            ss>>arcs[i].j;
            ss>>cost;
            ss>>arcs[i].capa;
            ss>>arcs[i].f;
            ss.clear();
            arcs[i].c.resize(ndemands, cost);
            arcs[i].b.resize(ndemands, arcs[i].capa);
        }
        
        for(int i=0; i<ndemands; ++i){
            getline(file,s);
            ss.str(s);
            ss>> d_k[i].O >> d_k[i].D >> d_k[i].quantity;
            ss.clear();
        }
        //calculate variable bounds
        for(int k=0;k<ndemands;++k){
            for(int e=0; e<narcs; ++e){
                if(d_k[k].quantity < arcs[e].b[k])
                    arcs[e].b[k] = d_k[k].quantity;
            }
        }
    }else if(style == 0){
        int pk=0;
        for(int i=0; i<narcs; ++i){
            getline(file,s);
            ss.str(s);
            ss>>arcs[i].i;
            ss>>arcs[i].j;
            ss>>arcs[i].f;
            ss>>arcs[i].capa;
            ss>>pk;
            ss.clear();
            arcs[i].c.resize(ndemands, 0.0);
            arcs[i].b.resize(ndemands, 0.0);
            for(int k=0;k<ndemands;++k){
                getline(file,s);
                ss.str(s);
                ss>>pk;
                ss>>arcs[i].c[pk-1];
                ss>>arcs[i].b[pk-1];
                ss.clear();
            }
        }
        // demands origin b+, destination b-
        
        int node; double b;
        for(int i=0; i<2*ndemands; ++i){
            getline(file,s);
            ss.str(s);
            ss>> pk >> node >> b;
            if(b<0) d_k[pk-1].D = node;
            else{d_k[pk-1].O = node; d_k[pk-1].quantity = b;}
            ss.clear();
        }
    }
    file.close();
    
}

//----------------------------------------------------------------------------------

MCND::~MCND() {
    
    try {
        arcs.clear();
        d_k.clear();

        //y.end();
        //model.end();
        cplex.clearModel();
        cplex.end();
        env.end();

    } catch (IloException& e) {
        std::cerr << "ERROR: " << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Error" << std::endl;
    }
}
