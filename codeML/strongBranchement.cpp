#include "modele1_lineaire.h"
#include <set> // Ajout de l'inclusion pour std::set
ILOSTLBEGIN
#define CPLEX_EPS 1e-6

/*
ILOBRANCHCALLBACK1(MyBranch, IloNumVarArray, vars_y) {
    if ( getBranchType() != BranchOnVariable )
       return;
 
    // Branch on var with largest objective coefficient
    // among those with largest infeasibility
 
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

// Ajouter en en-tête
#include <vector>
struct PseudocostData {
    double sigma_up = 0.0;
    double sigma_down = 0.0;
    int eta_up = 0;
    int eta_down = 0;
    double psi_up = 1.0;  // Valeur par défaut
    double psi_down = 1.0;
};

// Callback modifiée avec strong-branching et pseudocosts
ILOBRANCHCALLBACK1(MyBranch, IloNumVarArray, vars_y) {
    if (getBranchType() != BranchOnVariable) return;

    IloNumArray y(getEnv());
    IntegerFeasibilityArray feas(getEnv());
    getValues(y, vars_y);
    getFeasibilities(feas, vars_y);

    const double mu = 1.0/6.0; // Paramètre de score
    const int eta_reliable = 8; // Seuil de fiabilité
    IloInt bestj = -1;
    double best_score = -1e100;

    for (IloInt j = 0; j < vars_y.getSize(); ++j) {
        if (feas[j] != Infeasible) continue;

        double frac = y[j] - IloFloor(y[j]);
        if (frac > 0.5) frac = 1.0 - frac;

        // Récupérer les pseudocosts
        PseudocostData& pc = pseudocosts[j];
        double delta_up = pc.psi_up * (1 - frac);
        double delta_down = pc.psi_down * frac;

        // Si pseudocosts non fiables, effectuer strong-branching
        if (std::min(pc.eta_up, pc.eta_down) < eta_reliable) {
            IloCplex::BranchDirection dirs[2] = {IloCplex::BranchUp, IloCplex::BranchDown};
            double estimates[2];
            
            // Strong-branching avec limite d'itérations
            for (int d = 0; d < 2; ++d) {
                IloNumVar var = vars_y[j];
                IloNum bound = (dirs[d] == IloCplex::BranchUp) ? IloCeil(y[j]) : IloFloor(y[j]);
                estimates[d] = getObjValue() + getSlope(var, bound, dirs[d]);
            }
            
            delta_up = estimates[0] - getObjValue();
            delta_down = estimates[1] - getObjValue();
            
            // Mettre à jour les pseudocosts
            if (delta_up > 0) {
                pc.sigma_up += delta_up / (1 - frac);
                pc.eta_up++;
                pc.psi_up = pc.sigma_up / pc.eta_up;
            }
            if (delta_down > 0) {
                pc.sigma_down += delta_down / frac;
                pc.eta_down++;
                pc.psi_down = pc.sigma_down / pc.eta_down;
            }
        }

        // Calcul du score
        double score = (1 - mu) * std::min(delta_up, delta_down) + mu * std::max(delta_up, delta_down);
        if (score > best_score) {
            best_score = score;
            bestj = j;
        }
    }

    if (bestj >= 0) {
        makeBranch(vars_y[bestj], y[bestj], IloCplex::BranchUp, getObjValue());
        makeBranch(vars_y[bestj], y[bestj], IloCplex::BranchDown, getObjValue());
    }

    y.end();
    feas.end();
}
*/

ILOBRANCHCALLBACK1(MyBranch, IloNumVarArray, vars_y) {
    if (getBranchType() != BranchOnVariable) return;
    std::cout << "HI" << std::endl;

    IloEnv env = getEnv();
    IloNumArray    x(env);
    IntegerFeasibilityArray feas(env);
    IloInt n = vars_y.getSize();

    try {
        getValues(x, vars_y);
        getFeasibilities(feas, vars_y);

        const IloNum reliability = 4.0;   // Seuil de fiabilité
        const IloNum alpha       = 1.0/6; // Poids pour score
        IloInt bestVar    = -1;
        IloNum bestScore = -IloInfinity;

        for (IloInt j = 0; j < n; ++j) {
            if (feas[j] != Infeasible) continue;

            // fraction et pseudo-costs
            IloNum frac    = x[j] - IloFloor(x[j]);
            IloNum fDown   = frac;
            IloNum fUp     = 1.0 - frac;
            IloNum pcDown  = getDownPseudoCost(vars_y[j]);
            IloNum pcUp    = getUpPseudoCost(vars_y[j]);

            // si les deux pseudo-costs sont faibles, on laisse CPLEX décider
            // (on garde tout de même le score à zéro pour ne pas le privilégier)
            if (std::min(pcDown, pcUp) < reliability) {
                pcDown = pcUp = 0.0; 
            }

            // score final
            IloNum downImpact = fDown * pcDown;
            IloNum upImpact   = fUp   * pcUp;
            IloNum score = (1.0 - alpha) * std::min(downImpact, upImpact)
                          + alpha * std::max(downImpact, upImpact);

            if (score > bestScore) {
                bestScore = score;
                bestVar   = j;
            }
        }

        if (bestVar >= 0) {
            IloNum xj = x[bestVar];
            makeBranch(vars_y[bestVar], xj, IloCplex::BranchUp,   getObjValue());
            makeBranch(vars_y[bestVar], xj, IloCplex::BranchDown, getObjValue());
        }
    }
    catch (...) {
        x.end();
        feas.end();
        throw;
    }
    x.end();
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
    // Configuration du strong branching
    cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect,
        IloCplex::VariableSelect::Strong);
    cplex.setParam(IloCplex::Param::Simplex::Limits::Iterations,100);
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
    //cplex.use(MyBranch(env, y));
    cplex.use(MyBranch(env, y));
    
    
    // Résoudre le problème
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
