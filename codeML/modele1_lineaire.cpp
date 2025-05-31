#include "modele1_lineaire.h"
#include <set> // Ajout de l'inclusion pour std::set
#include <algorithm> 
#include <Python.h>

ILOSTLBEGIN

#define CPLEX_EPS 1e-6

// Utilitaire pour tester si un nombre est entier
bool isInteger(double value) {
    const double tol = 1e-5;
    return std::fabs(value - std::round(value)) <= tol;
}



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


ILOBRANCHCALLBACK3(RempliData, IloNumVarArray, vars_y, std::ofstream*, csvFile, IloCplex, cplex_obj){
    if (getBranchType() != BranchOnVariable) return;

    IloEnv env = getEnv();
    IloNumArray y_vals(env);
    IntegerFeasibilityArray feas(env);

    try {
        getValues(y_vals, vars_y);
        getFeasibilities(feas, vars_y);

        int status;
        CPXENVptr cpx_env = CPXopenCPLEX(&status);
        if (!cpx_env) { 
            env.warning() << "RempliDataAndBranch: Impossible d'ouvrir l'environnement CPLEX (CPXopenCPLEX)." << std::endl;
            y_vals.end(); feas.end(); return;
        }
        CPXLPptr cpx_lp = CPXcloneprob(cpx_env, cplex_obj.getImpl()->getCplexLp(), &status);
        if (!cpx_lp) {
            env.warning() << "RempliDataAndBranch: Impossible de cloner le problème (CPXcloneprob)." << std::endl;
            CPXcloseCPLEX(&cpx_env);
            y_vals.end(); feas.end(); return;
        }

        status = CPXchgprobtype(cpx_env, cpx_lp, CPXPROB_LP);
        if (status != 0) {
            env.warning() << "RempliDataAndBranch: Impossible de changer le type de problème en LP (CPXchgprobtype)." << std::endl;
            CPXfreeprob(cpx_env, &cpx_lp);
            CPXcloseCPLEX(&cpx_env);
            y_vals.end(); feas.end(); return;
        }
        std::vector<int> candidate_lp_indices;
        std::vector<IloInt> candidate_original_indices;

        for (IloInt j = 0; j < vars_y.getSize(); ++j) {
            if (feas[j] == Infeasible && !isInteger(y_vals[j])) {
                const char* var_name = vars_y[j].getName();
                if (!var_name || strlen(var_name) == 0) {
                    char namebuf[64]; sprintf(namebuf, "y_sb_idx_%d", (int)j);
                    vars_y[j].setName(namebuf); var_name = vars_y[j].getName();
                }
                int lp_idx;
                if (CPXgetcolindex(cpx_env, cpx_lp, var_name, &lp_idx) == 0) {
                    candidate_lp_indices.push_back(lp_idx);
                    candidate_original_indices.push_back(j);
                }
            }
        }

        if (candidate_lp_indices.empty()) {
            CPXfreeprob(cpx_env, &cpx_lp); CPXcloseCPLEX(&cpx_env);
            y_vals.end(); feas.end(); return; // CPLEX branchera par défaut
        }

        int num_candidates = candidate_lp_indices.size();
        std::vector<double> down_obj_values(num_candidates);
        std::vector<double> up_obj_values(num_candidates);

        CPXsetintparam(cpx_env, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
        CPXsetintparam(cpx_env, CPXPARAM_ScreenOutput, CPX_OFF);

        status = CPXstrongbranch(cpx_env, cpx_lp, candidate_lp_indices.data(), num_candidates,
                                 down_obj_values.data(), up_obj_values.data(), 10);

        IloInt best_var_for_sb_idx = -1;
        double best_actual_sb_score = -1e+100;
        double y_val_of_sb_chosen_var = 0.0;

        if (status == 0) {
            double current_obj_val_at_node = getObjValue();
            for (int i = 0; i < num_candidates; ++i) {
                IloInt original_idx = candidate_original_indices[i];
                double down_obj = (down_obj_values[i] < 1e+19) ? down_obj_values[i] : current_obj_val_at_node + 1e7;
                double up_obj = (up_obj_values[i] < 1e+19) ? up_obj_values[i] : current_obj_val_at_node + 1e7;
                double current_down_gain = down_obj - current_obj_val_at_node;
                double current_up_gain = up_obj - current_obj_val_at_node;
                double current_actual_sb_score = 0.4 * std::min(current_down_gain, current_up_gain) + 0.6 * std::max(current_down_gain, current_up_gain);

                if (current_actual_sb_score > best_actual_sb_score) {
                    best_actual_sb_score = current_actual_sb_score;
                    best_var_for_sb_idx = original_idx;
                    y_val_of_sb_chosen_var = y_vals[original_idx];
                }
            }
        }

        CPXfreeprob(cpx_env, &cpx_lp);
        CPXcloseCPLEX(&cpx_env);

        if (best_var_for_sb_idx >= 0) { // Si une meilleure variable a été trouvée par SB
            int node_depth_feature = getCurrentNodeDepth();
            if (csvFile && csvFile->is_open()) {
                // Enregistrer les features de la variable choisie et son score SB
                *csvFile << y_val_of_sb_chosen_var << ","  // y_value
                         << node_depth_feature << ","      // node_depth
                         << best_actual_sb_score << "\n";  // score_strong_branching (cible pour ML)
                csvFile->flush();
            }
            // Effectuer le branchement sur la variable choisie par strong branching
            makeBranch(vars_y[best_var_for_sb_idx], y_vals[best_var_for_sb_idx], IloCplex::BranchUp, getObjValue());
            makeBranch(vars_y[best_var_for_sb_idx], y_vals[best_var_for_sb_idx], IloCplex::BranchDown, getObjValue());
        }
        // Si aucune variable n'est choisie (par ex. SB échoue), CPLEX utilisera sa stratégie par défaut.

    } catch (const IloException& e) {
        std::cerr << "Exception Concert dans RempliDataAndBranch: " << e << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception std dans RempliDataAndBranch: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Exception inconnue dans RempliDataAndBranch." << std::endl;
    }

    y_vals.end();
    feas.end();
}


ILOBRANCHCALLBACK4(MyMLBranch, IloNumVarArray, vars_y, int, ml_max_depth, std::string, model_path_str, std::ofstream*, resultat_ml_file){
    if (getBranchType() != BranchOnVariable) return;

    int current_node_depth_val = getCurrentNodeDepth();
    if (ml_max_depth >= 0 && current_node_depth_val > ml_max_depth) return;

    IloEnv env = getEnv(); // Renommé pour éviter conflit avec cb_env plus bas si vous gardez ce nom
    IloNumArray y_values(env);
    IntegerFeasibilityArray feas(env);
    getValues(y_values, vars_y);
    getFeasibilities(feas, vars_y);

    // Initialisation des objets Python et modules
    static py::scoped_interpreter guard{}; // S'assure que Python est initialisé

    static bool python_objects_initialized_for_ml_branch = false; // Drapeau spécifique à cette callback
    static py::module_ joblib_module_ml;
    static py::module_ numpy_module_ml;
    static py::object static_ml_model_obj;  // Renommé pour éviter confusion potentielle
    static py::object static_ml_scaler_obj;

    try {
        py::gil_scoped_acquire acquire; // Acquiert le GIL

        if (!python_objects_initialized_for_ml_branch) {
            py::module_ sys_diag = py::module_::import("sys");
            // std::cout << "MyMLBranch Diag - Python Executable: " << sys_diag.attr("executable").cast<std::string>() << std::endl;
            py::list current_sys_path = sys_diag.attr("path").cast<py::list>();
            try {
                py::module_ os_path_diag = py::module_::import("os.path");
                py::module_ sysconfig_diag = py::module_::import("sysconfig");
                std::string python_version_str = sysconfig_diag.attr("get_python_version")().cast<std::string>();
                std::string venv_site_pkg_path_str = "./.venv/lib/python" + python_version_str + "/site-packages";
                std::string abs_venv_site_pkg_path = os_path_diag.attr("abspath")(venv_site_pkg_path_str).cast<std::string>();

                bool path_exists = false;
                for (py::handle p : current_sys_path) {
                    if (os_path_diag.attr("abspath")(p).cast<std::string>() == abs_venv_site_pkg_path) {
                        path_exists = true;
                        break;
                    }
                }
                if (!path_exists) {
                    current_sys_path.attr("insert")(0, abs_venv_site_pkg_path);
                    std::cout << "INFO MyMLBranch: Ajouté à sys.path (lors de l'init unique): " << abs_venv_site_pkg_path << std::endl;
                } else {
                    std::cout << "INFO MyMLBranch: Chemin venv site-packages déjà dans sys.path (lors de l'init unique)." << std::endl;
                }
            } catch (const py::error_already_set &e_path) {
                std::cerr << "MyMLBranch Erreur (init unique) modification sys.path: " << e_path.what() << std::endl;
                PyErr_Print();
            }

            // Importations de modules 
            try {
                py::module_::import("sklearn"); // Test d'import
                joblib_module_ml = py::module_::import("joblib");
                numpy_module_ml = py::module_::import("numpy");
            } catch (const py::error_already_set &e_import) {
                 std::cerr << "MyMLBranch Erreur (init unique) import modules (joblib/numpy/sklearn): " << e_import.what() << std::endl;
            }


            py::object loaded_obj_init;
            try {
                loaded_obj_init = joblib_module_ml.attr("load")(model_path_str.c_str());
            } catch (const py::error_already_set &e_load) {
                std::cerr << "MyMLBranch Erreur (init unique) chargement modèle '" << model_path_str << "': " << e_load.what() << std::endl;
                PyErr_Print();
                static_ml_model_obj = py::none(); // Important
            }

            if (!loaded_obj_init.is_none()) {
                if (py::isinstance<py::tuple>(loaded_obj_init)) {
                    py::tuple loaded_tuple = loaded_obj_init.cast<py::tuple>();
                    if (loaded_tuple.size() == 2) {
                        static_ml_model_obj = loaded_tuple[0];
                        static_ml_scaler_obj = loaded_tuple[1];
                    } else {
                        static_ml_model_obj = loaded_obj_init; // Fallback
                        static_ml_scaler_obj = py::none();
                    }
                } else { // Pas un tuple, suppose que c'est juste le modèle
                    static_ml_model_obj = loaded_obj_init;
                    static_ml_scaler_obj = py::none();
                }
            } else if (static_ml_model_obj.is_none()){ // Si loaded_obj_init était none ou si une exception a mis static_ml_model_obj à none
                 static_ml_model_obj = py::none(); // Assurer qu'il est None
            }
            python_objects_initialized_for_ml_branch = true;
            std::cout << "MyMLBranch: Initialisation unique des objets Python terminée." << std::endl;
            if(static_ml_model_obj.is_none()) std::cout << "MyMLBranch (init): Modèle ML est None." << std::endl;
            if(static_ml_scaler_obj.is_none()) std::cout << "MyMLBranch (init): Scaler ML est None." << std::endl;

        } 


        // Utiliser les objets statiques chargés : static_ml_model_obj et static_ml_scaler_obj
        if (static_ml_model_obj.is_none() || static_ml_model_obj.is(py::none())) {
            // Modèle non disponible (échec du chargement initial), laisser CPLEX brancher par défaut
            if (resultat_ml_file && resultat_ml_file->is_open()) {
                 *resultat_ml_file << current_node_depth_val << ",MODEL_NOT_LOADED_STATICALLY,-1,0.0\n";
                 resultat_ml_file->flush();
            }
            y_values.end(); feas.end(); return;
        }

        IloInt best_var_idx = -1;
        double max_predicted_score = -1e100;
        double node_depth_feature = static_cast<double>(current_node_depth_val);

        for (IloInt j = 0; j < vars_y.getSize(); ++j) {
            if (feas[j] == Infeasible && !isInteger(y_values[j])) {
                double y_value_feature = y_values[j];

                // Utiliser numpy_module_ml
                py::array temp_array = numpy_module_ml.attr("array")(py::make_tuple(py::make_tuple(y_value_feature, node_depth_feature)));
                py::array input_np = temp_array.attr("astype")("float64");
                py::object data_to_predict = input_np;

                if (!static_ml_scaler_obj.is_none() && !static_ml_scaler_obj.is(py::none())) {
                    try { data_to_predict = static_ml_scaler_obj.attr("transform")(input_np); }
                    catch (const py::error_already_set &e_scale) {
                        std::cerr << "MyMLBranch Erreur scaler.transform: " << e_scale.what() << std::endl; PyErr_Print();
                    }
                }
                
                py::object pred_obj;
                try { pred_obj = static_ml_model_obj.attr("predict")(data_to_predict); }
                catch (const py::error_already_set &e_pred) {
                    std::cerr << "MyMLBranch Erreur model.predict: " << e_pred.what() << std::endl; PyErr_Print();
                    continue;
                }
                double predicted_score = pred_obj.attr("__getitem__")(0).cast<double>();
                
                if (predicted_score > max_predicted_score) {
                    max_predicted_score = predicted_score;
                    best_var_idx = j;
                }
            }
        }

        if (best_var_idx >= 0){
            const char* chosen_var_name = vars_y[best_var_idx].getName();
            if (!chosen_var_name || strlen(chosen_var_name) == 0) chosen_var_name = "UNNAMED_VAR_ML";
            
            if (resultat_ml_file && resultat_ml_file->is_open()) {
                *resultat_ml_file << current_node_depth_val << ","
                                  << y_values[best_var_idx]<< ","
                                  << max_predicted_score << "\n";
                resultat_ml_file->flush();
            }
            
            makeBranch(vars_y[best_var_idx], y_values[best_var_idx], IloCplex::BranchUp,   getObjValue());
            makeBranch(vars_y[best_var_idx], y_values[best_var_idx], IloCplex::BranchDown, getObjValue());
        } else {
            if (resultat_ml_file && resultat_ml_file->is_open()) {
                 *resultat_ml_file << current_node_depth_val << ",NO_VAR_CHOSEN_BY_ML,-1,0.0\n";
                 resultat_ml_file->flush();
            }
        }
    } catch (const py::error_already_set &e) {
        if (resultat_ml_file && resultat_ml_file->is_open()) { *resultat_ml_file << current_node_depth_val << ",PYTHON_ERROR_ML_CB,-1,0.0\n"; resultat_ml_file->flush();}
        std::cerr << "Erreur Pybind11/Python MyMLBranch: " << e.what() << std::endl; PyErr_Print();
    } catch (const IloException& e_ilo) {
        if (resultat_ml_file && resultat_ml_file->is_open()) { *resultat_ml_file << current_node_depth_val << ",CPLEX_ILO_ERROR_ML_CB,-1,0.0\n"; resultat_ml_file->flush();}
        std::cerr << "Erreur Concert/CPLEX dans MyMLBranch: " << e_ilo << std::endl;
    } catch (const std::exception &e_std) {
        if (resultat_ml_file && resultat_ml_file->is_open()) { *resultat_ml_file << current_node_depth_val << ",CPP_STD_ERROR_ML_CB,-1,0.0\n"; resultat_ml_file->flush();}
        std::cerr << "Erreur C++ std::exception dans MyMLBranch: " << e_std.what() << std::endl;
    } catch (...) {
        if (resultat_ml_file && resultat_ml_file->is_open()) { *resultat_ml_file << current_node_depth_val << ",UNKNOWN_ERROR_ML_CB,-1,0.0\n"; resultat_ml_file->flush();}
        std::cerr << "Erreur inconnue MyMLBranch." << std::endl;
    }

    y_values.end();
    feas.end();
}
//----------------------------------------------------------------------------------

MCND::MCND(const std::string & instance) : cplex(env), model(env), y(env), x(env) { // Initialiser y et x ici aussi
    read_data(instance);
}
//----------------------------------------------------------------------------------

void MCND::set_parameters() {
    cplex.setParam(IloCplex::Param::ClockType, 1); // 1 pour temps CPU, 2 pour temps mur (wall clock)
    // cplex.setParam(IloCplex::Param::Threads, 1); // Pour la reproductibilité des tests, utile.
    // cplex.setOut(env.getNullStream()); // Pour désactiver les logs CPLEX
    cplex.setParam(IloCplex::Param::TimeLimit, 3600.0); // Limite de temps, ex: 2 heures
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
    std::cout << "Phase 1: Collecte des données pour l'instance " << instance << std::endl;

    std::ofstream csvFile("branching_data.csv");
    if (!csvFile.is_open()) {
        std::cerr << "Erreur: Impossible d'ouvrir branching_data.csv pour écriture." << std::endl;
        return false;
    }    csvFile << "y_value,node_depth,score_strong_branching\n"; // En-tête pour le fichier d'entraînement
    
    cplex.use(RempliData(env, y, &csvFile, cplex)); 

    std::cout << "Début de la résolution CPLEX pour la collecte de données..." << std::endl;
    bool solve_status = cplex.solve();
    
    csvFile.close();

    if(solve_status) {
        std::cout << "Collecte de données terminée. Statut: Succès. Objectif: " << cplex.getObjValue() << std::endl;
    } else {
        std::cout << "Collecte de données terminée. Statut: Échec ou pas de solution trouvée. CPLEX status: " << cplex.getStatus() << std::endl;
    }
    return solve_status;
}
bool MCND::solve_with_ml(const std::string & instance) {
        set_parameters();
        create_model();
        std::cout << "Phase 3: Résolution avec ML pour l'instance " << instance << std::endl;

        // Créer et ouvrir le fichier de log pour les choix du ML
        std::ofstream resultat_ml_file("resultat_ML.csv"); 
        if (!resultat_ml_file.is_open()) {
            std::cerr << "Erreur MCND::solve_with_ml: Impossible d'ouvrir resultat_ML.csv pour écriture." << std::endl;
        } else {
            // Écrire l'en-tête du fichier CSV pour les choix du ML
            resultat_ml_file << "NodeDepth,y_value,ScorePredit\n";
        }
        
        cplex.use(MyMLBranch(env, y, -1, "ml_model.pkl", &resultat_ml_file)); 
        
        std::cout << "Début de la résolution CPLEX avec guidage par le modèle ML..." << std::endl;
        bool solve_status = cplex.solve(); // Lancer la résolution

        // Fermer le fichier de log après la résolution
        if (resultat_ml_file.is_open()) {
            resultat_ml_file.close();
        }

        // Afficher les résultats de la résolution
        if (solve_status) {
            std::cout << "Résolution avec ML terminée. Statut: Succès." << std::endl;
            IloNumArray y_(env);
            IloNumArray x_(env);
            cplex.getValues(y_, y); // Récupérer les valeurs des variables y
            cplex.getValues(x_, x); // Récupérer les valeurs des variables x
    
            double val=0;
            for (int a = 0; a < narcs; ++a) {
                val += arcs[a].f * y_[a];
                for (int k = 0; k < ndemands; ++k) {
                    val += arcs[a].c[k] *d_k[k].quantity* x_[k * narcs + a];
                }
            }
            std::cout << "Coût total = " << std::setprecision(8)<<val<<std::endl;
            y_.end();
            x_.end(); 
            
        } else {
            std::cout << "Résolution avec ML: Pas de solution trouvée ou une erreur est survenue." << std::endl;
            std::cout << "  Statut de la solution CPLEX = " << cplex.getStatus() << std::endl;
        }
        return solve_status;
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