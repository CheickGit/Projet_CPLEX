#include "modele1_lineaire.h"
#include <set> // Ajout de l'inclusion pour std::set
ILOSTLBEGIN

#define CPLEX_EPS 1e-6

//----------------------------------------------------------------------------------
int callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *cbhandle) {
    if (contextid == CPX_CALLBACKCONTEXT_BRANCHING) {
        int status;
        int statind;

        status = CPXXcallbackgetrelaxationstatus(context, &statind, 0);
        if (status || (statind != CPX_STAT_OPTIMAL && statind != CPX_STAT_OPTIMAL_INFEAS)) return 0;

        int num_vars;
        CPXXcallbackgetnumvars(context, &num_vars);

        double *x = (double*) malloc(num_vars * sizeof(double));
        double objval;

        CPXXcallbackgetrelaxationpoint(context, x, 0, num_vars - 1, &objval);

        // Exemple : chercher la variable la plus fractionnaire
        int best_idx = -1;
        double best_frac = 1.0;
        for (int i = 0; i < num_vars; i++) {
            double frac = fabs(x[i] - round(x[i]));
            if (frac > 1e-6 && frac < best_frac) {
                best_frac = frac;
                best_idx = i;
            }
        }

        if (best_idx >= 0) {
            CPXCNT child1, child2;
            double lb = 0.0, ub = 1.0;
            char sense[1] = {'L'};
            int indices[1] = {best_idx};
            double coeffs[1] = {1.0};

            // y_i <= 0
            CPXXcallbackmakebranch(context, indices, coeffs, sense, &lb, 1, objval, NULL, &child1);

            // y_i >= 1
            sense[0] = 'G';
            CPXXcallbackmakebranch(context, indices, coeffs, sense, &ub, 1, objval, NULL, &child2);
        }

        free(x);
    }
    return 0;
}

// Classe de callback pour le branching personnalisé
/*
class MyBranchCallback : public IloCplex::BranchCallbackI {
    private:
        IloNumVarArray yVars;
    public:
        // Constructeur : on prend un NumVarArray
        MyBranchCallback(IloEnv env, const IloNumVarArray& y)
          : IloCplex::BranchCallbackI(env), yVars(env)
        {
            yVars = y;
        }
    
        // Doit renvoyer une copie de la callback pour chaque thread
        IloCplex::CallbackI* duplicateCallback() const override {
            return new (getEnv()) MyBranchCallback(getEnv(), yVars);
        }
    
        void main() override {
            std::cout << "On est la\n";
            // On ne branche que si la solution LP est fractionnaire
            if (!isIntegerFeasible()) {
                IloNumArray vals(getEnv());
                getValues(vals, yVars);
    
                int bestIdx = -1;
                double minDist = 1.0;
    
                // Trouver la var fractionnaire la plus proche de 0.5
                for (IloInt i = 0; i < yVars.getSize(); ++i) {
                    double v = vals[i];
                    double frac = fabs(v - std::round(v));
                    if (frac > CPLEX_EPS) {
                        double dist = fabs(v - 0.5);
                        if (dist < minDist) {
                            minDist = dist;
                            bestIdx = i;
                        }
                    }
                }
    
                if (bestIdx != -1) {
                    IloNumVar var = yVars[bestIdx];
                    double v = vals[bestIdx];
                    // Logging
                    getEnv().out()
                      << ">>> BranchCallback: " << var.getName()
                      << " = " << v
                      << " (dist=" << minDist << ")\n";
    
                    // Créer la branche var <= 0
                    makeBranch(
                      var,           // variable
                      0.0,           // nouvelle borne
                      IloCplex::BranchDown,
                      getObjValue(), // estimation d'obj
                      (NodeData*)0
                    );
                    getEnv().out() << "    BranchDown: " 
                                   << var.getName() << " <= 0\n";
    
                    // Créer la branche var >= 1
                    makeBranch(
                      var,
                      1.0,
                      IloCplex::BranchUp,
                      getObjValue(),
                      (NodeData*)0
                    );
                    getEnv().out() << "    BranchUp:   " 
                                   << var.getName() << " >= 1\n";
                }
    
                vals.end();
            }
        }
    };
    
class MyBranchCallback : public IloCplex::BranchCallbackI {
    private:
        IloNumVarArray yVars;
        
    public:
        MyBranchCallback(IloEnv env, const IloNumVarArray& y)
            : IloCplex::BranchCallbackI(env), yVars(env) 
        {
            yVars.add(y);
        }
        
        MyBranchCallback(const MyBranchCallback& other)
            : IloCplex::BranchCallbackI(other.getEnv()), yVars(other.getEnv()) 
        {
            yVars.add(other.yVars);
        }
        
        IloCplex::CallbackI* duplicateCallback() const override {
            return new (getEnv()) MyBranchCallback(*this);
        }
        
        void main() override {
            IloNumArray vals(getEnv());
            getValues(vals, yVars);
            
            int bestIndex = -1;
            double minDist = 1.0;
            
            // Logging du branchement dans le style CPLEX
            for (IloInt i = 0; i < yVars.getSize(); ++i) {
                if (!isIntegerFeasible(vals[i])) {
                    double dist = std::abs(vals[i] - 0.5);
                    if (dist < minDist) {
                        minDist = dist;
                        bestIndex = i;
                    }
                }
            }
            
            if (bestIndex != -1) {
                IloNumVar var = yVars[bestIndex];
                getEnv().out() << "Branching on y[" << bestIndex << "] = " 
                              << vals[bestIndex] << " (closest to 0.5)" << std::endl;
                
                // Création des branches avec logging CPLEX
                makeBranch(var, 0.0, IloCplex::BranchDown, getObjValue());
                getEnv().out() << "  Created branch: y[" << bestIndex << "] <= 0" << std::endl;
                
                makeBranch(var, 1.0, IloCplex::BranchUp, getObjValue());
                getEnv().out() << "  Created branch: y[" << bestIndex << "] >= 1" << std::endl;
            }else {
                abort();
            }
            vals.end();
        }
        
    private:
        bool isIntegerFeasible(double value) const {
            return (value < 1e-6 || value > 1 - 1e-6);
        }
};
    
IloCplex::Callback MyBranchCallbackFunc(IloEnv env, const IloNumVarArray& y) {
    return (IloCplex::Callback(new (env) MyBranchCallback(env, y)));
}*/

//----------------------------------------------------------------------------------

MCND::MCND(const std::string & instance) : cplex(env), model(env) {
    read_data(instance);
}

//----------------------------------------------------------------------------------

void MCND::set_parameters() {
    cplex.setParam(IloCplex::Param::ClockType, 1);
    cplex.setParam(IloCplex::Param::Threads, 1);
    // Utiliser une recherche arborescente traditionnelle pour les callbacks
    cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);
    cplex.setParam(IloCplex::Param::MIP::Strategy::MIPStart, 1);
    //cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Param::TimeLimit, 3600.0); // Time limit in seconds
}

//----------------------------------------------------------------------------------

void MCND::create_model() {
    x = IloNumVarArray(env, ndemands * narcs, 0, IloInfinity); // Variables de flux
    y = IloNumVarArray(env, narcs, 0, 1, ILOINT); // Variables entières

    /*IloExpr obj(env);
    for (int a = 0; a < narcs; ++a) {
        obj += arcs[a].f * y[a];
        for (int k = 0; k < ndemands; ++k) {
            obj += arcs[a].c[k] * x[k * narcs + a];
        }
    }
    model.add(IloMinimize(env, obj)); // Ajout de la fonction objectif
    obj.end();*/
    
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
                constraint += d_k[k].quantity; // Demande au nœud destination
            else if (i == d_k[k].O)
                constraint -= d_k[k].quantity; // Offre au nœud origine
            model.add(constraint == 0); // Conservation des flux
            constraint.end();
        }
    }

    // Contraintes de capacité globale
    for (int e = 0; e < narcs; ++e) {
        IloExpr constraint(env);
        for (int k = 0; k < ndemands; ++k) {
            constraint += x[k * narcs + e]; // Somme des flux sur l'arc
        }
        constraint -= arcs[e].capa * y[e]; // Capacité de l'arc activé
        model.add(constraint <= 0);
        constraint.end();
    }

    // Contraintes de capacité individuelle
    for (int k = 0; k < ndemands; ++k) {
        for (int e = 0; e < narcs; ++e) {
            IloExpr constraint(env);
            constraint += x[k * narcs + e]; // Flux de la demande k sur l'arc e
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
    setup_branch_callback();

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

        
        // Afficher les arcs activés
        std::cout << "Arcs activés :\n";
        for (int e = 0; e < narcs; ++e) {
            if (y_[e] > 0.5) {
                std::cout << "Arc activé : " << arcs[e].i << " -> " << arcs[e].j << std::endl;
            }
        }

        // Afficher les flots pour chaque demande
        std::cout << "\nDemandes transportées :\n";
        for (int k = 0; k < ndemands; ++k) {
            std::cout << "Demande " << k + 1 << " :\n";
            for (int e = 0; e < narcs; ++e) {
                double flow = x_[k * narcs + e];
                if (flow > 0) {
                    std::cout << "  Arc " << arcs[e].i << " -> " << arcs[e].j
                              << " : flot = " << flow << std::endl;
                }
            }
        }
        y_.end();
        x_.end();
    }
    return true;
    
}
//----------------------------------------------------------------------------------
void MCND::setup_branch_callback() {
    BranchCallback* branchCB = new BranchCallback(y);
    CPXLONG contextMask = IloCplex::Callback::Context::Id::Branching;
    cplex.use(branchCB, contextMask);
}

//----------------------------------------------------------------------------------

void
MCND::read_data(const std::string & instance) {
    
    std::ifstream file;
    file.open(instance.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open  datafile: %s\n "<<instance;
        abort();
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
