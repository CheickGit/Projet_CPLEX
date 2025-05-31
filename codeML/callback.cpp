#include "modele1_lineaire.h"

ILOSTLBEGIN

//----------------------------------------------------------------------------------
// Classe BranchingCallback pour contrôler le branching
class BranchingCallback : public IloCplex::BranchCallbackI {
public:
    // Constructeur
    BranchingCallback(IloEnv env, IloNumVarArray y) : IloCplex::BranchCallbackI(env), yVars(y) {}

    // Méthode principale de la callback
    void main() override {
        // Vérifier si le nœud est fractionnaire
        if (!isIntegerFeasible()) {
            IloNumArray yValues(getEnv());
            getValues(yValues, yVars); // Récupérer les valeurs des variables y

            // Trouver la variable y la plus proche de 0.5
            IloInt bestIndex = -1;
            double closestToHalf = 1.0; // Initialiser avec une valeur maximale
            for (IloInt i = 0; i < yVars.getSize(); ++i) {
                double fractionalPart = std::abs(yValues[i] - std::round(yValues[i]));
                if (fractionalPart > 1e-6 && fractionalPart < closestToHalf) {
                    closestToHalf = fractionalPart;
                    bestIndex = i;
                }
            }

            // Si une variable fractionnaire est trouvée, créer deux branches
            if (bestIndex != -1) {
                std::cout << "Branching on variable y[" << bestIndex << "] with value " 
                          << getValue(yVars[bestIndex]) << std::endl;
                makeBranch(yVars[bestIndex], 0.0, IloCplex::BranchDown, getObjValue());
                makeBranch(yVars[bestIndex], 1.0, IloCplex::BranchUp, getObjValue());
            }
            

            yValues.end(); // Libérer les ressources
        }
    }

    // Méthode pour créer une nouvelle instance de la callback
    IloCplex::CallbackI* duplicateCallback() const override {
        return new (getEnv()) BranchingCallback(getEnv(), yVars);
    }

private:
    IloNumVarArray yVars; // Variables binaires y
};

//----------------------------------------------------------------------------------

MCND::MCND(const std::string & instance) : cplex(env), model(env) {
    read_data(instance);
}

//----------------------------------------------------------------------------------

void MCND::set_parameters() {
    cplex.setParam(IloCplex::Param::ClockType, 1);
    cplex.setParam(IloCplex::Param::Threads, 1);
    //cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Param::TimeLimit, 3600.0); // Time limit in seconds
}

//----------------------------------------------------------------------------------

void MCND::create_model() {
    x = IloNumVarArray(env, ndemands * narcs, 0, IloInfinity); // Variables de flux
    y = IloNumVarArray(env, narcs, 0, 1, ILOFLOAT); // Variables continues pour relaxation continue

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
            if (i == d_k[k].O)
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

bool MCND::solve(const std::string & instance) {
    set_parameters();
    create_model();

    // Ajouter la callback pour le branching
    cplex.use(new (env) BranchingCallback(env, y));

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

        IloNumArray yValues(env);
        IloNumArray xValues(env);
        cplex.getValues(yValues, y);
        cplex.getValues(xValues, x);

        // Afficher les valeurs de y
        std::cout << "Valeurs des variables y (activation des arcs):" << std::endl;
        for (IloInt i = 0; i < yValues.getSize(); ++i) {
            std::cout << "y[" << i << "] = " << yValues[i] << std::endl;
        }

        // Afficher les valeurs de x non nulles
        std::cout << "Valeurs des variables x (flux transportés):" << std::endl;
        for (IloInt i = 0; i < xValues.getSize(); ++i) {
            if (xValues[i] > 1e-6) {
                std::cout << "x[" << i << "] = " << xValues[i] << std::endl;
            }
        }

        yValues.end();
        xValues.end();
    }

        /*IloNumArray y_(env);
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

        y_.end();
        x_.end();
    }
    */
    return true;
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
    } catch (IloException &e) {
        std::cerr << "ERROR: " << e.getMessage() << std::endl;
    } catch (...) {
        std::cerr << "Error" << std::endl;
    }
}


