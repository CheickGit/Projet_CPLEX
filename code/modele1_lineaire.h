#ifndef CG_H
#define CG_H

#include <ilcplex/ilocplex.h>

#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath> 
//#include <torch/script.h> // LibTorch
#include <fstream>


struct Arc{
    int i, j;
    double capa;
    double f;
    std::vector<double> c;
    std::vector<double> b;
    inline Arc(int ii, int jj, double cap, double ff, int nk, double k):
    i(ii), j(jj), capa(cap), f(ff){
        c.resize(nk,k);
    }
    inline Arc(int nk, double k){c.resize(nk,k);}
    inline Arc(int nk){c.resize(nk,0.0);}
    inline Arc(){}
    inline ~Arc(){c.clear();b.clear();}
};


struct Demand{
    int O, D;
    double quantity;
};
// Déclaration des pseudocosts
struct Pseudocost {
    double up_score = 0.0;
    double down_score = 0.0;
    int up_count = 0;
    int down_count = 0;
};

// Déclaration anticipée de la callback
class MyBranchCallback;

class MCND {
public:

    IloEnv env;
    IloCplex cplex;
    IloModel model;
    IloNumVarArray y;
    IloNumVarArray x; // Ajout de x comme variable membre
    IloIntVarArray yInt; // Variables entières pour le branching

    std::vector<Pseudocost> pseudocosts;  // ligne Ajouter 
    
    int narcs, nnodes, ndemands;
    double B0;
    double ub_;
    std::vector<Arc> arcs;
    std::vector<Demand> d_k;
    

	// construtores
    MCND(const std::string & instance);
    MCND():cplex(env), model(env), y(env){}
    ~MCND();
    
    void set_parameters();
	void read_data(const std::string & instance);
    void create_model();
	// resolve o problema
	bool solve(const std::string & instance);
    void ecrire_solution(const std::string & instance, const IloNumArray& x_, const IloNumArray& y_, int nbarcsActifs, double coutTotal, std::vector<int> taille_chemin);
    
    

};

#endif // CG_H
