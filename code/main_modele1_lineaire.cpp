#include <iostream>
#include <string>


#include "modele1_lineaire.h"



int main(int argc, char * argv[]) {

	if(argc != 2){
		std::cout << "Nombre de parametre non valide\n Ex: ./prog instance_base\n";
		return 0;
	}
	
    std::string instance(argv[1]);
	MCND mcnd(instance);
    mcnd.solve(instance);
	return 0;
}
