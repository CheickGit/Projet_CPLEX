#include <iostream>
#include <string>


#include "modele1_lineaire.h"



int main(int argc, char * argv[]) {

	if(argc != 2){
		std::cout << "Nombre de parametre non valide\n Ex: ./prog instance_base\n";
		return 0;
	}
	
    std::string instance(argv[1]);
	 // Phase 1: Collecte des données
	MCND mcnd(instance);
    mcnd.solve(instance);
	// Phase 2: Entraînement du modèle ML
	std::cout << "Phase 2: Entraînement du modèle ML..." << std::endl;
	int train_status = system("./.venv/bin/python entraine.py"); // Chemin relatif
	if (train_status != 0) {
		std::cerr << "Erreur lors de l'entraînement du modèle ML (code de sortie: " << train_status << "). Arrêt." << std::endl;
		return false;
	}
	std::cout << "Entraînement du modèle ML terminé." << std::endl;
	// Phase 3: Résolution avec ML
    MCND mcnd_ml(instance);
    mcnd_ml.solve_with_ml(instance);
	return 0;
}
