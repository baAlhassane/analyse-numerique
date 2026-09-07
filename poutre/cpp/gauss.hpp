#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

class Gauss {
public: 
    void gausse_methode(std::vector<std::vector<double>>& matrice) {
        size_t m = matrice.size(); 
        size_t n = matrice[0].size();

        // Une seule boucle principale pour chaque colonne i (le pivot actuel)
        for (size_t i = 0; i < m - 1; i++) {

            // --- ÉTAPE 1 : RECHERCHE DU PIVOT PARTIEL ---
            size_t ligne_max = i;
            for (size_t k = i + 1; k < m; k++) {
                if (std::abs(matrice[k][i]) > std::abs(matrice[ligne_max][i])) {
                    ligne_max = k;
                }
            }
            
            // --- ÉTAPE 2 : ÉCHANGE DES LIGNES ---
            if (ligne_max != i) {
                std::swap(matrice[i], matrice[ligne_max]);
            }

            // --- ÉTAPE 3 : TEST DU FAUX ZÉRO ---
            double pivot = matrice[i][i];
            if (std::abs(pivot) < 1e-9) {
                // Si toute la colonne sous le pivot est nulle, on passe à la colonne suivante
                continue; 
            }
            
            // --- ÉTAPE 4 : ÉLIMINATION DES LIGNES EN DESSOUS ---
            for (size_t k = i + 1; k < m; k++) {
                double multiplicateur = matrice[k][i] / pivot; 

                // On soustrait la ligne pivot de la ligne k
                for (size_t j = i; j < n; j++) {
                    matrice[k][j] = matrice[k][j] - multiplicateur * matrice[i][j]; 
                }
            }
        }
    }

    void afficher_matrice(const std::vector<std::vector<double>>& matrice) {
        std::cout << " %%%%%%%%%%%%%%%%%%% matrice de gauss %%%%%%%%%%%%%\n";
        // Configuration de l'affichage pour voir les nombres à virgule proprement
        std::cout << std::fixed << std::setprecision(5); 
        for (const auto& ligne : matrice) {
            for (double val : ligne) {
                std::cout << std::setw(10) << val << " ";
            }
            std::cout << "\n";
        }
        std::cout << "==================================================\n";
    }
};
 