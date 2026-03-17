#include <iostream>
#include <vector>
#include <functional>
#include "Poutre.h"
#include "Methode.h"
#include "Resolution.h"

#include "TypeEF.h"


//cd ../../../home/alhassaneba/document/analyse_numerique/poutre/cpp/

int main() {
    try {
        // 1. Paramètres communs
        int nb_noeuds = 10; 
        int n=10;
        double longueur = 1.0;
        auto f = [](double x) { return 1.0; };
        auto c = [](double x) { return 0.0; };

        // =========================================================
        // SECTION P1 : TRIDIAGONAL
        // =========================================================
        Poutre poutreP1(nb_noeuds, longueur, TypeEF::P1,  c, f);
        
        Methode::appliquerEF_P1(poutreP1);

        // Si c'est 0, on force pour voir si le solveur réagit
   
        // Condition aux limites : Bloqué à gauche (u[0] = 0)
        poutreP1.d_centrale[0] = 1.0;
        poutreP1.d_sup1[0] = 0.0; // Dans ton P1, d_sup devient d_sup1
        poutreP1.b[0] = 0.0;

        Resolution::tridiagonal(poutreP1);

        // =========================================================
        // SECTION P2 : PENTADIAGONAL (Ton nouveau code)
        // =========================================================
        Poutre poutreP2(nb_noeuds, longueur, TypeEF::P2, c, f);
        
        // 1. Assemblage de la matrice pentadiagonale
        Methode::appliquerEF_P2(poutreP2); 
        std::cout << "\n--- Matrice P2 assemblee (Pentadiagonale) ---" << std::endl;

        // 2. Conditions aux limites (Simple encastrement au noeud 0)
        // Pour P2, on doit annuler les DEUX sur-diagonales partant de la ligne 0
      // Dans le main, après appliquerEF_P2 :

// Blocage Noeud 0 (Gauche)
poutreP2.d_centrale[0] = 1.0;
poutreP2.d_sup1[0] = 0.0;
poutreP2.d_sup2[0] = 0.0;
poutreP2.b[0] = 0.0;

// Blocage Noeud n-1 (Droite) - INDISPENSABLE
int last = poutreP2.n - 1;
poutreP2.d_centrale[last] = 1.0;
poutreP2.d_inf1[last-1] = 0.0;
poutreP2.d_inf2[last-2] = 0.0;
poutreP2.b[last] = 0.0;
        // On assure aussi que les sous-diagonales n'impactent pas le noeud 0
        // (Géré par la sécurité dans ton solveur pentadiagonal)

        // VERIFICATION AVANT RESOLUTION
std::cout << "DEBUG: Diagonale milieu = " << poutreP2.d_centrale[n/2] << std::endl;
std::cout << "DEBUG: Force milieu = " << poutreP2.b[n/2] << std::endl;

        // 3. Résolution avec ton nouveau solveur
        Resolution::pentadiagonal(poutreP2);
        std::cout << "--- Systeme P2 resolu (Pentadiagonal) ---" << std::endl;

        // =========================================================
        // 5. Comparaison des résultats
        // =========================================================
  // Pour l'affichage final, utilise le h de chaque poutre
std::cout << "\nComparaison (x | u_P1 | u_P2) :" << std::endl;
int nb_points = std::min(poutreP1.n, poutreP2.n); // Pour éviter de sortir du tableau

for (int i = 0; i < nb_points; ++i) {
    std::cout << "x_P1=" << i * poutreP1.h << " | u_P1: " << poutreP1.u[i] 
              << " \t x_P2=" << i * poutreP2.h << " | u_P2: " << poutreP2.u[i] << std::endl;
}

    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << std::endl;
    }
    return 0;
}












/// main avec TypeEF::P1
/**   
int main() {
    try {
        // 1. Définition des paramètres
        int nb_noeuds = 10; // 10 segments
        double longueur = 1.0;
        
         TypeEF typePlolynome=TypeEF::P1;
    
        
        // Fonctions f(x) = 1 et c(x) = 0
        auto f = [](double x) { return 1.0; };
        auto c = [](double x) { return 0.0; };

        // 2. Initialisation de la Poutre
        Poutre maPoutre(nb_noeuds, longueur, typePlolynome, f, c);
        std::cout << "--- Poutre creee avec h = " << maPoutre.h << " ---" << std::endl;

        // 3. Construction de la matrice (Physique)
        Methode::appliquerEF_P1(maPoutre);
        std::cout << "--- Matrice assemblee ---" << std::endl;

        // --- IMPORTANT : Appliquer une condition aux limites simple ---
        // On bloque le noeud 0 (u[0] = 0)
        maPoutre.d_centrale[0] = 1.0;
        maPoutre.d_sup[0] = 0.0;
        maPoutre.b[0] = 0.0;
        // On doit aussi mettre p.d_inf[0] à 0 pour la symétrie
        maPoutre.d_inf[0] = 0.0; 

        // 4. Résolution (Mathématiques)
        Resolution::tridiagonal(maPoutre);
        std::cout << "--- Systeme resolu ---" << std::endl;

        // 5. Affichage des résultats
        std::cout << "\nResultats (Deplacements u) :" << std::endl;
        for (int i = 0; i < maPoutre.n; ++i) {
            std::cout << "Noeud " << i << " (x=" << i*maPoutre.h << ") : " << maPoutre.u[i] << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Erreur : " << e.what() << std::endl;
    }




    return 0;
}



***/