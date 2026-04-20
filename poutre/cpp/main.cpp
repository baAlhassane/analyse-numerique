#include <iostream>
#include <vector>
#include <functional>
#include "Poutre.h"
#include "Methode.h"
#include "Resolution.h"

#include "TypeEF.h"

#include <iomanip> // Indispensable pour setw et setprecision
#include <chrono>


//cd ../../../home/alhassaneba/document/analyse_numerique/poutre/cpp/
 //g++ -Wall -o simulation main.cpp Poutre.cpp Methode.cpp Resolution.cpp
 // ./simulation

int main() {
    try {
        // 1. Paramètres communs
        int nb_noeuds = 10; 
        int n=10;
        double longueur = 1.0;
        auto f = [](double x) { return 1.0; };
        auto c = [](double x) { return 0.0; };

   
   // =========================================================
// SECTION P1 : TRIDIAGONAL - BI-ENCASTRÉ
// =========================================================
Poutre poutreP1(nb_noeuds, longueur, TypeEF::P1, c, f);
Methode::appliquerEF_P1(poutreP1);

// --- Blocage GAUCHE (Noeud 0) ---
poutreP1.d_centrale[0] = 1.0;
poutreP1.d_sup1[0] = 0.0;
poutreP1.b[0] = 0.0;

// --- Blocage DROIT (Noeud n-1) ---
int lastP1_simple = poutreP1.n - 1;
poutreP1.d_centrale[lastP1_simple] = 1.0;
poutreP1.d_inf1[lastP1_simple - 1] = 0.0; // On coupe la liaison avec l'avant-dernier noeud
poutreP1.b[lastP1_simple] = 0.0;

auto start_thomas = std::chrono::high_resolution_clock::now();
Resolution::tridiagonal(poutreP1);
auto end_thomas = std::chrono::high_resolution_clock::now();
// Calcul de la durée
std::chrono::duration<double, std::micro> time_thomas = end_thomas - start_thomas;
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

      // Forçage manuel des conditions aux limites dans la structure de la poutre


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
        auto start_penta = std::chrono::high_resolution_clock::now();
        Resolution::pentadiagonal(poutreP2);
        auto end_penta = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::micro> time_penta = end_penta - start_penta;
        std::cout << "--- Systeme P2 resolu (Pentadiagonal) ---" << std::endl;

        // =========================================================
        // 5. Comparaison des résultats
        // =========================================================
  // Pour l'affichage final, utilise le h de chaque poutre



// =========================================================
// SECTION CHOLESKY P1 : BI-ENCASTRÉE (0 et L)
// =========================================================
Poutre poutre_cholesky_P1(nb_noeuds, longueur, TypeEF::P1, c, f);
Methode::appliquerEF_P1(poutre_cholesky_P1);

// --- CONDITION LIMITE GAUCHE (Noeud 0) ---
poutre_cholesky_P1.d_centrale[0] = 1.0;
poutre_cholesky_P1.d_sup1[0]     = 0.0; 
poutre_cholesky_P1.d_inf1[0]     = 0.0; // Symétrie pour Cholesky
poutre_cholesky_P1.b[0]          = 0.0;

// --- CONDITION LIMITE DROITE (Noeud n-1) ---
int lastP1 = poutre_cholesky_P1.n - 1;
poutre_cholesky_P1.d_centrale[lastP1] = 1.0;
poutre_cholesky_P1.d_inf1[lastP1-1]   = 0.0;
poutre_cholesky_P1.d_sup1[lastP1-1]   = 0.0; // Symétrie pour Cholesky
poutre_cholesky_P1.b[lastP1]          = 0.0;

// Résolution
auto start_cholesky_P1 = std::chrono::high_resolution_clock::now();
Resolution::cholesky_P1(poutre_cholesky_P1);
auto end_start_cholesky_P1 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::micro> time_chol_p1 = end_start_cholesky_P1 - start_cholesky_P1;

// =========================================================
// SECTION CHOLESKY P2
// =========================================================
Poutre poutre_cholesky_P2(nb_noeuds, longueur, TypeEF::P2, c, f);
Methode::appliquerEF_P2(poutre_cholesky_P2);

// CONDITIONS AUX LIMITES (Crucial pour éviter le -nan)
// Blocage Noeud 0
poutre_cholesky_P2.d_centrale[0] = 1.0;
poutre_cholesky_P2.d_sup1[0] = 0.0;
poutre_cholesky_P2.d_sup2[0] = 0.0;
poutre_cholesky_P2.d_inf1[0] = 0.0; // Symétrie
poutre_cholesky_P2.d_inf2[0] = 0.0; // Symétrie
poutre_cholesky_P2.b[0] = 0.0;

// Blocage Noeud n-1
int lastP2 = poutre_cholesky_P2.n - 1;
poutre_cholesky_P2.d_centrale[lastP2] = 1.0;
poutre_cholesky_P2.d_inf1[lastP2-1] = 0.0;
poutre_cholesky_P2.d_inf2[lastP2-2] = 0.0;
poutre_cholesky_P2.d_sup1[lastP2-1] = 0.0; // Symétrie
poutre_cholesky_P2.d_sup2[lastP2-2] = 0.0; // Symétrie
poutre_cholesky_P2.b[lastP2] = 0.0;

auto start_cholesky_P2 = std::chrono::high_resolution_clock::now();
Resolution::cholesky_P2(poutre_cholesky_P2);
auto end_cholesky_P2 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::micro> time_chol_p2 = end_cholesky_P2 - start_cholesky_P2;

// =========================================================
// SECTION JACOBI P1 & P2
// =========================================================
int max_iter = 100000;
double tol = 1e-7;

// --- JACOBI P1 ---
Poutre poutre_jac_P1(nb_noeuds, longueur, TypeEF::P1, c, f);
Methode::appliquerEF_P1(poutre_jac_P1);
// Dirichlet bords (0 et n-1)
poutre_jac_P1.d_centrale[0] = 1.0; 
poutre_jac_P1.d_sup1[0] = 0.0;
 poutre_jac_P1.b[0] = 0.0;
int lastJ1 = poutre_jac_P1.n - 1;
poutre_jac_P1.d_centrale[lastJ1] = 1.0; 
poutre_jac_P1.d_inf1[lastJ1-1] = 0.0;
 poutre_jac_P1.b[lastJ1] = 0.0;

 auto start_jacobi_P1 = std::chrono::high_resolution_clock::now();
Resolution::jacobi_P1(poutre_jac_P1, max_iter, tol);
 auto end_jacobi_P1 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> time_jac_p1 = end_jacobi_P1 - start_jacobi_P1;


// --- JACOBI P2 ---
Poutre poutre_jac_P2(nb_noeuds, longueur, TypeEF::P2, c, f);
Methode::appliquerEF_P2(poutre_jac_P2);
// Dirichlet bords
poutre_jac_P2.d_centrale[0] = 1.0;
 poutre_jac_P2.d_sup1[0] = 0.0; 
 poutre_jac_P2.d_sup2[0] = 0.0; 
 poutre_jac_P2.b[0] = 0.0;
int lastJ2 = poutre_jac_P2.n - 1;
poutre_jac_P2.d_centrale[lastJ2] = 1.0; 
poutre_jac_P2.d_inf1[lastJ2-1] = 0.0;
 poutre_jac_P2.d_inf2[lastJ2-2] = 0.0; 
 poutre_jac_P2.b[lastJ2] = 0.0;

  auto start_jacobi_P2 = std::chrono::high_resolution_clock::now();
Resolution::jacobi_P2(poutre_jac_P2, max_iter, tol);
 auto end_jacobi_P2 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> time_jac_p2   =  end_jacobi_P2 - start_jacobi_P2 ;


// =========================================================
// SECTION GAUSS-SEIDEL P1 & P2
// =========================================================

// --- GAUSS-SEIDEL P1 ---
Poutre poutre_gs_P1(nb_noeuds, longueur, TypeEF::P1, c, f);
Methode::appliquerEF_P1(poutre_gs_P1);
poutre_gs_P1.d_centrale[0] = 1.0; 
poutre_gs_P1.d_sup1[0] = 0.0; 
poutre_gs_P1.b[0] = 0.0;
poutre_gs_P1.d_centrale[poutre_gs_P1.n-1] = 1.0; 
poutre_gs_P1.d_inf1[poutre_gs_P1.n-2] = 0.0;
 poutre_gs_P1.b[poutre_gs_P1.n-1] = 0.0;

auto start_gs = std::chrono::high_resolution_clock::now();
Resolution::gaussSeidel_P1(poutre_gs_P1, max_iter, tol);
auto end_gs = std::chrono::high_resolution_clock::now();
// Pour GS, on utilise souvent les millisecondes car c'est plus lent
std::chrono::duration<double, std::milli> time_gs_p1 = end_gs - start_gs;

// --- GAUSS-SEIDEL P2 ---
Poutre poutre_gs_P2(nb_noeuds, longueur, TypeEF::P2, c, f);
Methode::appliquerEF_P2(poutre_gs_P2);
poutre_gs_P2.d_centrale[0] = 1.0; 
poutre_gs_P2.d_sup1[0] = 0.0; 
poutre_gs_P2.d_sup2[0] = 0.0; 
poutre_gs_P2.b[0] = 0.0;

int lastGS2 = poutre_gs_P2.n - 1;
poutre_gs_P2.d_centrale[lastGS2] = 1.0; 
poutre_gs_P2.d_inf1[lastGS2-1] = 0.0;
poutre_gs_P2.d_inf2[lastGS2-2] = 0.0; 
poutre_gs_P2.b[lastGS2] = 0.0;

auto start_gs_P2 = std::chrono::high_resolution_clock::now();
Resolution::gaussSeidel_P2(poutre_gs_P2, max_iter, tol);
auto end_gs_P2 = std::chrono::high_resolution_clock::now();
std::chrono::duration<double, std::milli> time_gs_p2 = end_gs_P2 -start_gs_P2 ;












std::cout << "\nComparaison (x | u_P1 | u_P2) :" << std::endl;
// 1. On utilise le nombre de points de P1 comme référence pour la boucle
int nb_points = poutreP1.n; 

std::cout << "\n" << std::string(160, '=') << std::endl;

// En-tête
std::cout << std::left << std::setw(10) << "x_val" 
          << " | " << std::setw(12) << "u_P1_Direct" 
          << " | " << std::setw(12) << "u_P2_Direct" 
          << " | " << std::setw(12) << "u_Chol_P1" 
          << " | " << std::setw(12) << "u_Chol_P2" 
          << " | " << std::setw(12) << "u_Jac_P1" 
          << " | " << std::setw(12) << "u_Jac_P2" 
          << " | " << std::setw(12) << "u_GS_P1" 
          << " | " << std::setw(12) << "u_GS_P2" << std::endl;

std::cout << std::string(160, '-') << std::endl;

// 2. Boucle d'affichage intelligente
for (int i = 0; i < nb_points; ++i) {
    // x est calculé sur la base de P1
    double x = i * poutreP1.h;
    
    // Indice correspondant pour P2 (car maillage 2x plus fin)
    int i2 = 2 * i; 

    std::cout << std::fixed << std::setprecision(5)
              << std::left << std::setw(10) << x 
              << " | " << std::setw(12) << poutreP1.u[i] 
              << " | " << std::setw(12) << (i2 < poutreP2.n ? poutreP2.u[i2] : 0.0) 
              << " | " << std::setw(12) << poutre_cholesky_P1.u[i] 
              << " | " << std::setw(12) << (i2 < poutre_cholesky_P2.n ? poutre_cholesky_P2.u[i2] : 0.0)
              << " | " << std::setw(12) << poutre_jac_P1.u[i] 
              << " | " << std::setw(12) << (i2 < poutre_jac_P2.n ? poutre_jac_P2.u[i2] : 0.0)
              << " | " << std::setw(12) << poutre_gs_P1.u[i] 
              << " | " << std::setw(12) << (i2 < poutre_gs_P2.n ? poutre_gs_P2.u[i2] : 0.0)
              << std::endl;
}

std::cout << std::string(160, '=') << std::endl;




std::cout << " " << std::endl;







std::cout << "\n" << std::string(60, '=') << std::endl;
std::cout << "          TABLEAU DES PERFORMANCES (n=" << nb_noeuds << ")" << std::endl;
std::cout << std::string(60, '-') << std::endl;
std::cout << std::left << std::setw(25) << "Methode" << " | " << "Temps de calcul" << std::endl;
std::cout << std::string(60, '-') << std::endl;

std::cout << std::setw(25) << "tridiagonal LU  (P1)"           << " | " << time_thomas.count() << " us" << std::endl;
std::cout << std::setw(25) << "Pentadiagonal LU (P2)"    << " | " << time_penta.count() << " us" << std::endl;
std::cout << std::setw(25) << "Cholesky (P1)"         << " | " << time_chol_p1.count() << " us" << std::endl;
std::cout << std::setw(25) << "Cholesky (P2)"         << " | " << time_chol_p2.count() << " us" << std::endl;
std::cout << std::setw(25) << "Jacobi (P1)"           << " | " << time_jac_p1.count() << " ms" << std::endl;
std::cout << std::setw(25) << "Jacobi (P2)"           << " | " << time_jac_p2.count() << " ms" << std::endl;
std::cout << std::setw(25) << "Gauss-Seidel (P1)"     << " | " << time_gs_p1.count() << " ms" << std::endl;
std::cout << std::setw(25) << "Gauss-Seidel (P2)"     << " | " << time_gs_p2.count() << " ms" << std::endl;





std::cout << std::string(60, '=') << std::endl;


// std::cout << "\n" << std::string(110, '=') << std::endl;
// std::cout << std::left << std::setw(15) << "x_val" 
//           << " | " << std::setw(12) << "u_P1" 
//           << " | " << std::setw(12) << "u_P2" 
//           << " | " << std::setw(18) << "u_chol_P1" 
//           << " | " << std::setw(18) << "u_chol_P2" << std::endl;
// std::cout << std::string(110, '-') << std::endl;

// for (int i = 0; i < nb_points; ++i) {
//     std::cout << std::fixed << std::setprecision(4)
//               << std::left << std::setw(15) << i * poutreP1.h 
//               << " | " << std::setw(12) << poutreP1.u[i] 
//               << " | " << std::setw(12) << poutreP2.u[i] 
//               << " | " << std::setw(18) << poutre_cholesky_P1.u[i] 
//               << " | " << std::setw(18) << poutre_cholesky_P2.u[i] 
//               << std::endl;
// }
std::cout << std::string(110, '=') << std::endl;

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