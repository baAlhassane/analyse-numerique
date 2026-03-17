#ifndef POUTRE_H
#define POUTRE_H
#include "TypeEF.h"
#include <vector>
#include <string>
#include <functional>

#include <vector>
#include <functional>



class Poutre {
public:
    int n;
    double L, h;
    TypeEF typePolynome;
    
    // Fonctions physiques
    std::function<double(double)> c; 
    std::function<double(double)> f;

    // Vecteurs de la matrice tridiagonale et du système
    std::vector<double> d_inf1;      // Taille n-1
    std::vector<double> d_centrale; // Taille n
    std::vector<double> d_sup1;      // Taille n-1
    std::vector<double> b;          // Taille n (Second membre)
    std::vector<double> u;          // Taille n (Solution)

    
    // Pour P2
    std::vector<double> d_inf2;
    std::vector<double> d_sup2;

    // Constructeur adapté
Poutre(int nb_elements_entree, double longueur, TypeEF p,
       std::function<double(double)> fonc_c, 
       std::function<double(double)> fonc_f) 
    : L(longueur), typePolynome(p), c(fonc_c), f(fonc_f) 
{
    // Calcul correct du nombre de noeuds totaux
    if(typePolynome == TypeEF::P1) {
        n = nb_elements_entree + 1; 
    }
    else {  
        // Pour P2, chaque élément a 2 sous-segments
        n = 2 * nb_elements_entree + 1;
    }
    
    h = L / (n - 1);

    // Initialisation rigoureuse des tailles
    d_centrale.assign(n, 0.0);
    b.assign(n, 0.0);
    u.assign(n, 0.0);

    // Diagonales 1 (voisins immédiats)
    d_inf1.assign(n-1 ,0.0); // On met n pour éviter les out_of_range
    d_sup1.assign(n-1, 0.0);

    // Diagonales 2 (saut de 1 noeud) - CRITIQUE POUR P2
    d_inf2.assign(n-2, 0.0); 
    d_sup2.assign(n-2, 0.0);
}


};

#endif