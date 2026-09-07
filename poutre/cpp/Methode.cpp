    
    #include"Methode.h"
    #include<functional>
    #include<iostream>


    void Methode::appliquerEF_P1(Poutre& p) {
    double h = p.h;
    int n = p.n;

    // On récupère la fonction c de la poutre
    // (Assure-toi que p.c est bien de type std::function dans Poutre.hpp)
    auto& c_func = p.c; 
    auto& f_func = p.f; 

    double resistance_flexion_intra = 1.0 / h;
    double resistance_flexion_extra = -1.0 / h;
   

    for (int e = 0; e < n - 1; ++e) {
        // 1. On calcule la valeur de c au milieu de l'élément (x_e + h/2)
        double x_milieu = (e + 0.5) * h;
        double ce = c_func(x_milieu); 

        
        ce=0;
        ///double fe = f_func(x_milieu); // Utilisation de f(x)s
        // 2. On utilise 'ce' (la valeur numérique) pour les calculsls
        double reaction_elastique_intra = (ce * h) / 3.0;
        double reaction_elastique_extra = (ce * h) / 6.0; // Attention au signe selon ta théorie

        // 3. Assemblage (on utilise 'e' comme indice)
        p.d_centrale[e]   += resistance_flexion_intra + reaction_elastique_intra;
        p.d_centrale[e+1] += resistance_flexion_intra + reaction_elastique_intra;
        
        p.d_inf1[e] += resistance_flexion_extra + reaction_elastique_extra; 
        p.d_sup1[e] += resistance_flexion_extra + reaction_elastique_extra; 

         //double force_repartie=(fe * h) / 2.0;
    double fe = f_func(x_milieu); 
    std::cout << "DEBUG ASSEMBLAGE: element " << e << " f=" << fe << " h=" << h << std::endl;
    double force_nodale = (fe * h) / 2.0;
    p.b[e]   += force_nodale;
    p.b[e+1] += force_nodale;
    }
}

void Methode::appliquerEF_P2(Poutre& p){
      
    // h_petit est la distance entre deux noeuds (L/(n-1))
    // H_element est la longueur d'un élément P2 (2 * h_petit)
    double h = p.h;
    double h_element = 2.0 * h;
    
    // Nombre d'éléments réels
    int n_elements = (p.n - 1) / 2;
    // On récupère la fonction c de la poutre
    // (Assure-toi que p.c est bien de type std::function dans Poutre.hpp)
    auto& c_func = p.c; 
    auto& f_func = p.f; 



    // terme u''
    double resistance_flexion_intra_noeud_extrem = 7*1.0 / (3*h_element);
    double resistance_flexion_intra_noeud_milieu = 16*1.0 / (3*h_element);
    double resistance_flexion_extra_1 = -8*1.0 / (3*h_element);
    double resistance_flexion_extra_2 = 1.0 / (3*h_element);

 

    for( int e=0; e<n_elements; ++e){


      double x_milieu = (2 * e + 1) * h;
      double ce = c_func(x_milieu); 

    double reaction_elastique_intra_noeud_extrem = 4.0*(ce * h_element) /  30.0;
    double reaction_elastique_intra_noeud_milieu = 16.0*(ce * h_element) /  30.0;
    double reaction_elastique_extra_1 = 2.0*(ce * h_element) / 30.0; // Attention au signe selon ta théorie
    double reaction_elastique_extra_2 = -1.0*(ce * h_element) / 30.0; 

        // int i0 = 2 * e;     // Gauche
        // int i1 = 2 * e + 1; // Milieu
        // int i2 = 2 * e + 2; // Droite

        int i = 2 * e;  

    // diagonal a 3 elements
    p.d_centrale[i]+= resistance_flexion_intra_noeud_extrem+reaction_elastique_intra_noeud_extrem ;
    p.d_centrale[i+1]+= resistance_flexion_intra_noeud_milieu + reaction_elastique_intra_noeud_milieu;
    p.d_centrale[i+2]+= resistance_flexion_intra_noeud_extrem + reaction_elastique_intra_noeud_extrem;


    // premeire  diagonal_ininerieur 1 a 2 elements  p.d_inf1 
    p.d_inf1[i]+= resistance_flexion_extra_1 +reaction_elastique_extra_1;
    p.d_inf1[i+1]+= resistance_flexion_extra_1 +reaction_elastique_extra_1;

     // premier  diagonal_ superieur a 2 element  p.d_sup 
    p.d_sup1[i]+= resistance_flexion_extra_1 +reaction_elastique_extra_1;
    p.d_sup1[i+1]+= resistance_flexion_extra_1 +reaction_elastique_extra_1;


    //   deuxieme  diagonal_ superieur a 2 element  p.d_sup 
    p.d_inf2[i]+= resistance_flexion_extra_2 +reaction_elastique_extra_2;
    p.d_sup2[i]+= resistance_flexion_extra_2 +reaction_elastique_extra_2;



// On suppose que h est la longueur de l'élément complet (contenant 3 noeuds)
    double fe = f_func(x_milieu); 
    if (e == 0) {
        std::cout << "DEBUG P2: x_milieu=" << x_milieu << " f=" << fe << " h_elem=" << h_element << std::endl;
    }

   double common_factor = (f_func(x_milieu) * h_element) / 6.0;

    p.b[i]     += 1.0 * common_factor; // Sommet gauche
    p.b[i + 1] += 4.0 * common_factor; // Milieu
    p.b[i+ 2] += 1.0 * common_factor; // Sommet droit


    }

    }

  


void Methode::appliquerDF(Poutre& p) {
    double h = p.h;
    int n = p.n;
    auto& c_func = p.c;
    auto& f_func = p.f;

    // Le pas au carré pour la dérivée seconde
    double h2 = h * h;

    for (int i = 0; i < n; ++i) {
        double x_i = i * h;
        double ci = c_func(x_i);
        double fi = f_func(x_i);

        // 1. Terme de la diagonale centrale
        // Provient de : -(-2/h^2) + ci = 2/h^2 + ci
        p.d_centrale[i] = (2.0 / h2) + ci;

        // 2. Termes des diagonales adjacentes (u_{i-1} et u_{i+1})
        // Provient de : -(1/h^2)
        if (i > 0) {
            p.d_inf1[i-1] = -1.0 / h2;
        }
        if (i < n - 1) {
            p.d_sup1[i] = -1.0 / h2;
        }

        // 3. Second membre
        p.b[i] = fi;
    }
}





   //determine la matrice A par ELF 
// void Methode::appliquerEF(Poutre& p) {
//     double h=p.h;
//     int n=p.n;
//     int L=p.L;
//     double c=p.c;
//     std::function<double(double)> c; 

//     // la valeutr de la rigidité(u'') pour (i=j) et (u !=j )  
//      double resistance_flexion_intra=1/h;
//      double resistance_flexion_extra=-1/h;

//       // la valeutr du ressort elastique( c(x)*u) pour (i=j) et (u != j )  
//      double reaction_elastique_intra=(c*h)/3;
//      double reaction_elastique_extra=-(c*h)/6;

//      for (int e = 0; e < n - 1; ++e) {
//     double ce = c((e + 0.5) * h) ;
    

     
//     p.d_centrale[i]+= resistance_flexion_intra + reaction_elastique_intra;
//     p.d_centrale[i+1]+= resistance_flexion_intra + reaction_elastique_intra;
//     p.d_inf[e]=resistance_flexion_extra +  reaction_elastique_extra; 
//     p.d_sup[e]= resistance_flexion_extra +  reaction_elastique_extra; 




// }

       
//     }


