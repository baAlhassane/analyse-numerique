#include"Resolution.h"
#include<iostream>
#include <cmath>


void Resolution::tridiagonal(Poutre& p) {
    int n = p.n;
    if (n <= 0) return;

    // L et y ont besoin de n éléments, U en a besoin de n-1
    std::vector<double> L(n); // Diagonal principal de L car les L_(i, j) avec i>j sont identique aux a
    std::vector<double> U(n - 1); // Les  premeier diagonal superieur de U. Car les U_(i,i) =1 pour Crout. 
    std::vector<double> y(n);

    // --- Étape 1 : Initialisation ---
    L[0] = p.d_centrale[0];
    if (std::abs(L[0]) < 1e-15) {
            std::cout << "ERREUR: Pivot nul P1 à i=" << 0 << std::endl;
            return;
        }
    if (L[0] == 0) return; // Sécurité : éviter division par zéro
    y[0] = p.b[0] / L[0];

    // --- Étape 2 : Descente (LU Decomposition + Ly = b) ---
    for (int i = 0; i < n - 1; ++i) {
        U[i] = p.d_sup1[i] / L[i];
        L[i+1] = p.d_centrale[i+1] - p.d_inf1[i] * U[i];
        if (std::abs(L[i+1]) < 1e-15) {
            std::cout << "ERREUR: Pivot nul P1 à i=" << i+1 << std::endl;
            return;
        }
        
        if (L[i+1] == 0) return; // Sécurité
        y[i+1] = (p.b[i+1] - p.d_inf1[i] * y[i]) / L[i+1];
    }

    // --- Étape 3 : Remontée (U u = y) ---
    // On écrit directement dans p.u qui est le vecteur de résultat de la poutre 
    p.u[n-1] = y[n-1];
    for (int i = n - 2; i >= 0; i--) {
        p.u[i] = y[i] - U[i] * p.u[i+1];
    }
}

void Resolution::pentadiagonal(Poutre& p) { //P2 

    int n = p.n;
    if (n <= 0) return;
        std::vector<double> L_diag1(n-1,0), L_diag2(n-1,0),  L_diag_principal(n,0), U_diag1(n-1,0), U_diag2(n-2,0), y(n,0), b(n,0),u(n,0);

    //std::vector<double> L_diag1(n,0), L_diag2(n,0),  L_diag_principal(n,0), U_diag1(n,0), U_diag2(n,0), y(n,0), b(n,0),u(n,0);

    L_diag_principal[0]=p.d_centrale[0];
    if (L_diag_principal[0] == 0) {
    std::cout << "!!! PIVOT NUL DETECTE A L'INDICE i = " << 0 << " !!!" << std::endl;
    return;
}
   
    //determination de U_(0,1), U_(0,2)
    U_diag1[0]=p.d_sup1[0]/L_diag_principal[0]; //  U_(i,i+1) = [A_(i,i+1) -( element gauche de de de la ligne i de L *  colonne i+1 de U )  /L_(i,i)  // Ic y a pas un element gauche de de la ligne 0  donc on utilise seulement le premier terme. 
    U_diag2[0]=p.d_sup2[0]/L_diag_principal[0]; //  U_(i,i+1) = [A_(i,i+2) / L_(i,i) 
   
    L_diag1[1]= p.d_inf1[1];  //   

    
    L_diag_principal[1]=p.d_centrale[1]-L_diag1[0]*U_diag1[0]; // L_(i,i)= A_(i,i)- (element de gauche de la ligne i de L * elements  colonne i de U)

 
    //determination de U_(1,2), U_(1,3) 
    U_diag1[1]=(p.d_sup1[1]- L_diag1[0]*U_diag2[0]) /L_diag_principal[1];
    U_diag2[1]=p.d_sup2[1]/L_diag_principal[1];


    y[0] = p.b[0] / L_diag_principal[0]; // descente initiale  
    y[1]=(p.b[1] - L_diag1[0]*y[0]) /L_diag_principal[1]   ;


   

    for (int i = 2; i < n; i++)
    {
     L_diag2[i-2]=p.d_inf2[i-2]; // y'a pas d'element a sa gauche pour corriger
     L_diag1[i-1]=p.d_inf1[i-1]  - L_diag2[i-2]*U_diag1[i-2]; // car prend la correction de l'element de sa gauche * u juste au dessus
     L_diag_principal[i]= p.d_centrale[i] -(  L_diag2[i-2]*U_diag2[i-2] + L_diag1[i-1]*U_diag1[i-1]);// prend en compte de tous les elements de sa gauche

    
   
     // des qu'on connais le L_(i,i) on determeine les  U de cette ligne et commencer la descente Ly=b
    
     y[i]= (p.b[i] - L_diag1[i-1]*y[i-1] - L_diag2[i-2]*y[i-2] )/L_diag_principal[i];


     if(i<n-2) {     
         U_diag2[i]=p.d_sup2[i]/L_diag_principal[i];
        
 }

 if(i<n-1){
     U_diag1[i]=( p.d_sup1[i] - L_diag1[i-1]*U_diag2[i-1])/L_diag_principal[i];

 }
     



     
    }

   // U*u=y; Calcul de upar remontee

   p.u[n-1]= y[n-1];

   p.u[n-2]= y[n-2] - U_diag1[n-2]*p.u[n-1];


   for (int i = n-3; i >=0 ; --i)
   {
    
    p.u[i]=y[i]- (U_diag1[i]*p.u[i+1] + U_diag2[i]* p.u[i+2]);

   }
   


// AJOUTE CECI à la fin de la fonction pentadiagonal pour déboguer
std::cout << "DEBUG SOLVEUR: u[5] final = " << p.u[5] << std::endl;


}




     void Resolution::cholesky_P1(Poutre& p) {
        // Ton code de décomposition tridiagonale (vu plus haut)

        int n = p.n;
        std::vector<double> L_diag1(n-1,0),   L_diag_principal(n,0), y(n,0), b(n,0),u(n,0);
        

        L_diag_principal[0] = std::sqrt(p.d_centrale[0]);
        if (L_diag_principal[0] == 0) {
            std::cout << "ERREUR: Pivot nul P1 à i=" << 0 << std::endl;
            return;
                }

        y[0]=p.b[0] / L_diag_principal[0];

        for (int i = 1; i < n; ++i) {
            
            L_diag1[i-1] = p.d_inf1[i-1] / L_diag_principal[i-1];
            L_diag_principal[i] = std::sqrt(p.d_centrale[i] - L_diag1[i-1] * L_diag1[i-1]);
            int argument_racine =  L_diag_principal[i];

            if (argument_racine <= 0) {
            std::cout << "ERREUR : Matrice non définie positive à i=" << i << std::endl;
            return;
        }
            y[i]=(p.b[i] - L_diag1[i-1]*y[i-1] )/ L_diag_principal[i];
    
        }
          
        p.u[n-1]= y[n-1]/L_diag_principal[n-1];
        for (int i = n-2; i >= 0; --i)
        {
            p.u[i]=(y[i]- L_diag1[i]*p.u[i+1])/L_diag_principal[i]; // Car ui 

        }
        

}



void Resolution::cholesky_P2(Poutre& p) {
    int n = p.n;
    if (n <= 0) return;

    // Déclaration des vecteurs de travail
    std::vector<double> L_diag1(n - 1, 0); // Première sous-diagonale
    std::vector<double> L_diag2(n - 2, 0); // Deuxième sous-diagonale
    std::vector<double> L_diag_principal(n, 0); 
    std::vector<double> y(n, 0); // Vecteur intermédiaire pour la descente

    // --- ÉTAPE 1 : Colonne 0 (Initialisation) ---
    L_diag_principal[0] = std::sqrt(p.d_centrale[0]);
    if (n > 1) L_diag1[0] = p.d_inf1[0] / L_diag_principal[0];
    if (n > 2) L_diag2[0] = p.d_inf2[0] / L_diag_principal[0];
    y[0] = p.b[0] / L_diag_principal[0];

    // --- ÉTAPE 2 : Colonne 1 (Une seule correction possible) ---
    if (n > 1) {
        L_diag_principal[1] = std::sqrt(p.d_centrale[1] - std::pow(L_diag1[0], 2));
        if (n > 2) L_diag1[1] = (p.d_inf1[1] - L_diag2[0] * L_diag1[0]) / L_diag_principal[1];
        if (n > 3) L_diag2[1] = p.d_inf2[1] / L_diag_principal[1];
        y[1] = (p.b[1] - L_diag1[0] * y[0]) / L_diag_principal[1];
    }

    // --- ÉTAPE 3 : Colonnes 2 à n-1 (Cas général, deux corrections) ---
    for (int j = 2; j < n; j++) {
double arg = p.d_centrale[j] - std::pow(L_diag1[j-1], 2) - std::pow(L_diag2[j-2], 2);
if (arg < 1e-15) { 
    // Si c'est proche de 0 ou négatif, c'est que la matrice est mal conditionnée
    // ou les conditions aux limites sont absentes.
    throw std::runtime_error("Erreur Cholesky : Pivot non positif");
}
 

        // Pivot
           L_diag_principal[j] = std::sqrt(arg);
        //L_diag_principal[j] = std::sqrt(p.d_centrale[j] - std::pow(L_diag1[j-1], 2) - std::pow(L_diag2[j-2], 2));

        // Sous-diagonale 1
        if (j < n - 1) {
            L_diag1[j] = (p.d_inf1[j] - L_diag2[j-1] * L_diag1[j-1]) / L_diag_principal[j];
        }
        // Sous-diagonale 2
        if (j < n - 2) {
            L_diag2[j] = p.d_inf2[j] / L_diag_principal[j];
        }
        
        // Descente Ly = b
        y[j] = (p.b[j] - L_diag1[j-1] * y[j-1] - L_diag2[j-2] * y[j-2]) / L_diag_principal[j];
    }

    // --- ÉTAPE 4 : Remontée L^T * u = y ---
    // On travaille directement sur p.u pour stocker le résultat final
    p.u[n-1] = y[n-1] / L_diag_principal[n-1];
    
    if (n > 1) {
        p.u[n-2] = (y[n-2] - L_diag1[n-2] * p.u[n-1]) / L_diag_principal[n-2];
    }

    for (int i = n - 3; i >= 0; --i) {
        p.u[i] = (y[i] - L_diag1[i] * p.u[i+1] - L_diag2[i] * p.u[i+2]) / L_diag_principal[i];
    }
}

    // Métho

  void Resolution::jacobi_P1(Poutre& p, int max_iter, double tol){

  }


    void Resolution::jacobi_P2(Poutre& p, int max_iter, double tol){

  }
   void Resolution::gaussSeidel_P1(Poutre& p, int max_iter, double tol) {
       
    }


   void Resolution::gaussSeidel_P2(Poutre& p, int max_iter, double tol) {
       
    }




    

    

// void Resolution::tridiagonal(Poutre& p) {

//         // Ton code de décomposition tridiagonale (vu plus haut)
        
//       int n= p.n;
//       std::vector<double> u(n), L(n), U(n) , b(n) , y(n);
//        L[0]=p.d_centrale[0];
//        b[0]=p.b[0];
//        y[0]= b[0]/L[0];
//         for(int i=0 ;i<n-1; ++i ){
//             //determiner  U de la ligne precedente . on connais U(i,i)=1 il reste juste u(i, i+1)
//           U[i]=p.d_sup[i]/ L[i]; //u[i,i+1]= A[i,i+1]
//           // Passer au pivot suivant
//           L[i+1]= p.d_centrale[i+1]-p.d_inf[i]*U[i];

//           //descente 
//           //y[i+1]=(p.b[i]-p.d_inf[i]*y[i])/L[i+1];
//           // CORRECTION : On utilise b[i+1] pour calculer y[i+1]
//         y[i+1] = (p.b[i+1] - p.d_inf[i] * y[i]) / L[i+1];

//         }

//         // La remontee U u=y

//         p.u[n-1]=y[n-1];
//         for(int i=n-2;i>=0; i--){

//            p.u[i]=y[i] - U[i]*p.u[i+1];




//         }
        
        
//     }

       // Méthode directe