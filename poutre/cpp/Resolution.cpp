#include"Resolution.h"
#include<iostream>


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

void Resolution::pentadiagonal(Poutre& p) {

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
    U_diag1[0]=p.d_sup1[0]/L_diag_principal[0]; //  U_(i,i+1) = [A_(i,i+1) -( element gauche de de de laligne i de L *  colonne i+1 de U )  /L(i,is)  // Ic y a pas un element gauche de de la ligne 0  donc on utilise seulement le premier terme. 
    U_diag2[0]=p.d_sup2[0]/L_diag_principal[0]; //  U_(i,i21) = [A_(i,i+2) / L_(i,i)
   
    L_diag1[1]= p.d_inf1[1];  //   

    
    L_diag_principal[1]=p.d_centrale[1]-L_diag1[0]*U_diag1[0]; // L_(i,i)= A_(i,i)- (element de gauche de deb la ligne i de L * elements  colonne i de U)

 
    //determination de U_(1,2), U_(1,3)
    U_diag1[1]=p.d_sup1[1]/L_diag_principal[1];
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
     void Resolution::cholesky(Poutre& p) {
        // Ton code de décomposition tridiagonale (vu plus haut)
        
    }

    // Métho

  void Resolution::jacobi(Poutre& p, int max_iter, double tol){

  }

   void Resolution::gaussSeidel(Poutre& p, int max_iter, double tol) {
       
    }