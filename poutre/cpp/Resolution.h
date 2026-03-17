#include "Poutre.h"
class Resolution {
public:
  // Méthode LU 
    static void tridiagonal(Poutre& p);

    static void  pentadiagonal(Poutre& p);

    


       // Méthode directe
    static void cholesky(Poutre& p) ;
    // Métho

     static void jacobi(Poutre& p, int max_iter, double tol);

    static void gaussSeidel(Poutre& p, int max_iter, double tol) ;

   
};


// class Resolution {
// public:
//     // --- Pour P1 et Différences Finies ---
//     static void thomasTridiag(Poutre& p);

//     // --- Pour P2 (Pentadiagonale) ---
//     static void thomasPenta(Poutre& p); 

//     // --- Méthodes Itératives (Marchent pour les deux si bien codées) ---
//     static void gaussSeidel(Poutre& p, int max_iter, double tol);
// };