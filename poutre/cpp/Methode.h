#include "Poutre.h"
class Methode {
public:
    // Remplit la matrice via le schéma des différences finies

    static void appliquerDF(Poutre& p) ;

       
    
    // Remplit la matrice via l'approche éléments finis
    static void appliquerEF_P1(Poutre& p) ;
        // Logique d'assemblage par éléments (matrices élémentaires)

      static void appliquerEF_P2(Poutre& p);
    
};