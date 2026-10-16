#ifndef VISCOSITY_DEBUG_HPP
#define VISCOSITY_DEBUG_HPP

#include <iostream>
#include <string>
#include <vector>

/**
 * @brief Affiche les contributions visqueuses (diffusion normale x, y, croisée et totale)
 *        aux 8 positions clés du domaine (maillage nx * ny).
 * 
 * Note : Cette fonction suppose que les tableaux/objets u_, v_ et nu_t_ 
 * sont accessibles dans le scope où cette fonction est appelée.
 */
inline void print_eight_key_positions(int nx, int ny, double dx, double dy, double nu_mol) {
    struct KeyPos {
        std::string name;
        int i;
        int j;
    };

    // Définition des 8 positions clés (indexation 1-based)
    std::vector<KeyPos> positions = {
        {"Sud-Ouest",    1,      1},
        {"Sud-Est",      nx,     1},
        {"Milieu Sud",   nx / 2, 1},
        {"Milieu Ouest", 1,      ny / 2},
        {"Milieu Est",   nx,     ny / 2},
        {"Nord-Ouest",   1,      ny},
        {"Milieu Nord",  nx / 2, ny},
        {"Nord-Est",     nx,     ny}
    };

    std::cout << "\n==================================================================" << std::endl;
    std::cout << "         BILAN VISQUEUX AUX 8 POSITIONS CLEFS (Smagorinsky)       " << std::endl;
    std::cout << "==================================================================" << std::endl;

    for (const auto& pos : positions) {
        int i = pos.i;
        int j = pos.j;

        // Calcul des viscosités effectives sur les faces
        double nu_E = 0.5 * (nu_t_(i+1, j) + nu_t_(i, j)) + nu_mol;
        double nu_W = 0.5 * (nu_t_(i-1, j) + nu_t_(i, j)) + nu_mol;
        double nu_N = 0.5 * (nu_t_(i, j+1) + nu_t_(i, j)) + nu_mol;
        double nu_S = 0.5 * (nu_t_(i, j-1) + nu_t_(i, j)) + nu_mol;

        // Calcul des termes de diffusion
        double diff_x1 = 2.0 * (nu_E * (u_(i+1, j) - u_(i, j)) - nu_W * (u_(i, j) - u_(i-1, j))) / (dx * dx); 
        double diff_x2 = (nu_N * (u_(i, j+1) - u_(i, j)) - nu_S * (u_(i, j) - u_(i, j-1))) / (dy * dy); 
        
        double dvdx_N = (v_(i+1, j+1) - v_(i-1, j+1) + v_(i+1, j) - v_(i-1, j)) / (4.0 * dx); 
        double dvdx_S = (v_(i+1, j) - v_(i-1, j) + v_(i+1, j-1) - v_(i-1, j-1)) / (4.0 * dx); 
        double diff_x_cross = (nu_N * dvdx_N - nu_S * dvdx_S) / dy; 

        double diff_x = diff_x1 + diff_x2 + diff_x_cross;

        // Affichage formaté console
        std::cout << "[" << pos.name << "] (i=" << i << ", j=" << j << ")\n"
                  << "  -> nu_t(P) = " << nu_t_(i, j) << "\n"
                  << "  -> diff_x1 = " << diff_x1 << "\n"
                  << "  -> diff_x2 = " << diff_x2 << "\n"
                  << "  -> cross   = " << diff_x_cross << "\n"
                  << "  -> Total   = " << diff_x << "\n"
                  << "------------------------------------------------------------------" << std::endl;
    }
}

#endif // VISCOSITY_DEBUG_HPP