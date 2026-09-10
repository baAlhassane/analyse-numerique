#ifndef SMAGORINSKY_MODEL_HPP
#define SMAGORINSKY_MODEL_HPP

#include <cmath>
#include <algorithm>
#include "TurbulenceModel.hpp" // Hérite de la classe abstraite

namespace drone { 
class SmagorinskyModel : public TurbulenceModel {
private:
    double Cs_; // Constante de Smagorinsky (ex: 0.17)

public:
    explicit SmagorinskyModel(double Cs = 0.17) : Cs_(Cs) {}

    /**
     * Calcule la viscosité turbulente nu_t au centre de chaque maille.
     * 
     * @param u Field2D contenant la vitesse horizontale (m/s)
     * @param v Field2D contenant la vitesse verticale (m/s)
     * @param nu_t Field2D dans lequel stocker la viscosité turbulente (m^2/s)
     * @param dx Pas d'espace selon x (m)
     * @param dy Pas d'espace selon y (m)
     */
    void computeTurbulentViscosity(
        const Field2D& u, 
        const Field2D& v, 
        Field2D& nu_t, 
        double dx, 
        double dy) override 
    {
        int nx = u.nx();
        int ny = u.ny();

        // 1. Taille du filtre LES (Moyenne géométrique de la cellule)
        double delta = std::sqrt(dx * dy);
        double l_mix_sq = std::pow(Cs_ * delta, 2.0); // (Cs * delta)^2

        // 2. Parallélisation OpenMP sur le domaine intérieur
        #pragma omp parallel for collapse(2)
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                
                // --- Gradients au centre de la maille (i, j) via différences centrées ---
                double dudx = (u(i + 1, j) - u(i - 1, j)) / (2.0 * dx);
                double dvdy = (v(i, j + 1) - v(i, j - 1)) / (2.0 * dy);
                
                double dudy = (u(i, j + 1) - u(i, j - 1)) / (2.0 * dy);
                double dvdx = (v(i + 1, j) - v(i - 1, j)) / (2.0 * dx);

                // --- Composantes du tenseur des taux de déformation S_ij ---
                double S_xx = dudx;
                double S_yy = dvdy;
                double S_xy = 0.5 * (dudy + dvdx);

                // --- Norme du tenseur |S| = sqrt(2 * S_ij * S_ij) ---
                double S_sq = 2.0 * (S_xx * S_xx + S_yy * S_yy + 2.0 * S_xy * S_xy);
                double S_norm = std::sqrt(std::max(0.0, S_sq));

                // --- Calcul de nu_t ---
                nu_t(i, j) = l_mix_sq * S_norm;
            }
        }

        // 3. Traitement des conditions aux limites pour nu_t (Ghost cells)
        applyBoundaryConditions(nu_t);
    }

private:
    void applyBoundaryConditions(Field2D& nu_t) {
        int nx = nu_t.nx();
        int ny = nu_t.ny();

        // Extrapolation de Neumann homogène aux parois (nu_t_frontiere = nu_t_interieur)
        for (int i = 0; i < nx; ++i) {
            nu_t(i, 0) = nu_t(i, 1);           // Paroi Sud
            nu_t(i, ny - 1) = nu_t(i, ny - 2); // Paroi Nord
        }
        for (int j = 0; j < ny; ++j) {
            nu_t(0, j) = nu_t(1, j);           // Paroi Ouest
            nu_t(nx - 1, j) = nu_t(nx - 2, j); // Paroi Est
        }
    }
};

}

#endif // SMAGORINSKY_MODEL_HPP