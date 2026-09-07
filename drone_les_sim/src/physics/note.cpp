// Viscosités aux 4 interfaces de la maille (i, j)
double nu_E = 0.5 * (nu_totale[i+1][j] + nu_totale[i][j]);
double nu_O = 0.5 * (nu_totale[i][j]   + nu_totale[i-1][j]);
double nu_N = 0.5 * (nu_totale[i][j+1] + nu_totale[i][j]);
double nu_S = 0.5 * (nu_totale[i][j]   + nu_totale[i][j-1]);

// 1. Calcul du terme Vy1 (interfaces Est / Ouest)
double dudy_E = (u[i+1][j+1] - u[i+1][j-1] + u[i][j+1] - u[i][j-1]) / (4.0 * dy);
double dudy_O = (u[i][j+1]   - u[i][j-1]   + u[i-1][j+1] - u[i-1][j-1]) / (4.0 * dy);

double flux_E = nu_E * ((v[i+1][j] - v[i][j]) / dx + dudy_E);
double flux_O = nu_O * ((v[i][j]   - v[i-1][j]) / dx + dudy_O);

double Vy1 = (flux_E - flux_O) / dx;

// 2. Calcul du terme Vy2 (interfaces Nord / Sud)
double flux_N = 2.0 * nu_N * ((v[i][j+1] - v[i][j]) / dy);
double flux_S = 2.0 * nu_S * ((v[i][j]   - v[i][j-1]) / dy);

double Vy2 = (flux_N - flux_S) / dy;

// 3. Somme visqueuse totale pour v
double Vy_total = Vy1 + Vy2;





#include "Field.hpp"

/**
 * Calcule le terme d'advection pour la vitesse v au centre de la maille (i, j).
 * 
 * @param u Champ de vitesse horizontale
 * @param v Champ de vitesse verticale
 * @param i Index x de la maille
 * @param j Index y de la maille
 * @param dx Pas d'espace x
 * @param dy Pas d'espace y
 * @param use_upwind Si true, utilise Upwind (ordre 1) ; si false, utilise Centré (ordre 2)
 */
double computeAdvectionV(
    const Field2D& u, 
    const Field2D& v, 
    int i, int j, 
    double dx, double dy, 
    bool use_upwind = false) 
{
    // --- 1. VITESSES TRANSPORTANTES AUX INTERFACES ---
    double u_E = 0.5 * (u(i + 1, j) + u(i, j));
    double u_O = 0.5 * (u(i, j)     + u(i - 1, j));
    
    double v_N = 0.5 * (v(i, j + 1) + v(i, j));
    double v_S = 0.5 * (v(i, j)     + v(i, j - 1));

    // --- 2. VALEURS ADVECTÉES (v) AUX INTERFACES ---
    double v_E_adv, v_O_adv, v_N_adv, v_S_adv;

    if (use_upwind) {
        // Schéma Upwind (Décentré amont)
        v_E_adv = (u_E >= 0.0) ? v(i, j)     : v(i + 1, j);
        v_O_adv = (u_O >= 0.0) ? v(i - 1, j) : v(i, j);
        
        v_N_adv = (v_N >= 0.0) ? v(i, j)     : v(i, j + 1);
        v_S_adv = (v_S >= 0.0) ? v(i, j - 1) : v(i, j);
    } else {
        // Schéma Centré (Ordre 2 - Recommandé pour LES)
        v_E_adv = 0.5 * (v(i + 1, j) + v(i, j));
        v_O_adv = 0.5 * (v(i, j)     + v(i - 1, j));
        
        v_N_adv = 0.5 * (v(i, j + 1) + v(i, j));
        v_S_adv = 0.5 * (v(i, j)     + v(i, j - 1));
    }

    // --- 3. CALCUL DES FLUX ET DU BILAN DE MAILLE ---
    double flux_E = u_E * v_E_adv;
    double flux_O = u_O * v_O_adv;
    
    double flux_N = v_N * v_N_adv;
    double flux_S = v_S * v_S_adv;

    double adv_x = (flux_E - flux_O) / dx;
    double adv_y = (flux_N - flux_S) / dy;

    return adv_x + adv_y;
}




#ifndef SMAGORINSKY_MODEL_HPP
#define SMAGORINSKY_MODEL_HPP

#include <cmath>
#include <algorithm>
#include "TurbulenceModel.hpp" // Hérite de la classe abstraite

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

#endif // SMAGORINSKY_MODEL_HPP




void solveur_NavierStokes_X(
    int nx, int ny, double dx, double dy, double dt, double rho,
    const std::vector<std::vector<double>>& u,
    const std::vector<std::vector<double>>& v,
    const std::vector<std::vector<double>>& P,
    const std::vector<std::vector<double>>& nu_totale,
    std::vector<std::vector<double>>& u_next) 
{
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {

            // --- A. Termes d'Advection ---
            double dudx = (u[i+1][j] - u[i-1][j]) / (2.0 * dx);
            double dudy = (u[i][j+1] - u[i][j-1]) / (2.0 * dy);
            double advection = u[i][j] * dudx + v[i][j] * dudy;

            // --- B. Terme de Pression ---
            double grad_p_x = (P[i+1][j] - P[i-1][j]) / (2.0 * rho * dx);

            // --- C. Viscosités aux interfaces ---
            double nu_E = nu_totale[i+1][j] + nu_totale[i][j]; // 2 * moyenne
            double nu_O = nu_totale[i-1][j] + nu_totale[i][j];
            double nu_N = 0.5 * (nu_totale[i][j+1] + nu_totale[i][j]);
            double nu_S = 0.5 * (nu_totale[i][j-1] + nu_totale[i][j]);

            // --- D. Premier terme visqueux V_x1 ---
            double Vx1 = (nu_E * (u[i+1][j] - u[i][j]) - nu_O * (u[i][j] - u[i-1][j])) / (dx * dx);

            // --- E. Second terme visqueux V_x2 ---
            // Gradients de v aux interfaces N et S
            double dvdx_N = (v[i+1][j+1] - v[i-1][j+1] + v[i+1][j] - v[i-1][j]) / (4.0 * dx);
            double dvdx_S = (v[i+1][j] - v[i-1][j] + v[i+1][j-1] - v[i-1][j-1]) / (4.0 * dx);

            double flux_N = nu_N * ((u[i][j+1] - u[i][j]) / dy + dvdx_N);
            double flux_S = nu_S * ((u[i][j] - u[i][j-1]) / dy + dvdx_S);

            double Vx2 = (flux_N - flux_S) / dy;

            // --- F. Mise à jour temporelle (Euler) ---
            u_next[i][j] = u[i][j] + dt * (-advection - grad_p_x + Vx1 + Vx2);
        }
    }
}








// --- VARIABLES REQUISES ---
// u[nx][ny], v[nx][ny] : vitesses au présent
// u_next[nx][ny], v_next[nx][ny] : vitesses au futur
// nu_t[nx][ny] : matrice pour stocker la viscosité turbulente
// dx, dy, dt, rho, nu_physique, Cs, delta

// 1. CALCUL DE LA VISCOSITÉ TURBULENTE PAR MAILLE
for (int i = 1; i < nx - 1; ++i) {
    for (int j = 1; j < ny - 1; ++j) {
        // Gradients au centre de la maille (schéma centré)
        double dudx = (u[i+1][j] - u[i-1][j]) / (2.0 * dx);
        double dvdy = (v[i][j+1] - v[i][j-1]) / (2.0 * dy);
        double dudy = (u[i][j+1] - u[i][j-1]) / (2.0 * dy);
        double dvdx = (v[i+1][j] - v[i-1][j]) / (2.0 * dx);

        // Magnitude du tenseur de déformation S
        double S_xy = 0.5 * (dudy + dvdx);
        double S_mag = sqrt(2.0 * (dudx*dudx + dvdy*dvdy + 2.0 * S_xy*S_xy));

        // Formule de Smagorinsky
        nu_t[i][j] = (Cs * delta) * (Cs * delta) * S_mag;
    }
}

// 2. MISE À JOUR DE NAVIER-STOKES (Exemple sur l'axe X : vitesse u)
for (int i = 1; i < nx - 1; ++i) {
    for (int j = 1; j < ny - 1; ++j) {
        
        // --- Terme d'Advection (Upwind ou Centré selon ta stabilité) ---
        double advection_x = u[i][j] * (u[i+1][j] - u[i-1][j]) / (2.0 * dx) +
                             v[i][j] * (u[i][j+1] - u[i][j-1]) / (2.0 * dy);

        // --- Terme de Pression ---
        double grad_p_x = (P[i+1][j] - P[i-1][j]) / (2.0 * rho * dx);

        // --- Termes Visqueux LES (Prise en compte de la variation de nu) ---
        // Viscosités totales aux interfaces des mailles (moyennes inter-mailles)
        double nu_est   = nu_physique + 0.5 * (nu_t[i+1][j] + nu_t[i][j]);
        double nu_ouest = nu_physique + 0.5 * (nu_t[i-1][j] + nu_t[i][j]);
        double nu_nord  = nu_physique + 0.5 * (nu_t[i][j+1] + nu_t[i][j]);
        double nu_sud   = nu_physique + 0.5 * (nu_t[i][j-1] + nu_t[i][j]);

        // Dérivées des contraintes (Discrétisation de la divergence)
        double div_sigma_x = (nu_est * (u[i+1][j] - u[i][j]) - nu_ouest * (u[i][j] - u[i-1][j])) / (dx * dx) +
                             (nu_nord * (u[i][j+1] - u[i][j]) - nu_sud * (u[i][j] - u[i][j-1])) / (dy * dy);

        // --- Assemblage Temporel (Euler explicite ici pour l'exemple) ---
        u_next[i][j] = u[i][j] + dt * (-advection_x - grad_p_x + div_sigma_x);
    }
}