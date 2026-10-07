#include "solvers/NavStokesSolver.hpp"

namespace drone {

void NavStokesSolver::computePredictorVelocity() {
    const double dx = mesh_.dx();
    const double dy = mesh_.dy();
    const int nx = mesh_.nx();
    const int ny = mesh_.ny();

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {

            // Viscosités effectives aux 4 faces
            double nu_E = 0.5 * (nu_t_(i+1, j) + nu_t_(i, j)) + nu_mol_;
            double nu_W = 0.5 * (nu_t_(i, j) + nu_t_(i-1, j)) + nu_mol_;
            double nu_N = 0.5 * (nu_t_(i, j+1) + nu_t_(i, j)) + nu_mol_;
            double nu_S = 0.5 * (nu_t_(i, j-1) + nu_t_(i, j)) + nu_mol_; // AJOUTÉ : nu_S

            // =================================================================
            // 1. COMPOSANTE X (u_star)
            // =================================================================

            // A. Advection de u
            double u_E = 0.5 * (u_(i, j) + u_(i+1, j));
            double u_W = 0.5 * (u_(i-1, j) + u_(i, j));
            double adv_EW = (u_E * u_E - u_W * u_W) / dx;

            double u_N = 0.5 * (u_(i, j) + u_(i, j+1));
            double v_N = 0.5 * (v_(i, j) + v_(i, j+1));
            double u_S = 0.5 * (u_(i, j-1) + u_(i, j));
            double v_S = 0.5 * (v_(i, j-1) + v_(i, j));
            double adv_NS = (u_N * v_N - u_S * v_S) / dy;

            double adv_x = adv_EW + adv_NS;
            // B. Diffusion de u
// CORRECTION : Facteur 2.0 sur la dérivée normale (contrainte normale 2*nu*du/dx)
double diff_x1 = 2.0 * (nu_E * (u_(i+1, j) - u_(i, j)) - nu_W * (u_(i, j) - u_(i-1, j))) / (dx * dx); 
double diff_x2 = (nu_N * (u_(i, j+1) - u_(i, j)) - nu_S * (u_(i, j) - u_(i, j-1))) / (dy * dy); 
 
// Termes croisés de cisaillement (d/dy (nu * dv/dx))
double dvdx_N = (v_(i+1, j+1) - v_(i-1, j+1) + v_(i+1, j) - v_(i-1, j)) / (4.0 * dx); 
double dvdx_S = (v_(i+1, j) - v_(i-1, j) + v_(i+1, j-1) - v_(i-1, j-1)) / (4.0 * dx); 
double diff_x_cross = (nu_N * dvdx_N - nu_S * dvdx_S) / dy; 

double diff_x = diff_x1 + diff_x2 + diff_x_cross;

            // // B. Diffusion de u
            // double diff_x1 = (nu_E * (u_(i+1, j) - u_(i, j)) - nu_W * (u_(i, j) - u_(i-1, j))) / (dx * dx);
            // double diff_x2 = (nu_N * (u_(i, j+1) - u_(i, j)) - nu_S * (u_(i, j) - u_(i, j-1))) / (dy * dy);
            
            // // Termes croisés de cisaillement
            // double dvdx_N = (v_(i+1, j+1) - v_(i-1, j+1) + v_(i+1, j) - v_(i-1, j)) / (4.0 * dx);
            // double dvdx_S = (v_(i+1, j) - v_(i-1, j) + v_(i+1, j-1) - v_(i-1, j-1)) / (4.0 * dx);
            // double diff_x_cross = (nu_N * dvdx_N - nu_S * dvdx_S) / dy;

            //double diff_x = diff_x1 + diff_x2 + diff_x_cross;

            // C. Mise à jour u_star
            u_star_(i, j) = u_(i, j) + dt_ * (-adv_x + diff_x);

            // =================================================================
            // 2. COMPOSANTE Y (v_star)
            // =================================================================

            // A. Advection de v
            double u_E_v = 0.5 * (u_(i, j) + u_(i+1, j));
            double v_E_v = 0.5 * (v_(i, j) + v_(i+1, j));
            double u_W_v = 0.5 * (u_(i-1, j) + u_(i, j));
            double v_W_v = 0.5 * (v_(i-1, j) + v_(i, j));
            double adv_y_EW = (u_E_v * v_E_v - u_W_v * v_W_v) / dx;

            double v_N_v = 0.5 * (v_(i, j) + v_(i, j+1));
            double v_S_v = 0.5 * (v_(i, j-1) + v_(i, j));
            double adv_y_NS = (v_N_v * v_N_v - v_S_v * v_S_v) / dy;

            double adv_y = adv_y_EW + adv_y_NS;

            // B. Diffusion de v
double diff_y1 = (nu_E * (v_(i+1, j) - v_(i, j)) - nu_W * (v_(i, j) - v_(i-1, j))) / (dx * dx); 
// CORRECTION : Facteur 2.0 sur la dérivée normale (contrainte normale 2*nu*dv/dy)
double diff_y2 = 2.0 * (nu_N * (v_(i, j+1) - v_(i, j)) - nu_S * (v_(i, j) - v_(i, j-1))) / (dy * dy); 

// Termes croisés de cisaillement (d/dx (nu * du/dy))
double dudy_E = (u_(i+1, j+1) - u_(i+1, j-1) + u_(i, j+1) - u_(i, j-1)) / (4.0 * dy); 
double dudy_W = (u_(i, j+1) - u_(i, j-1) + u_(i-1, j+1) - u_(i-1, j-1)) / (4.0 * dy); 
double diff_y_cross = (nu_E * dudy_E - nu_W * dudy_W) / dx; 

double diff_y = diff_y1 + diff_y2 + diff_y_cross;

            // // B. Diffusion de v
            // double diff_y1 = (nu_E * (v_(i+1, j) - v_(i, j)) - nu_W * (v_(i, j) - v_(i-1, j))) / (dx * dx);
            // double diff_y2 = (nu_N * (v_(i, j+1) - v_(i, j)) - nu_S * (v_(i, j) - v_(i, j-1))) / (dy * dy);

            // // Termes croisés de cisaillement
            // double dudy_E = (u_(i+1, j+1) - u_(i+1, j-1) + u_(i, j+1) - u_(i, j-1)) / (4.0 * dy);
            // double dudy_W = (u_(i, j+1) - u_(i, j-1) + u_(i-1, j+1) - u_(i-1, j-1)) / (4.0 * dy);
            // double diff_y_cross = (nu_E * dudy_E - nu_W * dudy_W) / dx;

            // double diff_y = diff_y1 + diff_y2 + diff_y_cross;

            // C. Mise à jour v_star
            v_star_(i, j) = v_(i, j) + dt_ * (-adv_y + diff_y);
        }
    }
    
    // Application des conditions aux limites
    // Application des conditions aux limites correctes
u_star_.fillBoundaryGhostCells(BCType::VelocityU, 1.0); // U avec couvercle mobile u = 1.0
v_star_.fillBoundaryGhostCells(BCType::VelocityV, 0.0); // V avec couvercle fixe v = 0.0
    // u_star_.fillBoundaryGhostCells();
    // v_star_.fillBoundaryGhostCells();
}


void NavStokesSolver::step() {
    // 1. Calcul de la viscosité turbulente (LES)
    if (les_model_) {
        les_model_->computeTurbulentViscosity(u_, v_, nu_t_, mesh_.dx(), mesh_.dy());
    }

    // 2. Étape de prédiction : calcul de u_star et v_star
    computePredictorVelocity();

    // 3 & 4. Projection complète (Poisson + Correction de vitesse + Ghost cells)
    projection_method_.project(u_star_, v_star_, u_, v_, p_, mesh_, poisson_solver_);
}


void NavStokesSolver::writeVTK(int step) {
    // Avant : "output_" + std::to_string(step) + ".vtk"
    // ✅ Après : ajouter "output/" au début du chemin
    std::string filename = "output/output_" + std::to_string(step) + ".vtk";
    
    vtk_writer_.write(filename, mesh_, u_, v_, p_);
}



} // namespace drone