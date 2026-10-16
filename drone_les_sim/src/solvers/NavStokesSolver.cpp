#include "solvers/NavStokesSolver.hpp"

namespace drone {

void NavStokesSolver::computePredictorVelocity() {
    const double dx = mesh_.dx();
    const double dy = mesh_.dy();
    const int nx = mesh_.nx();
    const int ny = mesh_.ny();


      if (les_model_) {
        les_model_->computeTurbulentViscosity(
            u_,
            v_,
            nu_t_,
            dx,
            dy
        );
    }

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




// std::cout
//     << "i=" << i << ", j=" << j
//     << " nu_t_P=" << nu_t_(i, j)
//     << " nu_E=" << nu_E
//     << " nu_W=" << nu_W
//     << " nu_N=" << nu_N
//     << " nu_S=" << nu_S
//     << " diff_x1=" << diff_x1
//     << " diff_x2=" << diff_x2
//     << " diff_x_cross=" << diff_x_cross
//     << " diff_x=" << diff_x
//     << '\n';

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
    

    //print_eight_key_positions( nx, ny, dx, dy);
    // Application des conditions aux limites
    // Application des conditions aux limites correctes
u_star_.fillBoundaryGhostCells(BCType::VelocityU, 1.0); // U avec couvercle mobile u = 1.0
v_star_.fillBoundaryGhostCells(BCType::VelocityV, 0.0); // V avec couvercle fixe v = 0.0
    // u_star_.fillBoundaryGhostCells();
    // v_star_.fillBoundaryGhostCells();
}




void NavStokesSolver::step(int step_number) {
    // 1. Calcul de la viscosité turbulente (LES)
    // if (les_model_) {
    //     les_model_->computeTurbulentViscosity(u_, v_, nu_t_, mesh_.dx(), mesh_.dy());
    // }

    // 2. Étape de prédiction : calcul de u_star et v_star
    computePredictorVelocity();
//     if(step_number==500){     
//         print_eight_key_positions(mesh_.nx(), mesh_.ny(), mesh_.dx(), mesh_.dy());
// }

    // 3 & 4. Projection complète (Poisson + Correction de vitesse + Ghost cells)
    projection_method_.project(u_star_, v_star_, u_, v_, p_, mesh_, poisson_solver_);

}
/* void NavStokesSolver::step(int step_number)
{
    // 1. Calcul du prédicteur
    computePredictorVelocity();

    // 2. Diagnostic uniquement au pas 500
    if (step_number == 500) {
        print_eight_key_positions(
            mesh_.nx(), mesh_.ny(),
            mesh_.dx(), mesh_.dy(),
            nu_mol_
        );
    }

    // 3. Projection de pression
    projection_method_.project(
        u_star_, v_star_,
        u_, v_, p_,
        mesh_, poisson_solver_
    );
}*/


void NavStokesSolver::writeVTK(int step) {
    // Avant : "output_" + std::to_string(step) + ".vtk"
    // ✅ Après : ajouter "output/" au début du chemin
    std::string filename = "output/output_" + std::to_string(step) + ".vtk";
    
    vtk_writer_.write(filename, mesh_, u_, v_, p_);
}

/**
 * @brief Affiche les contributions visqueuses (diffusion normale x, y, croisée et totale)
 *        aux 8 positions clés du domaine (maillage nx * ny).
 * 
 * Note : Cette fonction suppose que les tableaux/objets u_, v_ et nu_t_ 
 * sont accessibles dans le scope où cette fonction est appelée.
 */
void NavStokesSolver::print_eight_key_positions(int nx, int ny, double dx, double dy,double nu_mol) {
    struct KeyPos {
        std::string name;
        int i;
        int j;
    };

const int nx_phys = nx - 2;  // 48 cellules physiques
const int ny_phys = ny - 2;  // 48 cellules physiques

const int i_mid = nx_phys / 2;  // 24
const int j_mid = ny_phys / 2;  // 24

std::vector<KeyPos> positions = {
    {"Sud-Ouest",    1,       1},
    {"Sud-Est",      nx_phys, 1},
    {"Milieu Sud",   i_mid,   1},
    {"Milieu Ouest", 1,       j_mid},
    {"Milieu Est",   nx_phys, j_mid},
    {"Nord-Ouest",   1,       ny_phys},
    {"Milieu Nord",  i_mid,   ny_phys},
    {"Nord-Est",     nx_phys, ny_phys}
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




} // namespace drone