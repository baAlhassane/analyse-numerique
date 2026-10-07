
#ifndef PROJECTION_METHOD_HPP
#define PROJECTION_METHOD_HPP

#include "core/Field2D.hpp"
#include "core/Mesh.hpp"
#include "solvers/PoissonSolver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace drone {

class ProjectionMethod {
private:
    double dt_;
    double rho_;
    double max_div_u = 0.0;

public:
    ProjectionMethod(double dt, double rho = 1.0)
        : dt_(dt), rho_(rho) {}

    void project(const Field2D& u_star,
                 const Field2D& v_star,
                 Field2D& u_next,
                 Field2D& v_next,
                 Field2D& p,
                 const Mesh& mesh,
                 PoissonSolver& poisson_solver)
    {
        const int nx = mesh.nx();
        const int ny = mesh.ny();

        const double dx = mesh.dx();
        const double dy = mesh.dy();

        // =========================================================
        // 1. Calcul du RHS du Poisson
        //
        //     nabla² p = (rho / dt) * div(u*)
        //
        // =========================================================

        Field2D rhs(nx, ny, 0.0);

        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {

                // Vitesses u* aux faces
                const double u_e =
                    0.5 * (u_star(i, j) + u_star(i + 1, j));

                const double u_w =
                    0.5 * (u_star(i - 1, j) + u_star(i, j));

                // Vitesses v* aux faces
                const double v_n =
                    0.5 * (v_star(i, j) + v_star(i, j + 1));

                const double v_s =
                    0.5 * (v_star(i, j - 1) + v_star(i, j));

                // Divergence de la vitesse prédictive
                const double div_u_star =
                      (u_e - u_w) / dx
                    + (v_n - v_s) / dy;

                rhs(i, j) = (rho_ / dt_) * div_u_star;
            }
        }

        // =========================================================
        // 2. Résolution du Poisson
        //
        //     nabla² p = RHS
        //
        // =========================================================

        poisson_solver.solve(rhs, p, mesh);

        // =========================================================
        // 3. Calcul du résidu réel du Poisson
        // =========================================================

        double max_residual = 0.0;

        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {

                const double lap_p =
                      (p(i + 1, j) - 2.0 * p(i, j) + p(i - 1, j))
                        / (dx * dx)
                    + (p(i, j + 1) - 2.0 * p(i, j) + p(i, j - 1))
                        / (dy * dy);

                const double residual =
                    lap_p - rhs(i, j);

                max_residual =
                    std::max(max_residual, std::abs(residual));
            }
        }

        // =========================================================
        // 4. Correction de pression avec G_f
        // =========================================================

        

        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {

                // -------------------------------------------------
                // 4.1 Vitesses prédictives aux faces
                // -------------------------------------------------

                const double u_star_e =
                    0.5 * (u_star(i, j) + u_star(i + 1, j));

                const double u_star_w =
                    0.5 * (u_star(i - 1, j) + u_star(i, j));

                const double v_star_n =
                    0.5 * (v_star(i, j) + v_star(i, j + 1));

                const double v_star_s =
                    0.5 * (v_star(i, j - 1) + v_star(i, j));

                // -------------------------------------------------
                // 4.2 Gradient de pression G_f aux faces
                // -------------------------------------------------

                const double dp_dx_e =
                    (p(i + 1, j) - p(i, j)) / dx;

                const double dp_dx_w =
                    (p(i, j) - p(i - 1, j)) / dx;

                const double dp_dy_n =
                    (p(i, j + 1) - p(i, j)) / dy;

                const double dp_dy_s =
                    (p(i, j) - p(i, j - 1)) / dy;

                // -------------------------------------------------
                // 4.3 Correction de pression aux faces
                // -------------------------------------------------

                const double u_e =
                    u_star_e - (dt_ / rho_) * dp_dx_e;

                const double u_w =
                    u_star_w - (dt_ / rho_) * dp_dx_w;

                const double v_n =
                    v_star_n - (dt_ / rho_) * dp_dy_n;

                const double v_s =
                    v_star_s - (dt_ / rho_) * dp_dy_s;

                // -------------------------------------------------
                // 4.4 Reconstruction des vitesses au centre
                // -------------------------------------------------

                u_next(i, j) =
                    0.5 * (u_e + u_w);

                v_next(i, j) =
                    0.5 * (v_n + v_s);

                // -------------------------------------------------
                // 4.5 Divergence directement sur les flux corrigés
                // -------------------------------------------------

                const double div_corrected =
                      (u_e - u_w) / dx
                    + (v_n - v_s) / dy;

                max_div_u =
                    std::max(max_div_u,
                             std::abs(div_corrected));
            }
        }

        // =========================================================
        // 5. Conditions aux limites sur les vitesses finales
        // =========================================================

        u_next.fillBoundaryGhostCells(
            BCType::VelocityU,
            1.0
        );

        v_next.fillBoundaryGhostCells(
            BCType::VelocityV,
            0.0
        );

        // =========================================================
        // 6. Diagnostics
        // =========================================================

        // std::cout
        //     << "  Max Poisson residual = "
        //     << max_residual
        //     << std::endl;

        // std::cout
        //     << "  max_div_u face post projection = "
        //     << max_div_u
        //     << std::endl;
    }


     double getMaxDivU() const {
        return max_div_u;
    }
};

} // namespace drone

#endif // PROJECTION_METHOD_HPP

















// #ifndef PROJECTION_METHOD_HPP
// #define PROJECTION_METHOD_HPP

// #include "core/Field2D.hpp"
// #include "core/Mesh.hpp"
// #include "solvers/PoissonSolver.hpp"
// #include <iostream>


// namespace drone {
    
// class ProjectionMethod {
// private:
//     double dt_;
//     double rho_;

// public:
//     ProjectionMethod(double dt, double rho = 1.0)
//         : dt_(dt), rho_(rho) {}

//     // Exécute l'étape de correction par projection
//     void project(const Field2D& u_star, const Field2D& v_star,
//                  Field2D& u_next, Field2D& v_next,
//                  Field2D& p, const Mesh& mesh,
//                  PoissonSolver& poisson_solver) 
//     {
//         int nx = mesh.nx();
//         int ny = mesh.ny();
//         double dx = mesh.dx();
//         double dy = mesh.dy();

//         // 1. Instanciation du terme source RHS (un seul champ scalaire !)
//         Field2D rhs(nx, ny, 0.0);

//         // 2. Calcul du RHS : div(u*) * (rho / dt)
//         for (int j = 1; j < ny - 1; ++j) {
//             for (int i = 1; i < nx - 1; ++i) {
                
//                 // Divergence centrée de la vitesse prédictive u*
//                // double div_u_star = (u_star(i + 1, j) - u_star(i - 1, j)) / (2.0 * dx)
//                 //                  + (v_star(i, j + 1) - v_star(i, j - 1)) / (2.0 * dy);
//                 double u_e = 0.5 * (u_star(i, j) + u_star(i + 1, j));
//                double u_w = 0.5 * (u_star(i - 1, j) + u_star(i, j));

//               double v_n = 0.5 * (v_star(i, j) + v_star(i, j + 1));
//               double v_s = 0.5 * (v_star(i, j - 1) + v_star(i, j));

// double div_u_star = (u_e - u_w) / dx+ (v_n - v_s) / dy;

//                 rhs(i, j) = (rho_ / dt_) * div_u_star;  
//             }
//         }

//         // 3. Résolution du Laplacien de pression : nabla^2(p) = RHS
//         // On passe le terme source 'rhs', le champ de pression 'p' à mettre à jour, et le 'mesh'
//         poisson_solver.solve(rhs, p, mesh);

//         double max_residual = 0.0;

// for (int j = 1; j < ny - 1; ++j) {
//     for (int i = 1; i < nx - 1; ++i) {

//         double lap_p =
//             (p(i + 1, j) - 2.0 * p(i, j) + p(i - 1, j)) / (dx * dx)
//           + (p(i, j + 1) - 2.0 * p(i, j) + p(i, j - 1)) / (dy * dy);

//         double residual = lap_p - rhs(i, j);

//         max_residual =
//             std::max(max_residual, std::abs(residual)); 
//     }
// }

//  //std::cout << "Max Poisson residual = "         << max_residual         << std::endl;

//         // 4. Correction de la vitesse : u_next = u* - (dt / rho) * grad(p)
//         for (int j = 1; j < ny - 1; ++j) {
//             for (int i = 1; i < nx - 1; ++i) {
                
//                 double dp_dx = (p(i + 1, j) - p(i - 1, j)) / (2.0 * dx);
//                 double dp_dy = (p(i, j + 1) - p(i, j - 1)) / (2.0 * dy);

//                 u_next(i, j) = u_star(i, j) - (dt_ / rho_) * dp_dx;
//                 v_next(i, j) = v_star(i, j) - (dt_ / rho_) * dp_dy;
//             }
//         }

//         // 5. Mise à jour des mailles fantômes sur le champ final
//         // 5. Mise à jour correcte des mailles fantômes sur le champ final
// u_next.fillBoundaryGhostCells(BCType::VelocityU, 1.0);
// v_next.fillBoundaryGhostCells(BCType::VelocityV, 0.0);




// // --- NOUVEAU : Calcul et affichage de la divergence post-projection ---
//         double max_div_u = 0.0;
//         for (int j = 1; j < ny - 1; ++j) {
//             for (int i = 1; i < nx - 1; ++i) {
//                 double u_e = 0.5 * (u_next(i, j) + u_next(i + 1, j));
//                 double u_w = 0.5 * (u_next(i - 1, j) + u_next(i, j));
//                 double v_n = 0.5 * (v_next(i, j) + v_next(i, j + 1));
//                 double v_s = 0.5 * (v_next(i, j - 1) + v_next(i, j));

//                 double div_local = std::abs((u_e - u_w) / dx + (v_n - v_s) / dy);
//                 max_div_u = std::max(max_div_u, div_local);
//             }
//         }

//         std::cout<< "  max_div_u post projection  = "<<  max_div_u<<std::endl;
//         // Vous pouvez loguer ou retourner max_div_u ici
//         // u_next.fillBoundaryGhostCells();
//         // v_next.fillBoundaryGhostCells();
//     }
// };

// #endif // PROJECTION_METHOD_HPP

// } // namespace drone