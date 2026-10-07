

#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "core/Field2D.hpp"
#include "core/Mesh.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace drone {

class PoissonSolver {
private:
    int max_iter_;
    double tolerance_;
    double omega_;

public:
    PoissonSolver(int max_iter = 5000, double tolerance = 1e-6, double omega = 1.3)
        : max_iter_(max_iter), tolerance_(tolerance), omega_(omega) {}





void solve(const Field2D& rhs, Field2D& p, const Mesh& mesh)
{
    const int nx = mesh.nx();
    const int ny = mesh.ny();

    const double dx = mesh.dx();
    const double dy = mesh.dy();

    const double inv_dx2 = 1.0 / (dx * dx);
    const double inv_dy2 = 1.0 / (dy * dy);

    const double denominator =
        2.0 * inv_dx2 + 2.0 * inv_dy2;

    // ---------------------------------------------------------
    // 1. Correction de compatibilité du second membre
    // ---------------------------------------------------------

    double sum_rhs = 0.0;

    const int count_interior = (nx - 2) * (ny - 2);

    for (int j = 1; j < ny - 1; ++j)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            sum_rhs += rhs(i, j);
        }
    }

    const double avg_rhs =
        sum_rhs / static_cast<double>(count_interior);


    // ---------------------------------------------------------
    // 2. Itérations SOR
    // ---------------------------------------------------------

    for (int iter = 0; iter < max_iter_; ++iter)
    {
        // ---------------------------------------------
        // Gauss-Seidel / SOR
        // ---------------------------------------------

        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                const double p_old = p(i, j);

                // RHS réellement utilisé
                const double rhs_corr =
                    rhs(i, j) - avg_rhs;

                const double p_gs =
                    (
                        (p(i + 1, j) + p(i - 1, j)) * inv_dx2
                      + (p(i, j + 1) + p(i, j - 1)) * inv_dy2
                      - rhs_corr
                    )
                    / denominator;

                p(i, j) =
                    (1.0 - omega_) * p_old
                    + omega_ * p_gs;
            }
        }

        // ---------------------------------------------
        // Condition de jauge
        // ---------------------------------------------

        const double p_ref = p(1, 1);

        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                p(i, j) -= p_ref;
            }
        }

        // ---------------------------------------------
        // Conditions aux limites
        // ---------------------------------------------

        applyPressureBC(p, nx, ny);

        // ---------------------------------------------
        // Calcul du résidu ||L5 p - RHS_corr||inf
        // ---------------------------------------------

        double max_residual = 0.0;

        for (int j = 1; j < ny - 1; ++j)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                const double rhs_corr =
                    rhs(i, j) - avg_rhs;

                const double lap_p =
                    (p(i + 1, j)
                   - 2.0 * p(i, j)
                   + p(i - 1, j)) * inv_dx2

                  + (p(i, j + 1)
                   - 2.0 * p(i, j)
                   + p(i, j - 1)) * inv_dy2;

                const double residual =
                    lap_p - rhs_corr;

                max_residual =
                    std::max(
                        max_residual,
                        std::abs(residual)
                    );
            }
        }

        // ---------------------------------------------
        // Critère d'arrêt
        // ---------------------------------------------

        if (max_residual < tolerance_)
        {
            break;
        }
    }
}

    // void solve(const Field2D& rhs, Field2D& p, const Mesh& mesh) {
    //     const int nx = mesh.nx();
    //     const int ny = mesh.ny();
    //     const double dx = mesh.dx();
    //     const double dy = mesh.dy();

    //     const double inv_dx2 = 1.0 / (dx * dx);
    //     const double inv_dy2 = 1.0 / (dy * dy);
    //     const double denominator = 2.0 * inv_dx2 + 2.0 * inv_dy2;

    //     // Condition de compatibilité pour Neumann pur : moy(RHS) = 0
    //     double sum_rhs = 0.0;
    //     int count_interior = (nx - 2) * (ny - 2);

    //     for (int j = 1; j < ny - 1; ++j) {
    //         for (int i = 1; i < nx - 1; ++i) {
    //             sum_rhs += rhs(i, j);
    //         }
    //     }
        
    //      double avg_rhs = sum_rhs / count_interior; 
    //     // std::cout << "sum_rhs = " << sum_rhs
    //     //   << ", avg_rhs = " << avg_rhs << std::endl;

    //     // Boucle SOR
    //     for (int iter = 0; iter < max_iter_; ++iter) {
    //         double max_error = 0.0;

    //         for (int j = 1; j < ny - 1; ++j) {
    //             for (int i = 1; i < nx - 1; ++i) {
    //                 const double p_old = p(i, j);
    //                 const double rhs_corr = rhs(i, j) - avg_rhs;

    //                 const double p_gs = ((p(i + 1, j) + p(i - 1, j)) * inv_dx2
    //                                   + (p(i, j + 1) + p(i, j - 1)) * inv_dy2
    //                                   - rhs_corr) / denominator;

    //                 p(i, j) = (1.0 - omega_) * p_old + omega_ * p_gs;
    //                 max_error = std::max(max_error, std::abs(p(i, j) - p_old));
    //             }
    //         }

    //         // Fixer la référence de pression (ancrage à 0 sur p(1,1))
    //         double p_ref = p(1, 1);
    //         for (int j = 0; j < ny; ++j) {
    //             for (int i = 0; i < nx; ++i) {
    //                 p(i, j) -= p_ref;
    //             }
    //         }

    //         // Mettre à jour les mailles fantômes
    //         applyPressureBC(p, nx, ny);

    //         if (max_error < tolerance_) {
    //             break;
    //         }
    //     }
    // }

private:
    void applyPressureBC(Field2D& p, int nx, int ny) {
        for (int i = 0; i < nx; ++i) {
            p(i, 0) = p(i, 1);
            p(i, ny - 1) = p(i, ny - 2);
        }
        for (int j = 0; j < ny; ++j) {
            p(0, j) = p(1, j);
            p(nx - 1, j) = p(nx - 2, j);
        }
    }
};

} // namespace drone

#endif