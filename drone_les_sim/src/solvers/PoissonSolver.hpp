#ifndef POISSON_SOLVER_HPP
#define POISSON_SOLVER_HPP

#include "core/Field2D.hpp"
#include "core/Mesh.hpp"
#include <cmath>
#include <algorithm>


   namespace drone {
     
class PoissonSolver {
private:
    int max_iter_;
    double tolerance_;
    double omega_; // Facteur de sur-relaxation SOR (ex: 1.7)

public:
    PoissonSolver(int max_iter = 1000, double tolerance = 1e-5, double omega = 1.7)
        : max_iter_(max_iter), tolerance_(tolerance), omega_(omega) {}

    void solve(const Field2D& rhs, Field2D& p, const Mesh& mesh) {
        int nx = mesh.nx();
        int ny = mesh.ny();
        double dx = mesh.dx();
        double dx2 = dx * dx;

        for (int iter = 0; iter < max_iter_; ++iter) {
            double max_error = 0.0;

            for (int j = 1; j < ny - 1; ++j) {
                for (int i = 1; i < nx - 1; ++i) {

                    double p_old = p(i, j);

                    // Valeur de Gauss-Seidel
                    double p_gs = 0.25 * (p(i + 1, j) + p(i - 1, j) +
                                         p(i, j + 1) + p(i, j - 1) -
                                         dx2 * rhs(i, j));

                    // Mise à jour SOR (Sur-relaxation)
                    p(i, j) = (1.0 - omega_) * p_old + omega_ * p_gs;

                    // Calcul de l'erreur absolue pour le critère de convergence
                    max_error = std::max(max_error, std::abs(p(i, j) - p_old));
                }
            }

            // Conditions aux limites Neumann sur la pression (dp/dn = 0 sur les parois)
            applyPressureBC(p, nx, ny);

            // Test de convergence
            if (max_error < tolerance_) {
                break; // Le système a convergé !
            }
        }
    }

private:
    void applyPressureBC(Field2D& p, int nx, int ny) {
        // Neumann Homogène : dp/dn = 0 (les cellules fantômes prennent la valeur adjacente)
        for (int i = 0; i < nx; ++i) {
            p(i, 0) = p(i, 1);           // Paroi Sud
            p(i, ny - 1) = p(i, ny - 2); // Paroi Nord
        }
        for (int j = 0; j < ny; ++j) {
            p(0, j) = p(1, j);           // Paroi Ouest
            p(nx - 1, j) = p(nx - 2, j); // Paroi Est
        }
    }
};

   } // namespace drone 

#endif // POISSON_SOLVER_HPP