#ifndef PROJECTION_METHOD_HPP
#define PROJECTION_METHOD_HPP

#include "core/Field2D.hpp"
#include "core/Mesh.hpp"
#include "solvers/PoissonSolver.hpp"



class ProjectionMethod {
private:
    double dt_;
    double rho_;

public:
    ProjectionMethod(double dt, double rho = 1.0)
        : dt_(dt), rho_(rho) {}

    // Exécute l'étape de correction par projection
    void project(const Field2D& u_star, const Field2D& v_star,
                 Field2D& u_next, Field2D& v_next,
                 Field2D& p, const Mesh& mesh,
                 PoissonSolver& poisson_solver) 
    {
        int nx = mesh.nx();
        int ny = mesh.ny();
        double dx = mesh.dx();
        double dy = mesh.dy();

        // 1. Instanciation du terme source RHS (un seul champ scalaire !)
        Field2D rhs(nx, ny, 0.0);

        // 2. Calcul du RHS : div(u*) * (rho / dt)
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                
                // Divergence centrée de la vitesse prédictive u*
                double div_u_star = (u_star(i + 1, j) - u_star(i - 1, j)) / (2.0 * dx)
                                  + (v_star(i, j + 1) - v_star(i, j - 1)) / (2.0 * dy);

                rhs(i, j) = (rho_ / dt_) * div_u_star;
            }
        }

        // 3. Résolution du Laplacien de pression : nabla^2(p) = RHS
        // On passe le terme source 'rhs', le champ de pression 'p' à mettre à jour, et le 'mesh'
        poisson_solver.solve(rhs, p, mesh);

        // 4. Correction de la vitesse : u_next = u* - (dt / rho) * grad(p)
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                
                double dp_dx = (p(i + 1, j) - p(i - 1, j)) / (2.0 * dx);
                double dp_dy = (p(i, j + 1) - p(i, j - 1)) / (2.0 * dy);

                u_next(i, j) = u_star(i, j) - (dt_ / rho_) * dp_dx;
                v_next(i, j) = v_star(i, j) - (dt_ / rho_) * dp_dy;
            }
        }

        // 5. Mise à jour des mailles fantômes sur le champ final
        u_next.fillBoundaryGhostCells();
        v_next.fillBoundaryGhostCells();
    }
};

#endif // PROJECTION_METHOD_HPP