#include "core/Mesh.hpp"
#include "solvers/NavStokesSolver.hpp"
#include "physics/SmagorinskyModel.hpp"
#include "solvers/NavStokesSolver.hpp"

#include <iostream>
using namespace drone; 
int main() {
    std::cout << "=== Simulation CFD - Cavité entraînée ===" << std::endl;

    const int nx = 50, ny = 50;
    const double Lx = 1.0, Ly = 1.0;
    Mesh mesh(nx, ny, Lx, Ly);

    const double dt = 0.001;
    const double nu_mol = 0.001;
    const int max_steps = 500;
    const int save_every = 50;

    auto les_model = std::make_unique<SmagorinskyModel>(0.18); 
    NavStokesSolver solver(mesh, dt, nu_mol, std::move(les_model));

    for (int step = 0; step <= max_steps; ++step) {
        solver.step();
        if (step % save_every == 0) {
            std::cout << "Pas " << step << " / " << max_steps << " terminé." << std::endl;
        }
    }

    std::cout << "=== Simulation terminée ===" << std::endl;
    return 0;
}