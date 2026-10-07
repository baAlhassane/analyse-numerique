#include "core/Mesh.hpp"
#include "solvers/NavStokesSolver.hpp"
#include "physics/SmagorinskyModel.hpp"
#include <memory>
#include <iostream>
#include <filesystem>
#include <cmath>
#include <algorithm>

using namespace drone;

int main() {
    std::cout << "=== Simulation CFD - Cavité entraînée (export VTK) ==="
              << std::endl;

    std::filesystem::create_directory("output");

    const int nx = 50, ny = 50;
    const double Lx = 1.0, Ly = 1.0;

    Mesh mesh(nx, ny, Lx, Ly);

    const double dt = 0.001;
    const double nu_mol = 0.001;
    const int max_steps = 500;
    const int save_every = 50;

    auto les_model = std::make_unique<SmagorinskyModel>(0.18);

    NavStokesSolver solver(
        mesh,
        dt,
        nu_mol,
        std::move(les_model)
    );

    // ============================================================
    // DIAGNOSTIC : calcul de la divergence maximale
    // ============================================================


//    auto calculerMaxDivergence = [&](const auto& u, const auto& v) {
//     double max_div = 0.0;

//     for (int j = 1; j < ny - 1; ++j) {
//         for (int i = 1; i < nx - 1; ++i) {
//             double u_e = 0.5 * (u(i, j) + u(i + 1, j));
//             double u_w = 0.5 * (u(i - 1, j) + u(i, j));
//             double v_n = 0.5 * (v(i, j) + v(i, j + 1));
//             double v_s = 0.5 * (v(i, j - 1) + v(i, j));

//             double div = (u_e - u_w) / mesh.dx() + (v_n - v_s) / mesh.dy();
//             max_div = std::max(max_div, std::abs(div));
//         }
//     }

//      return max_div;
//  };

    auto calculerMaxDivergence = [&](const auto& u, const auto& v) {

    const int nx = mesh.nx();
    const int ny = mesh.ny();
    const double dx = mesh.dx();
    const double dy = mesh.dy();

    double max_div = 0.0;
    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            // Reconstruction des flux aux faces
            double u_e = 0.5 * (u(i, j) + u(i + 1, j));
            double u_w = 0.5 * (u(i - 1, j) + u(i, j));
            double v_n = 0.5 * (v(i, j) + v(i, j + 1));
            double v_s = 0.5 * (v(i, j - 1) + v(i, j));

            double div_cell = std::abs((u_e - u_w) / dx + (v_n - v_s) / dy);
            max_div = std::max(max_div, div_cell);
        }
    }
    //std::cout << "  Max |div(U)|: " << max_div << std::endl;


        return max_div;
    };


    // ============================================================
    // SAUVEGARDE DE L'ÉTAT INITIAL
    // ============================================================

    solver.writeVTK(0);


    // ============================================================
    // FONCTION D'INSPECTION DES CHAMPS
    // ============================================================

    auto inspecterChamps = [&](int step) {

        const auto& u = solver.getU();
        const auto& v = solver.getV();
        const auto& p = solver.getP();
        const auto& u_star = solver.getUStar();
        const auto& v_star = solver.getVStar();
    

        double min_u = 1e9;
        double max_u = -1e9;
        double max_u_star= -1e9;
        double sum_u = 0.0;

        double min_v = 1e9;
        double max_v = -1e9;
        double max_v_star= -1e9;
        double sum_v = 0.0;

        double min_p = 1e9;
        double max_p = -1e9;


        // --------------------------------------------------------
        // IMPORTANT :
        // On exclut les ghost cells des statistiques.
        // --------------------------------------------------------

        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {

                double val_u = u(i, j);
                double val_v = v(i, j);
                double val_p = p(i, j);

               double  val_u_max= u_star(i,j);
               double  val_v_max= v_star(i,j);

                min_u = std::min(min_u, val_u);
                max_u = std::max(max_u, val_u);
                max_u_star=std::max(max_u_star, val_u_max);
                sum_u += std::abs(val_u);

                min_v = std::min(min_v, val_v);
                max_v = std::max(max_v, val_v);
                max_v_star=std::max(max_v_star, val_v_max);
                sum_v += std::abs(val_v);

                min_p = std::min(min_p, val_p);
                max_p = std::max(max_p, val_p);
            }
        }


        int total_points = (nx - 2) * (ny - 2);


        // --------------------------------------------------------
        // DIVERGENCE MAXIMALE
        // --------------------------------------------------------
        double  max_div_ustar= calculerMaxDivergence(u_star, v_star);

        double max_div =  solver.getMaxDivU();//calculerMaxDivergence(u, v);


        // --------------------------------------------------------
        // AFFICHAGE
        // --------------------------------------------------------

        std::cout << "\n----------------------------------------"
                  << std::endl;

        std::cout << "--- INSPECTION PAS "
                  << step
                  << " ---"
                  << std::endl;

        std::cout << "  U (Vx) -> Min: "
                  << min_u
                  << " | Max: "
                  << max_u
                  << " | Moyenne |U|: "
                  << sum_u / total_points
                  << std::endl;

        std::cout << "  V (Vy) -> Min: "
                  << min_v
                  << " | Max: "
                  << max_v
                  << " | Moyenne |V|: "
                  << sum_v / total_points
                  << std::endl;

        std::cout << "  Pression -> Min: "
                  << min_p
                  << " | Max: "
                  << max_p
                  << std::endl;

         std::cout << "  Max |div(U_Star)|: "
                  << max_div_ustar
                  << std::endl;

        std::cout << "  Max |div(U)|: "
                  <<  max_div 
                  << std::endl;


        // --------------------------------------------------------
        // POINT CENTRAL
        // --------------------------------------------------------

        int cx = nx / 2;
        int cy = ny / 2;

        std::cout << "  [Centre ("
                  << cx << "," << cy
                  << ")]       u: "
                  << u(cx, cy)
                  << " | v: "
                  << v(cx, cy)
                  << " | p: "
                  << p(cx, cy)
                  << std::endl;


        // --------------------------------------------------------
        // PAROI SUPÉRIEURE
        // --------------------------------------------------------

        int top_x = nx / 2;

        // Cellule physique juste sous le mur
        int top_cell_y = ny - 2;

        // Ghost cell au-dessus du mur
        int top_ghost_y = ny - 1;


        // Reconstruction de la vitesse exactement à la paroi
        double u_wall =
            0.5 * (
                u(top_x, top_cell_y)
                +
                u(top_x, top_ghost_y)
            );

        double v_wall =
            0.5 * (
                v(top_x, top_cell_y)
                +
                v(top_x, top_ghost_y)
            );


        std::cout << "  [Paroi Haut]       "
                  << "u_wall: "
                  << u_wall
                  << " | v_wall: "
                  << v_wall
                  << " | p: "
                  << p(top_x, top_cell_y)
                  << std::endl;

        std::cout << "----------------------------------------\n"
                  << std::endl;
    };


    // ============================================================
    // ÉTAT INITIAL
    // ============================================================

    std::cout << "\n>>> État Initial (t = 0) :"
              << std::endl;

    inspecterChamps(0);


    // ============================================================
    // BOUCLE TEMPORELLE
    // ============================================================

    for (int step = 1; step <= max_steps; ++step) {

        solver.step();

        if (step % save_every == 0) {

            solver.writeVTK(step);

            std::cout << "Pas "
                      << step
                      << " / "
                      << max_steps
                      << " exporté vers output/"
                      << std::endl;

            inspecterChamps(step);
        }
    }


    std::cout << "=== Simulation terminée avec succès ==="
              << std::endl;

    return 0;
}






























// #include "core/Mesh.hpp"
// #include "solvers/NavStokesSolver.hpp"
// #include "physics/SmagorinskyModel.hpp"
// #include <memory>
// #include <iostream>
// #include <filesystem>
// #include <cmath>
// #include <algorithm>
// #include <functional>

// using namespace drone;






// int main() {
//     std::cout << "=== Simulation CFD - Cavité entraînée (export VTK) ===" << std::endl;

//     std::filesystem::create_directory("output");

//     const int nx = 50, ny = 50;
//     const double Lx = 1.0, Ly = 1.0;
//     Mesh mesh(nx, ny, Lx, Ly);

//     const double dt = 0.001;
//     const double nu_mol = 0.001;
//     const int max_steps = 500;
//     const int save_every = 50;

//     auto les_model = std::make_unique<SmagorinskyModel>(0.18);
//     NavStokesSolver solver(mesh, dt, nu_mol, std::move(les_model));

//     // Sauvegarde de l'état initial (t = 0)
//     solver.writeVTK(0);

//     // Fonction d'inspection des champs Field2D
//     auto inspecterChamps = [&](int step) {
//         const auto& u = solver.getU();
//         const auto& v = solver.getV();
//         const auto& p = solver.getP();

//         double min_u = 1e9, max_u = -1e9, sum_u = 0.0;
//         double min_v = 1e9, max_v = -1e9, sum_v = 0.0;
//         double min_p = 1e9, max_p = -1e9;

//         for (int j = 0; j < ny; ++j) {
//             for (int i = 0; i < nx; ++i) {
//                 // Utilisation de l'opérateur (i, j) classique de Field2D
//                 double val_u = u(i, j);
//                 double val_v = v(i, j);
//                 double val_p = p(i, j);

//                 min_u = std::min(min_u, val_u);
//                 max_u = std::max(max_u, val_u);
//                 sum_u += std::abs(val_u);

//                 min_v = std::min(min_v, val_v);
//                 max_v = std::max(max_v, val_v);
//                 sum_v += std::abs(val_v);

//                 min_p = std::min(min_p, val_p);
//                 max_p = std::max(max_p, val_p);
//             }
//         }

//         int total_points = nx * ny;

//         std::cout << "\n----------------------------------------" << std::endl;
//         std::cout << "--- INSPECTION PAS " << step << " ---" << std::endl;
//         std::cout << "  U (Vx) -> Min: " << min_u << " | Max: " << max_u << " | Moyenne |U|: " << sum_u / total_points << std::endl;
//         std::cout << "  V (Vy) -> Min: " << min_v << " | Max: " << max_v << " | Moyenne |V|: " << sum_v / total_points << std::endl;
//         std::cout << "  Pression -> Min: " << min_p << " | Max: " << max_p << std::endl;

//         // Points de contrôle : centre et paroi supérieure (Lid)
//         int cx = nx / 2, cy = ny / 2;
//         int top_x = nx / 2, top_y = ny - 1;

//         std::cout << "  [Centre (" << cx << "," << cy << ")]       u: " << u(cx, cy) << " | v: " << v(cx, cy) << " | p: " << p(cx, cy) << std::endl;
//         std::cout << "  [Paroi Haut (" << top_x << "," << top_y << ")] u: " << u(top_x, top_y) << " | v: " << v(top_x, top_y) << " | p: " << p(top_x, top_y) << std::endl;
//         std::cout << "----------------------------------------\n" << std::endl;
//     };

//     std::cout << "\n>>> État Initial (t = 0) :" << std::endl;
//     inspecterChamps(0);

//     for (int step = 1; step <= max_steps; ++step) {
//         solver.step();
        
//         if (step % save_every == 0) {
//             solver.writeVTK(step);
//             std::cout << "Pas " << step << " / " << max_steps << " exporté vers output/" << std::endl;
//             inspecterChamps(step);
//         }
//     }

//     std::cout << "=== Simulation terminée avec succès ===" << std::endl;




    
//     return 0;
// }





























// // #include "core/Mesh.hpp"
// // #include "solvers/NavStokesSolver.hpp"
// // #include "physics/SmagorinskyModel.hpp"
// // #include <memory>
// // #include <iostream>
// // #include <filesystem>

// // using namespace drone;

// // int main() {
// //     std::cout << "=== Simulation CFD - Cavité entraînée (export VTK) ===" << std::endl;

// //     // Création du dossier d'export s'il n'existe pas
// //     std::filesystem::create_directory("output");

// //     const int nx = 50, ny = 50;
// //     const double Lx = 1.0, Ly = 1.0;
// //     Mesh mesh(nx, ny, Lx, Ly);

// //     const double dt = 0.001;
// //     const double nu_mol = 0.001;
// //     const int max_steps = 500;
// //     const int save_every = 50;

// //     auto les_model = std::make_unique<SmagorinskyModel>(0.18);
// //     NavStokesSolver solver(mesh, dt, nu_mol, std::move(les_model));

// //     // Sauvegarde de l'état initial (t = 0)
// //     solver.writeVTK(0);

// //     for (int step = 1; step <= max_steps; ++step) {
// //         solver.step();
        
// //         if (step % save_every == 0) {
// //             solver.writeVTK(step);
// //             std::cout << "Pas " << step << " / " << max_steps << " exporté vers output/" << std::endl;
// //         }
// //     }

// //     std::cout << "=== Simulation terminée avec succès ===" << std::endl;
// //     return 0;
// // }