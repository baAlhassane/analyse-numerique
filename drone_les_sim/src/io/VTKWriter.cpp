#include "io/VTKWriter.hpp"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace drone {

void VTKWriter::write(const std::string& filename,
                      const Mesh& mesh,
                      const Field2D& u,
                      const Field2D& v,
                      const Field2D& p)
{
    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Erreur : impossible d'ouvrir le fichier "
                  << filename << std::endl;
        return;
    }

    const int nx = mesh.nx();
    const int ny = mesh.ny();

    const double dx = mesh.dx();
    const double dy = mesh.dy();

    // -------------------------------------------------------------------------
    // Cellules physiques :
    //
    // i = 1 ... nx-2
    // j = 1 ... ny-2
    //
    // Les indices 0 et nx-1 / ny-1 sont les ghost cells.
    // -------------------------------------------------------------------------

    const int nx_phys = nx - 2;
    const int ny_phys = ny - 2;

    if (nx_phys <= 0 || ny_phys <= 0) {
        std::cerr << "Erreur : dimensions physiques invalides."
                  << std::endl;
        return;
    }

    // Coordonnées du premier centre physique.
    const double x0 = mesh.x(1);
    const double y0 = mesh.y(1);

    const int n_points = nx_phys * ny_phys;

    // -------------------------------------------------------------------------
    // VTK
    //
    // Ici, chaque point VTK représente directement le centre d'une cellule
    // physique du solveur.
    //
    // Cela permet de conserver exactement les coordonnées :
    //
    // x = mesh.x(i)
    // y = mesh.y(j)
    //
    // sans exporter les ghost cells.
    // -------------------------------------------------------------------------

    file << "# vtk DataFile Version 3.0\n";
    file << "DroneCFD - Cavity Flow\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    file << "DIMENSIONS "
         << nx_phys << " "
         << ny_phys << " 1\n";

    file << "ORIGIN "
         << std::setprecision(12)
         << x0 << " "
         << y0 << " 0\n";

    file << "SPACING "
         << dx << " "
         << dy << " 1\n";

    file << "POINT_DATA " << n_points << "\n";

    // -------------------------------------------------------------------------
    // 1. PRESSION
    // -------------------------------------------------------------------------

    file << "\n";
    file << "SCALARS pressure double 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {
            file << p(i, j) << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // 2. VITESSE
    // -------------------------------------------------------------------------

    file << "\n";
    file << "VECTORS velocity double\n";

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {

            file << u(i, j) << " "
                 << v(i, j) << " "
                 << 0.0 << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // 3. NORME DE LA VITESSE
    // -------------------------------------------------------------------------

    file << "\n";
    file << "SCALARS velocity_magnitude double 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {

            const double ui = u(i, j);
            const double vi = v(i, j);

            const double magnitude =
                std::sqrt(ui * ui + vi * vi);

            file << magnitude << "\n";
        }
    }

    // -------------------------------------------------------------------------
    // 4. VORTICITE
    //
    // omega = dv/dx - du/dy
    //
    // Différences centrées.
    // -------------------------------------------------------------------------

    file << "\n";
    file << "SCALARS vorticity double 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int j = 1; j < ny - 1; ++j) {
        for (int i = 1; i < nx - 1; ++i) {

            const double dv_dx =
                (v(i + 1, j) - v(i - 1, j))
                / (2.0 * dx);

            const double du_dy =
                (u(i, j + 1) - u(i, j - 1))
                / (2.0 * dy);

            const double omega = dv_dx - du_dy;

            file << omega << "\n";
        }
    }

    file.close();

    std::cout << "VTK écrit : " << filename
              << " (" << nx_phys << " x "
              << ny_phys << " points physiques)"
              << std::endl;
}

} // namespace drone