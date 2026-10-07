#include "io/VTKWriter.hpp"
#include <fstream>
#include <iostream>

namespace drone {

void VTKWriter::write(const std::string& filename, 
                      const Mesh& mesh, 
                      const Field2D& u, 
                      const Field2D& v, 
                      const Field2D& p) 
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << filename << std::endl;
        return;
    }

    int nx = mesh.nx();
    int ny = mesh.ny();

    file << "# vtk DataFile Version 3.0\n";
    file << "DroneCFD Output\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "ORIGIN 0 0 0\n";
    file << "SPACING " << mesh.dx() << " " << mesh.dy() << " 1\n";

    file << "POINT_DATA " << nx * ny << "\n";

    // Pression
    file << "SCALARS pressure double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << p(i, j) << "\n";
        }
    }

    // Vitesse (u, v, 0)
    file << "VECTORS velocity double\n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            file << u(i, j) << " " << v(i, j) << " 0.0\n";
        }
    }
}

} // namespace drone