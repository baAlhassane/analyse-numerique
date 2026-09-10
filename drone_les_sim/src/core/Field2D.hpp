// Core/Field.hpp

#ifndef FIELD2D_HPP
#define FIELD2D_HPP
#include <vector>

namespace drone {


class Field2D {
private:
    int nx_, ny_; // n_x longueur  suivant l'horizontale
    std::vector<double> data_;
public:
    Field2D(int nx, int ny, double init_val = 0.0) 
        : nx_(nx), ny_(ny), data_(nx * ny, init_val) {}

    inline double& operator()(int i, int j) { return data_[j * nx_ + i]; }
    inline double operator()(int i, int j) const { return data_[j * nx_ + i]; }
    inline int nx() const { return nx_; }
    inline int ny() const { return ny_; }
    
    //void fillBoundaryGhostCells(); // Applique les CL_


void fillBoundaryGhostCells() {
    // 1. Paroi SUD (j = 0) : Paroi fixe (Dirichlet u = 0 -> ghost = -interieur)
    for (int i = 0; i < nx_; ++i) {
        (*this)(i, 0) = -(*this)(i, 1);
    }

    // 2. Paroi NORD (j = ny - 1) : Couvercle mobile ou paroi fixe
    for (int i = 0; i < nx_; ++i) {
        (*this)(i, ny_ - 1) = -(*this)(i, ny_ - 2); 
        // Note : Si c'est un couvercle mobile u_wall, ce serait : 2.0 * u_wall - (*this)(i, ny_ - 2)
    }

    // 3. Bord OUEST (i = 0) : Entrée (Inflow Dirichlet) ou Paroi
    for (int j = 0; j < ny_; ++j) {
        (*this)(0, j) = (*this)(1, j); // Exemple : Neumann (dp/dn = 0 ou du/dx = 0)
    }

    // 4. Bord EST (i = nx - 1) : Sortie libre (Outflow Neumann du/dx = 0)
    for (int j = 0; j < ny_; ++j) {
        (*this)(nx_ - 1, j) = (*this)(nx_ - 2, j);
    }

    
}
   
};

}

#endif // FIELD2D_HPP