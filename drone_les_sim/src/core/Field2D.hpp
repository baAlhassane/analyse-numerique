// Core/Field.hpp
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
    
    void fillBoundaryGhostCells(); // Applique les CL_
};