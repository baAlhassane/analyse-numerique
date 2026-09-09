#ifndef MESH_HPP
#define MESH_HPP


class Mesh {
private:
    int nx_, ny_;
    double Lx_, Ly_;
    double dx_, dy_;

public:
    Mesh(int nx, int ny, double Lx, double Ly)
        : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly) 
    {
        dx_ = Lx_ / static_cast<double>(nx_);
        dy_ = Ly_ / static_cast<double>(ny_);
    }

    // Getters
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    double dx() const { return dx_; }
    double dy() const { return dy_; }
    double Lx() const { return Lx_; }
    double Ly() const { return Ly_; }

    // Coordonnées des centres de mailles (0.5 * dx pour le centre)
    double x(int i) const { return (i + 0.5) * dx_; }
    double y(int j) const { return (j + 0.5) * dy_; }
};


#endif // MESH_HPP