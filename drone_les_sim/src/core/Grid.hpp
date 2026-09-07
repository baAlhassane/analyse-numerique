#ifndef GRID_HPP
#define GRID_HPP

class Grid {
private:
    int nx_, ny_;       // Nombre de mailles (x et y)
    double Lx_, Ly_;    // Longueur physique du domaine (mètres)
    double dx_, dy_;    // Pas d'espace (mètres)

public:
    // Constructeur : calcule automatiquement dx et dy
    Grid(int nx, int ny, double Lx, double Ly)
        : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly) 
    {
        dx_ = Lx_ / static_cast<double>(nx_);
        dy_ = Ly_ / static_cast<double>(ny_);
    }

    // --- GETTERS (Constants et rapides) ---
    inline int nx() const { return nx_; }
    inline int ny() const { return ny_; }

    inline double Lx() const { return Lx_; }
    inline double Ly() const { return Ly_; }

    inline double dx() const { return dx_; }
    inline double dy() const { return dy_; }

    // Obtenir la position x réelle au centre de la maille i
    inline double x(int i) const {
        return (i + 0.5) * dx_;
    }

    // Obtenir la position y réelle au centre de la maille j
    inline double y(int j) const {
        return (j + 0.5) * dy_;
    }
};

#endif // GRID_HPP