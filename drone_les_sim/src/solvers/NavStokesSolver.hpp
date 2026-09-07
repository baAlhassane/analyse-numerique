// Solvers/NavStokesSolver.hpp
class NavStokesSolver {
private:
    Grid grid_;
    Field2D u_, v_, p_;
    Field2D u_star_, v_star_;
    Field2D nu_t_;
    
    std::unique_ptr<TurbulenceModel> les_model_;
    PoissonSolver poisson_solver_;
    VTKWriter vtk_writer_;

    double dt_, nu_mol_;

public:
    NavStokesSolver(const Grid& grid, double nu_mol, std::unique_ptr<TurbulenceModel> model)
        : grid_(grid), nu_mol_(nu_mol), les_model_(std::move(model)) {}

    void step() {
        // 1. Calculer nu_t avec le modèle LES choisi
        les_model_->computeTurbulentViscosity(u_, v_, nu_t_, grid_.dx(), grid_.dy());

        // 2. Calculer les vitesses intermédiaires u* et v* (Advection + Viscosité)
        computePredictorVelocity();

        // 3. Résoudre l'équation de Poisson pour P : div(grad P) = div(u*) / dt
        poisson_solver_.solve(u_star_, v_star_, p_, dt_);

        // 4. Corriger les vitesses u et v avec le gradient de P
        correctVelocity();
    }

private:
    void computePredictorVelocity();
    void correctVelocity();
};