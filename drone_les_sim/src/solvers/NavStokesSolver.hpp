#ifndef NAV_STOKES_SOLVER_HPP
#define NAV_STOKES_SOLVER_HPP

#include "core/Mesh.hpp"
#include "core/Field2D.hpp"
#include "solvers/PoissonSolver.hpp"
#include "physics/ProjectionMethod.hpp"
#include "physics/TurbulenceModel.hpp"
#include "io/VTKWriter.hpp"
#include <memory>
#include <string>

namespace drone {

class NavStokesSolver {
private:
    // 1. Géométrie et maillage
    Mesh mesh_;

    // 2. Champs physiques (Vitesse u, v et Pression p)
    Field2D u_, v_, p_;
    
    // 3. Champs intermédiaires de travail
    Field2D u_star_, v_star_;
    Field2D nu_t_;

    // 4. Moteurs algorithmiques et modèles
    std::unique_ptr<TurbulenceModel> les_model_;
    PoissonSolver poisson_solver_;
    ProjectionMethod projection_method_;
    VTKWriter vtk_writer_;

    // 5. Paramètres physiques et temporels
    double dt_;
    double nu_mol_;

    

public:
    // Constructeur
    NavStokesSolver(const Mesh& mesh, double dt, double nu_mol, 
                    std::unique_ptr<TurbulenceModel> model)
        : mesh_(mesh),
          u_(mesh.nx(), mesh.ny(), 0.0),
          v_(mesh.nx(), mesh.ny(), 0.0),
          p_(mesh.nx(), mesh.ny(), 0.0),
          u_star_(mesh.nx(), mesh.ny(), 0.0),
          v_star_(mesh.nx(), mesh.ny(), 0.0),
          nu_t_(mesh.nx(), mesh.ny(), 0.0),
          les_model_(std::move(model)),
          poisson_solver_(1000, 1e-5, 1.7),
          projection_method_(dt, 1.0),
          dt_(dt),
          nu_mol_(nu_mol)
    {}

    // Effectue un pas de temps complet
     double getMaxDivU() const {
        return projection_method_.getMaxDivU();
    }
    void step();

   // ✅ Après :
void writeVTK(int step_number) ;

// --- ACCESSEURS CORRIGÉS ---
    const Field2D& getU() const { return u_; }
    const Field2D& getV() const { return v_; }
    const Field2D& getP() const { return p_; }
    const Field2D& getUStar() const { return u_star_; }
    const Field2D& getVStar() const { return v_star_; }
 

private:
    // Étape 1 : Calcul de la vitesse prédictive
    void computePredictorVelocity();



};



} // namespace drone

#endif // NAV_STOKES_SOLVER_HPP