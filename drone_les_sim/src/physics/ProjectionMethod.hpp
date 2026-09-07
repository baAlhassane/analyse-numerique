#include "Field2D.hpp"
class ProjectionMethod {

public:
void projectionMethod(const Field2D & u_star,  const Field2D & v_star , double rho, double dt){
    int ny= u_star.nx();
    int nx= v_star.ny();
    Field2D div_u_star( u_star.nx(), v_star.ny()); 
    Field2D div_v_star( u_star.nx(), v_star.ny()); 
     Field2D  div_cellule(u_star.nx(), v_star.ny());
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx i++){

             div_u_star(i,j)= (u_star(i+i,j)- u_star(i-1,j))/(2*dx );
             div_v_star(i,j)= (v_star(i,j+1)- v_star(i,j-1))/(2*dy ) ;
             
             div_cellule(i,j)= (rho/dt )*(div_u_star(i,j)+  div_v_star(i,j));
        }

    }


   PoissonSolver.solve( pression, div_cellule, nx, ny ); 
   

}

};