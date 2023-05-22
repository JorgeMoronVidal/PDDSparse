#include "GMSolver.hpp"
#include "Monegros_Poisson.hpp"
#include "rectangle.hpp"
int main(int argc, char *argv[]){
    BVP bvp;
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, pfscalarN> scalarN_init;
    std::map<std::string, std::string> string_init;
    std::vector<double> parameters;
    parameters.resize(4);
    parameters[0] = -10.0;
    parameters[1] = -10.0;
    parameters[2] = 10.0;
    parameters[3] = 10.0;
    //BVP initialization 
    scalar_init["f"] = Equation_f;
    vector_init["F"] = Equation_F;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping);
    bvp.BVP_init(2,scalar_init,scalarN_init,vector_init, matrix_init,string_init, Equation_RBF);
    GMSolver solver;
    Eigen::VectorXd X0;
    X0.resize(2);
    X0[0] = 5.0; X0[1] = 6.0;
    for(double h = 1; h > 0.005; h = 0.5*h){
        solver.Configure(bvp,parameters,h,78945);
        solver.Solve(X0,0.2);
        std::cout << h << " " << solver.err << " "<< solver.std << std::endl;
    }
    return 0;
}