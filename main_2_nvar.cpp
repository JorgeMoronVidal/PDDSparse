#include "PDDSparseGM_old.hpp"
#include "poisson_3.hpp"
#include "rectangle.hpp"
//mpiexec -np 4 xterm -e gdb ./main -ex run
int main(int argc, char *argv[]){
    std::string config("configuration.txt");
    PDDSparseGM PDDS(argc,argv,config);
    BVP bvp;
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, pfscalarN> scalarN_init;
    std::map<std::string, std::string> string_init;
    //BVP initialization 
    scalar_init["f"] = Equation_f;
    //vector_init["F"] = Equation_F;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping);
    bvp.BVP_init(2,scalar_init,scalarN_init,vector_init, matrix_init,string_init, Equation_RBF);
    PDDS.Solve(bvp);
    return 0;
}