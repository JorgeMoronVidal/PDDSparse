#include <iostream>
#include "BVP.hpp"
#include "edepaco.hpp"
#include "sphere.hpp"

int main() {
    BVP bvp;

    //We initialize the variables to charge
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, std::string> string_init;

    scalar_init["f"] = fpaco;
    scalar_init["c"] = cpaco;
    scalar_init["u"] = upaco;
    scalar_init["g"] = gpaco;
    scalar_init["varphi"] = varpsipaco;
    scalar_init["psi"] = psipaco;
    vector_init["b"] = bpaco;
    vector_init["F"] = Fpaco;
    matrix_init["sigma"] = sigmapaco;
    bvp.Boundary_init(Sphere, Stopping);
    bvp.BVP_init(2,scalar_init, vector_init, matrix_init,string_init);
   
    return 0;
}