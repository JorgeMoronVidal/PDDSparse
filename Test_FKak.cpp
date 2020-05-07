#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <vector>
#include "BVP.hpp"
#include "node.hpp"
#include "equation.hpp"
#include "rectangle.hpp"
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <map>
#include <iterator>
#include <eigen3/Eigen/Dense>

int main(int argc, char *argv[]) {
    std::vector<Eigen::VectorXf> points;
    Eigen::VectorXf vaux;
    vaux.resize(2);
    vaux(0) = 0.0f; vaux(1) = 0.0f;
    points.push_back(vaux);
    vaux(0) = 0.0f; vaux(1) = 0.5f;
    points.push_back(vaux);
    vaux(0) = 0.5f; vaux(1) = 0.5f;
    points.push_back(vaux);
    vaux(0) = 0.0f; vaux(1) = -0.5f;
    points.push_back(vaux);
    vaux(0) = -0.5f; vaux(1) = -0.5f;
    points.push_back(vaux);
    vaux(0) = -0.5f; vaux(1) = 0.5f;
    points.push_back(vaux);
    vaux(0) = 0.5f; vaux(1) = 0.0f;
    points.push_back(vaux);
    vaux(0) = 0.5f; vaux(1) = -0.5f;
    points.push_back(vaux);
     Node node;
    const gsl_rng_type * T;
    gsl_rng * rng;
    unsigned long mySeed;
    gsl_rng_env_setup();
    //Constant seed
    mySeed = (unsigned int) 1;
    T = gsl_rng_default; //Generator setup
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, mySeed);
    //Boundary Value problem
    BVP bvp;
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, std::string> string_init;
    scalar_init["f"] = Equation_f;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    vector_init["F"] = Equation_F;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping);
        bvp.BVP_init(2,scalar_init, vector_init, matrix_init,string_init, Equation_RBF);
    float global_p[4];
    global_p[0] = -1.0f;
    global_p[1] = -1.0f;
    global_p[2] = 1.0f;
    global_p[3] = 1.0f;
    for(int i = 0; i < (int)points.size(); i++){
        std::cout << "Solving [" << points[i](0) 
        << "," << points[i](1) << "]";
        node.init(points[i], 0.25f, 0.001f);
        node.Solve_FKAK(bvp, global_p, rng);
        std::cout << " obtaining " << node.solution << 
        " with analytic solution " << bvp.u.Value(points[i],points[i]) <<
        std::endl;
    }

    std::cout << "here";
}
