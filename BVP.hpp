#ifndef BoundaryValueProblem
#define BoundaryValueProblem

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "scalarfunction.hpp"
#include "vectorfunction.hpp"
#include "matrixfunction.hpp"
#include "boundary.hpp"
typedef float (*pfscalar)(Eigen::VectorXf, Eigen::VectorXf);
typedef Eigen::VectorXf (*pfvector)(Eigen::VectorXf, Eigen::VectorXf);
typedef Eigen::MatrixXf (*pfmatrix)(Eigen::VectorXf, Eigen::VectorXf);
typedef float (*pfbound)(float* params, 
                             Eigen::VectorXf & position, 
                             Eigen::VectorXf & exitpoint,
                             Eigen::VectorXf & normal);
typedef bool (*pfstop)(Eigen::VectorXf);
/*
 BVP stores the functions that can be used to compute the solution of 
 such kind of problems using the Feynman-Kak formula.
*/
class BVP
{  

    private:

        std::map<std::string, bool> control;

    public:

        //Functions whose image is scalar s.t. c(X), f(X), u(X), varphi(X)
        ScalarFunction f,c, u, varphi;

        //Functions whose image is a vector s.t b(X), F(x), g(X), psi(X)
        VectorFunction F, mu, b, psi;

        //sigma is the /sigma matrix
        MatrixFunction sigma;
        
        //The boundary is stored in surf
        Boundary boundary;

        //All the entries of the maps are set as Default. bvplat remains empty
        BVP(void);

        //The object is properly initialized given the functions or the directorys where the look up tables are stored.
        void BVP_init(int dim,
                    std::map<std::string, pfscalar> map_fscalar,
                    std::map<std::string, pfvector> map_fvector,
                    std::map<std::string, pfmatrix> sigma,
                    std::map<std::string, std::string> map_lut);

        /*Initialization of the surface variable with analytic boundary
        -boundary is a function s.t. (*float)(float* params, Eigen::VectorXf & position, 
                             Eigen::VectorXf & exitpoint, Eigen::VectorXf & normal)
        -stopf is a function s.t. (*pfstop)(Eigen::MatrixXf)
        */
        void Boundary_init(pfbound boundary, pfstop stopf);
        
        /*Initialization of the surface variable with analytic boundary
        -dim is the dimension of the problem
        -boundary is the directory where the lup of the distance is storedAh per
        -stopf is a function s.t. (*pfstop)(Eigen::MatrixXf)
        */
        void Boundary_init(int dim, std::string boundary , pfstop stopf);

};




#endif