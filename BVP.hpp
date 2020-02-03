#ifndef BoundaryValueProblem
#define BoundaryValueProblem

#include <iostream>
#include <map>
#include <string>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "surface.hpp"

/*
 BVP stores the functions that can be used to compute the solution of 
 such kind of problems using the Feynman-Kak formula.
*/

class BVP
{  
    typedef float (*pfscalar)(Eigen::VectorXf);
    typedef Eigen::VectorXf (*pfvector)(Eigen::VectorXf);
    typedef Eigen::VectorXf (*pfmatrix)(Eigen::MatrixXf);
    typedef float (*pfbound)(float* params, Eigen::VectorXf &, Eigen::VectorXf &);

    private:

        //Default functions which always return 0.0f
        float Default_Scalar(Eigen::VectorXf);
        Eigen::VectorXf Default_Vector(Eigen::VectorXf);

        //LookAtTable stores the LookAtTable of those functions that are not defined analytically
        std::map<std::string, LookUpTable> bvplat;

    public:

        //bvpscalar stores the functions whose image is scalar s.t. c(X), f(X), u(X), varphi(X)
        std::map<std::string, pfscalar> bvpscalar;

        //bvpvector stores the functions whose image is a vector s.t b(X), F(x), g(X)
        std::map<std::string, pfvector> bvpvector;

        //Sigma is the /sigma matrix
        pfmatrix Sigma;

        //All the entries of the maps are set as Default. bvplat remains empty
        BVP(void);

        //The object is properly initialized given the functions or the directorys where the
        // Look At Tables are stored.
        void BVP_init(std::map<std::string, pfscalar>,
                    std::map<std::string, pfvector>,
                    std::map<std::string, pfmatrix>,
                    std::map<std::string, std::string>);
        //Initialization of the surface variable
        void Surf_init(pfbound);
        void Surf_init(std::string);

};
#endif