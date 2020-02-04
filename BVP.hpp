#ifndef BoundaryValueProblem
#define BoundaryValueProblem

#include <iostream>
#include <map>
#include <string>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "boundary.hpp"

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
    typedef bool (*pfstop)(Eigen::MatrixXf);

    private:

        //Default functions which always return 0.0f
        float Default_Scalar(Eigen::VectorXf);
        Eigen::VectorXf Default_Vector(Eigen::VectorXf);
        Eigen::MatrixXf Default_Matrix(Eigen::VectorXf);
        //LookAtTable stores the LookAtTable of those functions that are not defined analytically
        std::map<std::string, LookUpTable> bvplat;

    public:

        //Functions whose image is scalar s.t. c(X), f(X), u(X), varphi(X)
        pfscalar f,c, u, varphi;

        //Functions whose image is a vector s.t b(X), F(x), g(X), psi(X)
        pfvector F, mu, b, psi;

        //sigma is the /sigma matrix
        pfmatrix sigma;
        
        //The boundary is stored in surf
        Boundary bound;

        //The type of BC's is stored in stop
        pfstop stop;

        //All the entries of the maps are set as Default. bvplat remains empty
        BVP(void);

        //The object is properly initialized given the functions or the directorys where the look up tables are stored.
        void BVP_init(std::map<std::string, pfscalar> map_fscalar,
                    std::map<std::string, pfvector> map_fvector,
                    std::map<std::string, pfmatrix> map_fmatrix,
                    std::map<std::string, std::string> map_lut);

        //Initialization of the surface variable
        void Surf_init(pfbound boundary, pfstop stopf);
        void Surf_init(std::string boundary , pfstop stopf);

};
#endif