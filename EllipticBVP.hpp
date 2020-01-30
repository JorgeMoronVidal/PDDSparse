#ifndef EBVP
#define EBVP

#include <iostream>
#include <map>
#include <string>
#include <eigen3/Eigen/Core>

class EllipticBVP
{   

    typedef float (*pfscalar)(Eigen::VectorXf);
    typedef Eigen::VectorXf (*pfvector)(Eigen::VectorXf);

    private:

        

    public:

        std::map<std::string, pfscalar> bvpscalar;
        std::map<std::string, pfvector> bvpvector;
        EllipticBVP(void);

};
#endif