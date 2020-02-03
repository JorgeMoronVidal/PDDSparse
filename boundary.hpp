#ifndef BOUND
#define BOUND

#include <iostream>
#include <map>
#include <string>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"

class Boundary{

    typedef float (*pfbound)(float* params, Eigen::VectorXf &, Eigen::VectorXf &);

    private:

        pfbound Distance_Analytic;
        LookUpTable Distance_Numeric;
        bool analytic;

    public:

        Boundary(void);
        void _init_(pfbound);
        void _init_(std::string);
        float dist(float* params, Eigen::VectorXf &, Eigen::VectorXf &);

};
#endif