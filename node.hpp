#ifndef NODE
#define NODE
//If we want to debug our code
#define DEBUG
#include <iostream>
#include <map>
#include <random>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "BVP.hpp"

class Node{

    private:
        //x0 stores the point where the solution is going to be computed 
        Eigen::VectorXf x0;

        //UpdateStat upgrades the statistic variables involved in the computation of the solution
        void UpdateStat(float &sum , float &sqsum, float &, int &n, int &sumsteps, int steps);


    public:
    /*
        -h stores the value of the time discretization
        -solution stores the value of the solution in x0
        -tolerance stores the tolerable statistic error 
        (relative to the absolute error if solution available)
        -std is the standard deviation
        -var is the variance
        -pearcoef measures the correlation between the solution
         and ji;
        -
    */
    float h, solution, tolerance, var, std, covar, pearcoef, apl;
    Node(void);
    void init (Eigen::VectorXf X_init, float tol, float discrezitation);
    void Solve(BVP bvp);

};
#endif