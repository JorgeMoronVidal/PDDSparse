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
#include "stat.hpp"
class Node{

    private:
        //x0 stores the point where the solution is going to be computed 
        Eigen::VectorXf x0;
        Stat statistics;

    public:
    /*  -h stores the value of the time discretization
        -solution stores the value of the solution in x0
        -tolerance stores the tolerable statistic error 
        (relative to the absolute error if solution available)
    */
    float h, solution, tolerance;

    //initialization of the class
    Node(void);
    void init (Eigen::VectorXf X_init, float tol, float discrezitation);
    
    //Obtains node solution
    void Solve(BVP bvp);
};
#endif