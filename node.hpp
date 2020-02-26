#ifndef NODE
#define NODE
//If we want to debug our code
#define DEBUG
#include <iostream>
#include <map>
#include <random>
#include <thread>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "BVP.hpp"
class Node{

    private:
        //x0 stores the point where the solution is going to be computed 
        Eigen::VectorXf x0,X, N, E_P;
        float Y,Z;
        unsigned int seed;
        /*Those  variables are used by the normal random number generator 
        with the polar Marsaglia method*/
        float spare, u, v, s;
	    bool generate;  
        const unsigned int max_mt = std::mt19937::max();
    public:
    /*  -h stores the value of the time discretization
        -solution stores the value of the solution in x0
        -tolerance stores the tolerable statistic error 
        (relative to the absolute error if solution available)
        -var is the variance of the result
        -std is the standard deviation of the result
        -covar is the covariance between xi and the solution 
    */
    float h, solution, tolerance,   
          var, std, covar;
    //Solved indicates if the solution for the node has been computed
    bool solved;
    //Initialization of the class
    Node(void);
    void init (Eigen::VectorXf X_init, 
               float tol, 
               float discrezitation, 
               unsigned int random_seed);
    
    //Obtains node solution using Feynmann Kak Formula
    void Solve_FKAK(BVP bvp);

    //Obtains node solution for the Meshless Algorithm
    void Solve_PDDSparse(BVP bvp);

    //Updates the statistical variables
    void Update_Stat(float sol_0, float xi, float & summ, float & mean,
                 float & sqsumm, float & var, float & std, float & summ_0, 
                 float & sqsumm_0, float & summ_xi, float & sqsumm_xi, 
                 float & var_xi, float & crossumm, int & counter, int counterN, 
                 int & summN);

    //Draws a trayectory of the brownian path
    void Draw_trayectory(void);

    //Prints the result in an output file
    void Print_Result(void);

    //Turns a uniform random number distribution into a normal one with avg = 0 and var = 1
    float Random_Normal(int random_int);
};
#endif