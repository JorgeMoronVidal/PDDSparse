#ifndef NODE
#define NODE
#define EIGEN_SPARSEVECTOR_PLUGIN
//If we want to debug our code
#define DEBUG
#include <iostream>
#include <map>
#include <random>
#include <thread>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include "lookuptable.hpp"
#include "BVP.hpp"

class Node{

    private:
        //x0 stores the point where the solution is going to be computed 
        Eigen::VectorXf X, N, E_P;
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
    Eigen::VectorXf x0;
    float h, solution, tolerance,   
          var, std, covar,pearson_c;
    //Solved indicates if the solution for the node has been computed
    bool solved;
    /*
    -Subdomain.
    -i_node stores the index of the node
    -i_interface stores the index of the interface which stores the node0
    */
    int i_node;
    std::vector<int> i_interface;
    std::vector<int> subdomains;


    //i_node row of the inverted Psi matrix
    Eigen::SparseMatrix<float> psim_i;

    //Initialization of the class
    Node(void);
    void init (Eigen::VectorXf X_init, 
               float tol, 
               float discrezitation, 
               unsigned int random_seed);

    void init (Eigen::VectorXf X_init, 
               Eigen::MatrixXf H_mtrx,
               float tol, 
               float discretization, 
               unsigned int random_seed,
               int node_index,
               std::vector<int> interface_index,
               std::vector<int> subdomains,
               Eigen::SparseMatrix<float> psi_m );

    //G using Feynmann Kak Formula
    void Solve_FKAK(BVP bvp, float * params);
    //Updates the statistical variables
    void Update_Stat(float sol_0, float xi, float & summ, float & mean,
                 float & sqsumm, float & summ_0, 
                 float & sqsumm_0, float & summ_xi, float & sqsumm_xi, 
                 float & var_xi, float & crossumm, int & counter, int & counterN, 
                 int & summN);

    //Obtains i_node row of the G matrix for the Meshless Algorithm
    Eigen::SparseMatrix<float> Solve_PDDSparse(BVP bvp, float * params);

    //It returns a row of the Hij matrix
    Eigen::SparseMatrix<float>  H_Update(std::vector<Eigen::VectorXf>,
                                         std::vector<int>);

    //Draws a trayectory of the brownian path
    void Draw_trayectory(void);

    //Prints the result in an output file
    void Print_Result(void);

    //Turns a uniform random number distribution into a normal one with avg = 0 and var = 1
    float Random_Normal(std::mt19937 & random_int);
};
#endif