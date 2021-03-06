#ifndef NODE
#define NODE
#include <iostream>
#include <fstream>
#include <map>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <sys/time.h>
#include <stdio.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include "lookuptable.hpp"
#include "BVP.hpp"
#include "rectangle.hpp"
#include "stencil.hpp"
class Node{

  private:
    //x0 stores the point where the solution is going to be computed 
    Eigen::VectorXf X, N, E_P;
    float Y,Z,xi;
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
    float h, t, sqh, solution, tolerance,   
          var, std, covar,pearson_c;
    //Solved indicates if the solution for the node has been computed
    bool solved;
    /*
    -Subdomain.
    -i_node stores the index of the node
    -N stores the number of MC realizations the algorithm has to do 
    -i_interface stores the index of the interface which stores the node0
    */
    int i_node, N_trayectories;
    std::vector<int> i_interface;
    std::vector<int> subdomains;

    //Initialization of the class
    Node(void);
    void init (Eigen::VectorXf X_init, 
               float tol, 
               float discrezitation);

    void init (Eigen::VectorXf X_init, 
               float discretization,
               int N_tray, 
               int node_index,
               std::vector<int> interface_index,
               std::vector<int> subdomain);

    //G using Feynmann Kak Formula
    void Solve_FKAK(BVP bvp, float * params, gsl_rng *rng);
    //Updates the statistical variables
    void Update_Stat(float sol_0, float xi, float & summ, float & mean,
                 float & sqsumm, float & summ_0, 
                 float & sqsumm_0, float & summ_xi, float & sqsumm_xi, 
                 float & var_xi, float & crossumm, int & counter, int & counterN, 
                 int & summN);

    /*
    Elliptic solver for the node
    BVP bvp The boundary Value problem
    gsl_rng *rng GSL random number generator variable
    int N_tray Number of trayectories per node
    float c2 Constant for the meshless interpolator
    Stencil stencil Node stencil
    std::vector<int> & G_j G matrix column indexes
    std::vector<float> & G G matrix values
    float & B B vector value for the node
    */
    void  Solve_PDDSparse(BVP bvp,
                       gsl_rng *rng,
                       int N_tray,
                       float c2,
                       Stencil stencil,
                       std::vector<int> & G_j,
                       std::vector<float> & G,
                       float & B);
    /*
    Parabolic solver for the node
    float T is the time for the trayectories to start
    BVP bvp The boundary Value problem
    gsl_rng *rng GSL random number generator variable
    int N_tray Number of trayectories per node
    float c2 Constant for the meshless interpolator
    Stencil stencil Node stencil
    std::vector<int> & G_j G matrix column indexes
    std::vector<float> & G G matrix values
    float & B B vector value for the node
    */
    void  Solve_PDDSparse(float T,
                       BVP bvp,
                       gsl_rng *rng,
                       int N_tray,
                       float c2,
                       Stencil stencil,
                       std::vector<int> & G_j,
                       std::vector<float> & G,
                       float & B);
    /*
    Test of the mesheless interpolator for a given node (Old version)
    */
    void Test_Interpolator(BVP bvp,
                       gsl_rng *rng, 
                       float * parameters_stencil,
                       float * parameters_surface,
                       int N_tray,
                       float c2,
                       Stencil stencil);
    /*
    Test of the G storage algorithm for a given node (Old version)*/
    void Test_G(BVP bvp,
                       gsl_rng *rng, 
                       float * parameters_stencil,
                       float * parameters_surface,
                       int N_tray,
                       float c2,
                       Stencil stencil,
                       std::vector<int> & stencil_index,
                       std::vector<float> & G,
                       float & B);
    //Draws a trayectory of the brownian path
    void Draw_trayectory(void);

    //Prints the result in an output file
    void Print_Result(void);
};
#endif