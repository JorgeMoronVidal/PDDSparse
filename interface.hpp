#ifndef INTERFACE
#define INTERFACE

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include "node.hpp"
#include "stencil.hpp"
#define DEBUG


class Interface
{   friend class Node;
    private:

        //Nodes that are in the interface
        std::vector<Node> node_array;
        //Iterator over all the nodes in the interface
        std::vector<Node>::iterator Interface_it; 
        //True if the domain is interior
        bool interior;
        //0 if horizontal 1 if vertical
        int direction;
        //Index Interface, subdomain index;
        std::vector<int> interface_index, subdomain_index;
        //Parameters of the surface
        float params[4];
        /*The two coordinates that set the Sstd::vector<char> v;outh Est 
        and North west corner of the subdomain*/
        Eigen::VectorXf SE,NW;

    public:
        //Default initialization so we can define arrays of Interfaces
        Interface(void);
        //Proper initialization of the class

        void Init (Eigen::VectorXf start, 
                   Eigen::VectorXf end,
                   std::vector<int>  i_index,
                   std::vector<int>  n_index,
                   int direc,
                   std::vector< int > subd_index,
                   bool inter,
                   bool chebyshev,
                   int N_Steps,
                   float discretization);

        void Init (std::vector<Eigen::VectorXf> node_pos,
                   std::vector<int>  i_index,
                   std::vector<int>  n_index,
                   int direc,
                   std::vector< int > subd_index,
                   bool inter,
                   bool chebyshev,
                   int N_Steps,
                   float discretization);

        //The nodes are solved
        void Solve(BVP bvp,
                  gsl_rng *rng,
                  float * parameters_stencil,
                  float * parameters_surface,
                  int N_tray,
                  float c2,
                  std::vector< std::vector<float> > & iPsi,
                  std::vector< Eigen::VectorXf > & stencil_position,
                  std::vector< int > & stencil_index,
                  std::vector<float> & G,
                  std::vector<float> & B,
                  std::vector<int> & G_j,
                  std::vector<int> & G_i,
                  std::vector<int> & B_i);

        void Solve(BVP bvp,
                  gsl_rng *rng,
                  int N_tray,
                  float c2,
                  Stencil stencil,
                  std::vector<float> & G,
                  std::vector<float> & B,
                  std::vector<int> & G_j,
                  std::vector<int> & G_i,
                  std::vector<int> & B_i);
        void Test_Interpolator(BVP bvp,
                  gsl_rng *rng,
                  int N_tray,
                  float c2,
                  Stencil stencil);
        void Test_G(BVP bvp,
                  gsl_rng *rng,
                  int N_tray,
                  float c2,
                  Stencil stencil,
                  std::vector<float> & G,
                  std::vector<float> & B,
                  std::vector<int> & G_j,
                  std::vector<int> & G_i,
                  std::vector<int> & B_i);
        //Prints on screen interface atributes 
        void Print_Interface(void);

};
#endif