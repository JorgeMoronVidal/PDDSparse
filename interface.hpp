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
        /*Proper initialization of the class:
        std::vector<Eigen::VectorXf> node_pos Vector with the position of every node,
        std::vector<int>  i_index Vector with the 2D index of the interface.
        std::vector<int>  n_index Vector with the indexes of every node.
        int direc 1 if the einterface is vertical 0 if horizontal
        std::vector< int > subd_index 2D index of the subdomain.
        bool inter True if the interface is interior false if it is not.
        bool chebyshev True if we wanna the nodes to follow a Chevyshev distribution flase if not.
        int N_Steps Number of steps per node of the interface
        float discretization discretization of the FKAK algorithm)
        */
        void Init (std::vector<Eigen::VectorXf> node_pos,
                   std::vector<int>  i_index,
                   std::vector<int>  n_index,
                   int direc,
                   std::vector< int > subd_index,
                   bool inter,
                   bool chebyshev,
                   int N_Steps,
                   float discretization);

        /*
        Solves all the nodes of the interface

        BVP bvp Boundary  value problem
        gsl_rng *rng GSL random number generator
        int N_tray Number of trayectories per node.
        float c2 constant for the Meshelss interpolator
        Stencil stencil stencil of the interface
        std::vector<float> & G G Matrix values
        std::vector<float> & B B Matrix values 
        std::vector<int> & G_j G matrix column index
        std::vector<int> & G_i G matrix  row index
        std::vector<int> & B_i B index
        */
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
        /*
        Test the meshless interpolator
        
        BVP bvp Boundary  value problem
        gsl_rng *rng GSL random number generator
        int N_tray Number of trayectories per node.
        float c2 constant for the Meshelss interpolator
        Stencil stencil stencil of the interface
        */
        void Test_Interpolator(BVP bvp,
                  gsl_rng *rng,
                  int N_tray,
                  float c2,
                  Stencil stencil);
        /*
        Test the G matrix construction

        BVP bvp Boundary  value problem
        gsl_rng *rng GSL random number generator
        int N_tray Number of trayectories per node.
        float c2 constant for the Meshelss interpolator
        Stencil stencil stencil of the interface
        std::vector<float> & G G Matrix values
        std::vector<float> & B B Matrix values 
        std::vector<int> & G_j G matrix column index
        std::vector<int> & G_i G matrix  row index
        std::vector<int> & B_i B index
        */
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