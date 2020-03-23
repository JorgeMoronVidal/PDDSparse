#ifndef INTERFACE
#define INTERFACE

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Sparse>
#include "node.hpp"
#include "rectangle.hpp"


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
        //Index Interface, column and row number of the stencil
        std::vector<int> interface_index, subdomain_index, stencil_j;
        //Positions of the stencil of the interface
        std::vector<Eigen::VectorXf> stencil_position;
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
                   std::vector<int> i_index,
                   std::vector<int> subd_index,
                   std::vector<int> n_index,
                   int direc,
                   Eigen::SparseMatrix<float> Psim,
                   bool inter,
                   bool chebyshev,
                   float tol,
                   float discretization);
        //Stencils are created 
        void Set_Stencil(std::vector<Interface> Sten, float *parameters);
        //The nodes are solved
        Eigen::SparseMatrix<float> Solve (BVP bvp);
        //Prints on screen interface atributes 
        void Print_Interface(void);

};
#endif