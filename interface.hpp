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
{  
    private:

        //Nodes that are in the interface
        std::vector<Node> node_array;
        //Iterator over all the nodes in the interface
        std::vector<Node>::iterator Interface_it; 
        //True if the domain is interior
        bool interior;
        //True if the interface is horizontal
        bool horizontal;
        //Index Interface, column and row number of the stencil
        std::vector<int> interface_index, stencil_i, stencil_j;
        //Positions of the stencil of the interface
        std::vector<Eigen::VectorXf> stencil_position;
        //Parameters of the surface
        float params[4];
        /*The two coordinates that set the South Est 
        and North west corner of the subdomain*/
        Eigen::VectorXf SE,NW;

    public:
        
        //Default initialization so we can define arrays of Interfaces
        Interface(void);
        //Proper initialization of the class
        void Init (Eigen::VectorXf start, 
                   Eigen::VectorXf end,
                   int number_nodes,
                   std::vector<int> i_index,
                   bool horizontal,
                   bool chebyshev);
        //Stencils are created 
        void Set_Stencil(std::vector<Interface> Sten);
        //The nodes are solved
        Eigen::SparseMatrix<float> Solve (BVP bvp, float * params);

};
#endif