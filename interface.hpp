#ifndef INTERFACE
#define INTERFACE

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <eigen3/Eigen/Core>
#include "node.hpp"
#include "rectangle.hpp"


class Interface
{  
    private:
        //Nodes that are in the interface
        std::vector<Node> node_array;
        //Iterator over all the nodes in the interface
        std::vector<Node>::iterator Interface_it; 
        bool interior;
        float params[4];
        //The two coordinates that set the South Est and North west corner of the 
        Eigen::VectorXf SE,NW;

    public:

        Interface(void);
        void Init (Eigen::VectorXf start, 
                   Eigen::VectorXf end,
                   int number_nodes);
        void Solve (void);

};
#endif