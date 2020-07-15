#ifndef STENCIL 
#define STENCIL
#include <iostream>
#include <vector>
#include <unistd.h>
#include <stdio.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <fstream>
#include <algorithm>
#include "boundary.hpp"
#include "BVP.hpp"
#define EPSILON 1e-5
class Stencil
{
private:
       //Inverse psi matrix for all the four stencils
	Eigen::MatrixXf ipsi_north, ipsi_south, ipsi_east, ipsi_west;
       //Position of the nodes of the stencils
	std::vector<Eigen::VectorXf> pos_north, pos_south, pos_east, pos_west;
       //2D indexes of the interfaces that originate the stencils 
       std::vector<int> index_north, index_south, index_east, index_west;
       //Cumulatives of the G matrix values associated to the different stencils
       std::vector<float> G_north, G_south, G_east, G_west;
       //Counter of the trayectoies that ends in each stencil
       uint32_t counter_north, counter_south, counter_east, counter_west;
       //Computes the inverse matrix of the psi matrix
       Eigen::MatrixXf Compute_ipsi(std::vector<Eigen::VectorXf> & sten_position, BVP bvp, float c2);
       
public:
       //Parameters of the stencil boundary
       float stencil_parameters[4];
       //Parameters of the boundary of the problem
       float global_parameters[4];
       //Deafault initializer
	Stencil(void);
       /*Float comparison with epsilon tolerance
       float a One of the float numbers to compare
       float b the ther float number
       */
       bool AreSame(float a, float b);
       /*Proper intializer of the Stencil class when the node belongs to one interface
       int dir Direction of the interface its stencil we are computing. 0 for horizontal 1 for vertical.
       bool interior  True if the interface is interior false if it is not
       std::vector<int> current_index Index of the interface
       std::vector<Eigen::VectorXf> & node_position Position of the nodes in the problem
       std::vector<std::vector<int> > & node_index Indexes of the nodes in the problem
       std::map<std::vector<int>, int> & interface_map Map which maps the interfaces 2D index into the interfaces 1D index
       std::vector<int> & number_interfaces Number of interfaces in each direction
       float g_parameters Parameters of the stencil boudnary 
       */
	void Init(int dir,
              bool interior, 
              std::vector<int> current_index,
              std::vector<Eigen::VectorXf> & node_position,
              std::vector<std::vector<int> > & node_index,
              std::map<std::vector<int>, int> & interface_map,
              std::vector<int> & number_interfaces, 
              float *g_parameters);
       /*Proper intializer of the Stencil class when the node belongs to 4 interfaces 
       std::vector<int> interface_east East interface the node belongs to
       std::vector<int> interface_south South interface the node belongs to
       std::vector<int> interface_west West interface the node belongs to
       std::vector<int> interface_north North interface the node belongs to 
       std::vector<int> current_index Index of the interface
       std::vector<Eigen::VectorXf> & node_position Position of the nodes in the problem
       std::vector<std::vector<int> > & node_index Indexes of the nodes in the problem
       std::map<std::vector<int>, int> & interface_map Map which maps the interfaces 2D index into the interfaces 1D index
       std::vector<int> & number_interfaces Number of interfaces in each direction
       float g_parameters Parameters of the stencil boudnary 
       */
       void Init(std::vector<int> interface_east,
              std::vector<int> interface_south,
              std::vector<int> interface_west,
              std::vector<int> interface_north,
              std::vector<Eigen::VectorXf> & node_position,
              std::vector<std::vector<int> > & node_index,
              std::map<std::vector<int>, int> & interface_map,
              std::vector<int> & number_interfaces, 
              float *g_parameters);
       /*Resets the stencil object*/
       void Reset(void);
       /*Computes the inverse psi matrices given
       BVP bvp The boundary value problem 
       float c2 The constan for the mesheless interpolator
       */
       void Compute_ipsi(BVP bvp, float c2);
       /* Updates G matrix 
       Eigen::VectorXf X position where the trayectory ended
       float Y Y FKAC value where the trayectoy ended
       BVP bvp Boundary value problem
       float c2 Constant for the mesheless interpolator*/
       int G_update(Eigen::VectorXf X, float Y, BVP bvp, float c2);
       /* Updates G matrix for the test fucntion
       Eigen::VectorXf X position where the trayectory ended*/
       int G_Test_update(Eigen::VectorXf X);
       /* Test the Meshless interpolator in the stencils
       Eigen::VectorXf X position where the trayectory ended
       BVP bvp Boundary value problem
       float c2 Constant for the mesheless interpolator*/
       float Test_Interpolator(Eigen::VectorXf X, BVP bvp, float c2);
       /* Return G matrix if no nodes are repeated in the stencil
       std::vector<int> & stencil_index The indexes of the nodes that compound the stencil
       std::vector<float> & G The values of the G matrix
       */
       void G_return(std::vector<int> & stencil_index, std::vector<float> & G);
       /* Return G matrix taking care of the repeated nodes of the stencil
       std::vector<int> & stencil_index The indexes of the nodes that compound the stencil
       std::vector<float> & G The values of the G matrix
       int N_tray Total number of trayectories 
       */
       void G_return_withrep(std::vector<int> & stencil_index, std::vector<float> & G, int N_tray);
       /*Returns G function when the meshless interpolator is being tested:
       std::vector<int> & stencil_index Indexes of the stencil nodes
       std::vector<float> & G Values od the G matrix
       */
       void G_Test_return(std::vector<int> & stencil_index, std::vector<float> & G);
       /*Prints the stencil parameters
       int node_index index of the node whose stencil we are trying to compute
       */
       void Print(int node_index);
	
};
#endif