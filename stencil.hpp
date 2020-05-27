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
	Eigen::MatrixXf ipsi_north, ipsi_south, ipsi_east, ipsi_west;
	std::vector<Eigen::VectorXf> pos_north, pos_south, pos_east, pos_west;
       std::vector<int> index_north, index_south, index_east, index_west;
       std::vector<float> G_north, G_south, G_east, G_west;
       uint32_t counter_north, counter_south, counter_east, counter_west;
       Eigen::MatrixXf Compute_ipsi(std::vector<Eigen::VectorXf> & sten_position, BVP bvp, float c2);
       
public:
       float stencil_parameters[4];
       float global_parameters[4];
	Stencil(void);
       bool AreSame(float a, float b);
	void Init(int dir,
              bool interior, 
              std::vector<int> current_index,
              std::vector<Eigen::VectorXf> & node_position,
              std::vector<std::vector<int> > & node_index,
              std::map<std::vector<int>, int> & interface_map,
              std::vector<int> & number_interfaces, 
              float *g_parameters);
       void Init(std::vector<int> interface_east,
              std::vector<int> interface_south,
              std::vector<int> interface_west,
              std::vector<int> interface_north,
              std::vector<Eigen::VectorXf> & node_position,
              std::vector<std::vector<int> > & node_index,
              std::map<std::vector<int>, int> & interface_map,
              std::vector<int> & number_interfaces, 
              float *g_parameters);
       void Reset(void);
       void Compute_ipsi(BVP bvp, float c2);
       int G_update(Eigen::VectorXf X, float Y, BVP bvp, float c2);
       int G_Test_update(Eigen::VectorXf X);
       float Test_Interpolator(Eigen::VectorXf X, BVP bvp, float c2);
       void G_return(std::vector<int> & stencil_index, std::vector<float> & G);
       void G_return_withrep(std::vector<int> & stencil_index, std::vector<float> & G, int N_tray);
       void G_Test_return(std::vector<int> & stencil_index, std::vector<float> & G);
       void Print(int node_index);
	
};
#endif