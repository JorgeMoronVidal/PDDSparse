#ifndef SOLVERPDDSPARSE
#define SOLVERPDDSPARSE
#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <map>
#include <stdio.h>
#include <fstream>
#include <iterator>
#include <eigen3/Eigen/Dense>
#include "stencil.hpp"
#include "BVP.hpp"
#include "node.hpp"
#include "equation.hpp"
#include "rectangle.hpp"

#define N_tray 10000
#define REQUEST 1
#define REPLY 2
#define TAG_Gi 10
#define TAG_Gj 11
#define TAG_Gval 12
#define TAG_Bi 20
#define TAG_Bval 21

class SolverPDDS{

  private:
    /*
    - h is the time discretization of the numerical stimator for the FKAK formula
    - fac is the factor that multiplies internodal separation in order to compute the c
      of the meshless interpolator of the interfaces
    */
    float h, fac;
    /*
    - workers is the comunicator of those processes that compute the 
    */
    MPI_Comm world, workers;
    MPI_Group world_group, worker_group;
    MPI_Status status;
    int ranks[1], numprocs, myid, server, workerid;
    /*
    - interface_indexes stores the indexes of the interfaces
    - node_indexes stores the indexes of the nodes
    - node_interface associates each node with its interfaces
    - subd_index stores the indexes of the subdomains
    - i_N stores the number of interfaces in each direction index
      0 for horizontal 1 for vertical
    - n_N stores the number of nodes in each direction index
      0 for horizontal 1 for vertical
    */
    std::vector<std::vector<int>> interface_indexes, node_indexes, 
                                  node_interface, subd_index,
                                  i_N, n_N;
    void interface_metadata(std::vector<Eigen::VectorXf> & node_positions,
                        std::vector<std::vector<int> > & node_indexes,
                        std::vector<std::vector<int> > & inter_indexes,
                        std::map<std::vector<int>, int> & interface_map,
                        std::vector<int> n_inter,
                        std::vector<int> n_node,
                        Eigen::VectorXf SW,
                        Eigen::VectorXf NE);

    unsigned int fullfill_node_index(std::vector <std::vector<int> >& node_indexes, 
                         std::vector<int> n_inter, 
                         std::vector<int> n_node, 
                         std::map<std::vector<int>, int> interface_map);

    void fullfill_node_interface(unsigned int N_node,
                             std::vector <std::vector<int> >& node_interface, 
                             std::vector <std::vector<int> > node_indexes);

    void fullfill_node_pos(std::vector<Eigen::VectorXf> & node_pos,
                       std::vector<std::vector<int> > i_index,
                       std::vector<std::vector<int> > n_index,
                       std::map<std::vector<int>, int> in_index,
                       std::vector<Eigen::VectorXf> start,
                       std::vector<Eigen::VectorXf> end,
                       unsigned int N_nodes);

    void set_direction_interior(bool & interior,
                            int & dir,
                            int control,
                            std::vector<std::vector <int> > & i_index,
                            std::vector<int> & i_N);

    void Print_interface(int node,
                    std::vector< std::vector<int> > node_interface,
                    std::vector< std::vector<int> > node_index,
                    std::vector< Eigen::VectorXf > node_pos);
  public:

    SolverPDDS(int argc, char *argv[]);

};
#endif