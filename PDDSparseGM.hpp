//mpiexec --mca orte_base_help_aggregate 0  -np 4 ./main
#ifndef PDDSPARSEGM
#define PDDSPARSEGM
//Fraction of the interface the stencil is elongated
#define STEN_ELONG 0.1
//Minimum value of N
#define N_min 500
//Minimum value of h
#define h_min 1E-05
//MPI TAGS
#define REQUEST_NODE 1
#define REPLY_NODE 2
#define NODE_POSITION 3
#define PARAMETERS 4
#define SIZES 5
#define IND_NORTH 6
#define IND_SOUTH 7
#define IND_EAST 8
#define IND_WEST 9
#define X_NORTH 10
#define X_SOUTH 12
#define X_EAST 14
#define X_WEST 16
#define Y_NORTH 11
#define Y_SOUTH 13
#define Y_EAST 15
#define Y_WEST 17
#define INTER_NUMBER 18
#define INTER_LABELS 19
#define TAG_Gi 20
#define TAG_Gj 21
#define TAG_Gval 22
#define TAG_Gvar 23
#define TAG_Bi 30
#define TAG_Bval 31
#define TAG_Bvar 32
#define TODO_JOB_UINT  101
#define TODO_JOB_INT  102
#define TODO_JOB_DOUBLE 103
#define DONE_JOB_UINT 111
#define DONE_JOB_INT 112
#define DONE_JOB_DOUBLE 113
//Interfaces
#define INTERSECTIONS_YES
#include "GMSolver.hpp"
#include "interface.hpp"
#include "stencil.hpp"
#include "subdomain.hpp"
#include<eigen3/Eigen/SparseLU>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <fstream>
#include <forward_list>
//enum direction {North, South, East, West};
typedef Eigen::Triplet<double,int> T;
class PDDSJob {
    public:
    int index[1];
    unsigned int N[2];
    double h[1];
    PDDSJob(void);
    void Send_To_Worker(MPI_Status & status, MPI_Comm & world);
    void Recieve_From_Server(int server, MPI_Comm & world);
    //void Send_To_Server(int server, MPI_Comm & world);
    //void Recieve_From_Worker(MPI_Status status, MPI_Comm & world);
};

class PDDSparseGM{
    friend class GMSolver;
    private:
        /*
        MPI-related variables
        */
        MPI_Comm world, workers;
        MPI_Group world_group, worker_group;
        MPI_Status status;
        int ranks[1],work_control[2],numprocs, myid, server, workerid;
        /*
        Interfaces of the problem
        */
        std::vector<Interface> interfaces;
        /*
        Subdomains of the problem
        */
        std::vector<Subdomain> subdomains;
        /*-maps each interface 2D index in its 1D index*/
        std::map<std::vector<int>, int> interface_map;
        /*
         -Number of interfaces (And subdomains) per direction.
         -Nuber of nodes per interface.
        */
        std::vector<int> iN,nN;
        /*File were the parameters of the problem are written*/
        std::string configuration_file;
        /*-h0 is the value of time discretization for the lowest level.
          -T_start is the time value where we want to compute dhe solution
          -eps represents the rse objective of the stimator
          -parameters stores the boundary parameters
          -factor that multiplies the internodal distance
          -constant of the RBF interpolator 
          */
        double h0, T_start, eps, fac, c2;
        std::vector<int> h_vec;
        std::vector<double> parameters;
        /*-SW is the south west point of the domain.
          -NE is the north east point of the domain.
        */
       Eigen::VectorXd SW, NE;
       /*G and B storage vector*/
       std::vector<double> G, B, G_var, B_var;
       std::vector<int>  G_j, G_i, B_i;
       /*Triplet's vector*/
       std::vector<T> T_vec_G, T_vec_Gvar, T_vec_B, T_vec_Bvar;
        /*
          -N is the initial number of trayectories 
        */
       int N, N_job, nNodes;
       std::vector<int> N_vec;
       /*Stencil of the node*/
       Stencil node_stencil;
       /*Position of the Node being solver*/
       Eigen::VectorXd position;
       /*Filename of the Flux DEBUG File*/
       char debug_fname[256];
       /*Starts the MPI values and structures*/
       void MPI_Configuration(int argc, char *argv[]);
       /*Returns the interfaces the node belongst to*/
       std::vector<std::vector<int> > Get_Interfaces(int index);
       /*Returns the stencils labels*/
       std::map<direction, std::vector<std::vector <int> > > Labels_Stencil(int index);
       /*Send node information*/
       void Send_Node(PDDSJob job);
       /*Send node information*/
       void Send_Node_Loop(PDDSJob job);
       /*Receive node information*/
       bool Receive_Node(PDDSJob & job);
       /*Sends the stencil data to the worker*/
       void Send_Stencil_Data(int index);
       void Send_Stencil_Data_Loop(int index);
       /*Receives the stencil data from the server*/
       Stencil Recieve_Stencil_Data(void);
       Stencil Recieve_Stencil_Data_Loop(void);
       /*Send G matrix and B Matrix*/
       void Send_G_B(void);
       /*Receive G matrix and B Matrix*/
       void Receive_G_B(void);
       /*Having G and B as triplets, computes the solution for the PDDS problem*/
       void Compute_Solution(BVP bvp);
       /*Position of a given node*/
       Eigen::VectorXd Node_Position(int node_index);
       /*Prints metadata file*/
       void Print_metadata(void);
       /*Prints solution file*/
       void Print_solution(BVP bvp);
       /*Prints Domains files*/
       void Print_subdomain(BVP bvp);
       /*Reading subdomain files routines*/
       std::vector<double> Read_File(char fname[256]);
    public:
       /*Reads the problem information from a file.*/
       void ReadFromFile(std::string file);
       /*Interface and node structure intialization*/
       void Fullfill_interfaces(void);
       /*Interface and node structure intialization*/
       void Fullfill_subdomains(BVP bvp);
       /*Initialization of the class*/
       PDDSparseGM(int argc, char *argv[]);
       /*Initialization of the class with a configuration file*/
       PDDSparseGM(int argc, char *argv[], std::string file);
       /*Solves the problem given by the configuration file*/
       void Solve(BVP bvp, std::string file);
       /*Solves the problem given by the configuration file*/
       void Solve(BVP bvp); 
       /*Prints the problem parameters on screen*/
       void Print_Problem(void); 
       /*Prints the content of the Interface vector on screen*/
       void Print_Interface(void);
       /*Reads solution file*/
       void Read_Solution(void);
       /*Solves subdomains*/
       void Solve_Subdomains(BVP bvp);
       /*Solves with numerical Variance Reduction after a warm-up fase*/
       void Solve_NumVR(BVP bvp, std::vector<double> h_vec, std::vector<unsigned int> N_vec);
};
#endif