#include <mpi.h>
#include <iostream>
#include <map>
#include <stdio.h>
#include <fstream>
#include <iterator>
#include "BVP.hpp"
#include "node.hpp"
#include "equation.hpp"
#include "rectangle.hpp"



#define N 101

int main() {
    BVP bvp;

    //We initialize the variables to charge
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, std::string> string_init;

    scalar_init["f"] = Equation_f;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    vector_init["F"] = Equation_F;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping);
    bvp.BVP_init(2,scalar_init, vector_init, matrix_init,string_init);

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    Node node_array[world_size];

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    float params[4] = {-1.0f, -1.0f, 1.0f, 1.0f};
    Eigen::VectorXf Pos, Norm, E_P;
    Pos.resize(2);
    Pos(0) = 0.1 * world_rank;
    Pos(1) = 0.1 * world_rank;
    Norm = Pos;
    E_P = Pos;
    node_array[world_rank].init(Pos, 0.2f, 0.00025f, world_rank*120);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d\n",
           processor_name, world_rank, world_size);
    node_array[world_rank].Solve_FKAK(bvp, params);
    printf("Numerical solution of node %d at position (%f, %f) is: %f +/- %f \n Analytical solution is %f \n", 
    world_rank, node_array[world_rank].x0(0), node_array[world_rank].x0(1), 
    node_array[world_rank].solution, 2*node_array[world_rank].std,bvp.u.Value(node_array[world_rank].x0, Norm));

    // Finalize the MPI environment.
    MPI_Finalize();
    
}