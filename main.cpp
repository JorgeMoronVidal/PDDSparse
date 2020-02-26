#include <mpi.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <iterator>
#include "BVP.hpp"
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
    
    float params[4] = {-1.0f, -1.0f, 1.0f, 1.0f};
    Eigen::VectorXf Pos, Norm, E_P;
    Pos.resize(2);
    Pos(0) = 0.1;
    Pos(1) = 0.1;
    Norm = Pos;
    E_P = Pos;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
   printf("Hello world from processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    std::cout << bvp.boundary.Dist(params, Pos, E_P,Norm) << "\n";
    // Finalize the MPI environment.
    MPI_Finalize();
    
}