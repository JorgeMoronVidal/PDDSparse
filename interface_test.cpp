#include <mpi.h>
#include <iostream>
#include <map>
#include <stdio.h>
#include <fstream>
#include <iterator>
#include "BVP.hpp"
#include "node.hpp"
#include "interface.hpp"
#include "equation.hpp"
#include "rectangle.hpp"

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
    Eigen::VectorXf aux_vec;
    aux_vec.resize(2);
    std::vector<Eigen::VectorXf> start_array, end_array;

    aux_vec(0) = 0.05f; aux_vec(1) = -0.95;
    start_array.push_back(aux_vec);
    aux_vec(0) = 0.05f; aux_vec(1) = -0.05;
    end_array.push_back(aux_vec);

    aux_vec(0) = -0.95f; aux_vec(1) = 0.0;
    start_array.push_back(aux_vec);
    aux_vec(0) = -0.05f; aux_vec(1) = 0.0;
    end_array.push_back(aux_vec);

    aux_vec(0) = 0.05f; aux_vec(1) = 0.0;
    start_array.push_back(aux_vec);
    aux_vec(0) = 0.95f; aux_vec(1) = 0.0;
    end_array.push_back(aux_vec);

    aux_vec(0) = 0.05f; aux_vec(1) = 0.05;
    start_array.push_back(aux_vec);
    aux_vec(0) = 0.05f; aux_vec(1) = 0.95;
    end_array.push_back(aux_vec);

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    Interface interface_array[world_size];

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    std::vector<int> indexes;
    
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d out of %d\n",
           processor_name, world_rank, world_size);
    
    for(int i = 0; i < 9; i++){

        indexes.push_back(9*world_rank + i);
    }

    std::vector<int> subd_i, inter_i;
    int direc;
    switch (world_rank) {
    
        case 0 : 
            subd_i.push_back(0); subd_i.push_back(1);
            inter_i.push_back(0); inter_i.push_back(0);
            direc = 1;

            break;
    
        case 1:
            subd_i.push_back(0); subd_i.push_back(2);
            inter_i.push_back(0); inter_i.push_back(1);
            direc = 0;
            break;

        case 2:
            subd_i.push_back(1); subd_i.push_back(3);
            inter_i.push_back(1); inter_i.push_back(1);
            direc = 0;
            break;

        case 3:
            subd_i.push_back(2); subd_i.push_back(3);
            inter_i.push_back(0); inter_i.push_back(2);
            direc = 1;
            break;
        

        default : 
            std::cout << "Drunk emoji" << '\n';
    
    }
    Eigen::SparseMatrix<float> Pm;
    Pm.resize(9*4, 9*4);
    bool inter = false, chebyshev = false;
    float tol = 0.01, dis = 0.001;
    interface_array[world_rank].Init(start_array[world_rank],
                                     end_array[world_rank],
                                     inter_i,
                                     subd_i,
                                     indexes,
                                     direc,
                                     Pm,
                                     inter, 
                                     chebyshev,
                                     tol,
                                     dis
                                   );
    
    interface_array[0].Print_Interface();
    // Finalize the MPI environment.
    MPI_Finalize();
}