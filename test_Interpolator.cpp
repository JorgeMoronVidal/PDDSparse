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
#include "interface.hpp"
#include "equation.hpp"
#include "rectangle.hpp"

#define N_tray 10
#define REQUEST 1
#define REPLY 2
#define TAG_Gi 10
#define TAG_Gj 11
#define TAG_Gval 12
#define TAG_Bi 20
#define TAG_Bval 21
#define h 0.001f
#define fac 2.0f
#define DEBUG

void interface_metadata(std::vector<Eigen::VectorXf> & start, 
                        std::vector<Eigen::VectorXf> & end,
                        std::vector<std::vector<int> > & inter_indexes,
                        std::vector<std::vector<int> > & node_indexes,
                        std::vector<int> n_inter,
                        std::vector<int> n_node,
                        int & N_node,
                        Eigen::VectorXf SW,
                        Eigen::VectorXf NE,
                        std::map<std::vector<int>, int> & interface_index);

void interface_metadata(std::vector<Eigen::VectorXf> & node_positions,
                        std::vector<std::vector<int> > & node_indexes,
                        std::vector<std::vector<int> > & inter_indexes,
                        std::map<std::vector<int>, int> & interface_map,
                        std::vector<int> n_inter,
                        std::vector<int> n_node,
                        Eigen::VectorXf SW,
                        Eigen::VectorXf NE);
void stencil_metadata(std::vector<std::vector<int> > & stencil,
                      std::vector<Eigen::VectorXf> & stencil_pos,
                      std::vector<int> & stencil_index,
                      std::vector<std::vector<int> > & n_index,
                      int dir,
                      bool interior, 
                      std::vector<int> current_index,
                      std::vector<Eigen::VectorXf> & start, 
                      std::vector<Eigen::VectorXf> & end,
                      std::vector<int> & n_inter,
                      std::vector<int> & n_node,
                      std::map<std::vector<int>, int> & interface_index,  
                      float *params_stencil);

void set_direction_interior(bool & interior,
                            int & dir,
                            int control,
                            std::vector<std::vector <int> > & i_index,
                            std::vector<int> & i_N);

void compute_ipsi(std::vector<std::vector<float> > & ipsi_val,
                  std::vector<Eigen::VectorXf> & sten_position,
                  BVP bvp,
                  float c2);

void fullfill_node_pos(std::vector<Eigen::VectorXf> & node_pos,
                       std::vector<std::vector<int> > i_index,
                       std::vector<std::vector<int> > n_index,
                       std::map<std::vector<int>, int> in_index,
                       std::vector<Eigen::VectorXf> start,
                       std::vector<Eigen::VectorXf> end);

int main(int argc, char *argv[]) {
    /*  -World is the communicator with all the process involved in the program
        -Workers is the comunitar with all the process that are going to be computing stuff
        -Status stores the status of a given process
        -ranks stores the rank of one process we are insterested in
        -numprocs stores the total number of processes
        -myid gives the rank of the current process
        -server stores the rnak of the process in charge of organizing workers work
    */
    MPI_Comm world, workers;
    MPI_Group world_group, worker_group;
    MPI_Status status;
    int ranks[1], numprocs, myid, server, workerid;
    MPI_Init(&argc, &argv);
    world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &numprocs);
    MPI_Comm_rank(world, &myid);
    server = numprocs-1;
    MPI_Comm_group(world, &world_group);
    ranks[0] = server;
    //We exclude a group in world to create workers
    MPI_Group_excl(world_group, 1, ranks, &worker_group);
    //Workers comunicator is created
    MPI_Comm_create(world, worker_group, &workers);
    //We free both comunicators groups
    MPI_Group_free(&worker_group);
    MPI_Group_free(&world_group);
    
    //Start and end stores the first and last node positions of one interface
    std::vector<Eigen::VectorXf> start, end;
    //Indexes of interfaces, nodes and subdomains
    std::vector<std::vector<int>> i_index, n_index, subd_index, stencil;
    //Number of interfaces and nodes in each direction.
    std::vector<int> i_N, n_N;
    std::vector<std::vector<float> > ipsi_val;
    //Initialization of these variables
    i_N.push_back(3); n_N.push_back(10);i_N.push_back(3);n_N.push_back(10);

    //SW and NE points of the domain 
    Eigen::VectorXf SW,NE;
    SW.resize(2); NE.resize(2);
    SW(0) = -1.0f; SW(1) = -1.0f;
    NE(0) = 1.0f;  NE(1) = 1.0f;

    //direction of the interface: 0 if horizontal 1 if vertical
    int dir;

    //0 if the interface is not interior 1 if it is
    bool interior;

    //true if the nodes of the interfaces are distirbuted by chebysehv distribution
    bool chebyshev = false;

    //Parameters of the surfaces
    float sten_p[4], global_p[4];
    sten_p[0] = global_p[0] = SW(0);
    sten_p[1] = global_p[1] = SW(1);
    sten_p[2] = global_p[2] = NE(0);
    sten_p[3] = global_p[3] = NE(1);

    //It gives the index where is stored the information of a given interface in a vector
    std::map<std::vector<int>, int> interface_map;
    std::vector<Eigen::VectorXf> Node_pos;
    interface_metadata(Node_pos, n_index, i_index, interface_map, i_N, n_N, SW, NE);
    //c2 is defined
    float c2 = fac/Node_pos.size();
    //Subdomain is resized
    subd_index.resize(i_index.size());
    //work_control controls the flow of the main program
    int work_control[2] = {0, 0}, aux = 0;
    if(myid == server){
    //I am the server who distributes task
        do {
            aux = work_control[0];
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
                     REQUEST, world, &status);

            work_control[0] = aux;
            MPI_Send(work_control, 2, MPI_INT, 
                    status.MPI_SOURCE, REPLY, world);

            work_control[0]++;

        }while(work_control[0] < (int)i_index.size());
        uint32_t N_node = Node_pos.size();
        Eigen::SparseMatrix<float> G, I, B;
        I.resize(N_node, N_node);
        I.setIdentity();
        G.resize(N_node, N_node);
        B.resize(N_node, 1);
        for(int i = 0; i < server; i++){
            //It ends worker process solving
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
                     REQUEST, world, &status);
            work_control[0] = 0;
            work_control[1] = 1;
            MPI_Send(work_control, 2, MPI_INT, 
                status.MPI_SOURCE, REPLY, world);
        }

    } else {
        //Interface object
        Interface interface;
        //Randon number generator 
        const gsl_rng_type * T;
        gsl_rng * rng;
        unsigned long mySeed;
        gsl_rng_env_setup();
        //Constant seed
        mySeed = (unsigned int) myid +1;
        //NonConstan seed
        /*
        struct timeval tv;
        gettimeofday(&tv,0);
        mySeed = tv.tv_sec * myid + tv.tv_usec;
        */
        T = gsl_rng_default; //Generator setup
        rng = gsl_rng_alloc(T);
        gsl_rng_set(rng, mySeed);

        //Boundary Value problem
        BVP bvp;
        std::map<std::string, pfscalar> scalar_init;
        std::map<std::string, pfvector> vector_init;
        std::map<std::string, pfmatrix> matrix_init;
        std::map<std::string, std::string> string_init;
        //Stencil indexes
        std::vector<int> sten_index;
        //Stencil Positions
        std::vector<Eigen::VectorXf> sten_pos;
        char buffer[100];
        std::ofstream ofile;
        sprintf(buffer,"Output/Debug/Test_Interpolator.txt");
        while(work_control[1] == 0){

            MPI_Send(work_control, 2, MPI_INT, server, REQUEST, world);
            MPI_Comm_rank(workers, &workerid);

            MPI_Recv(work_control, 2, MPI_INT, server, REPLY, world, 
            MPI_STATUS_IGNORE);

            if(work_control[1] == 0){

                set_direction_interior(interior, dir, work_control[0], i_index, i_N);
                ofile.open(buffer, std::ios::app);
                ofile <<"Interface " << i_index[work_control[0]][0] << " " << i_index[work_control[0]][1] << "\n";
                ofile.close();
                std::cout << "Process " << workerid << " is solving itfc "
                << work_control[0] << '\n';
                 
                //BVP initialization 
                scalar_init["f"] = Equation_f;
                scalar_init["c"] = Equation_c;
                scalar_init["u"] = Equation_u;
                scalar_init["g"] = Equation_g;
                vector_init["b"] = Equation_b;
                vector_init["F"] = Equation_F;
                matrix_init["sigma"] = Equation_sigma;
                bvp.Boundary_init(Rectangle2D, Stopping);

                bvp.BVP_init(2,scalar_init, vector_init, matrix_init,string_init, Equation_RBF);
                 //Parameters of the stencil boundary
                sten_p[0] = global_p[0];
                sten_p[1] = global_p[1];
                sten_p[2] = global_p[2];
                sten_p[3] = global_p[3];
                //Stencil 
                Stencil stencil;
                stencil.Init(dir, interior, i_index[work_control[0]], Node_pos, n_index, interface_map, i_N, global_p);
                //Interface initialization
                interface.Init(Node_pos, i_index[work_control[0]],
                               n_index[work_control[0]], dir, subd_index[work_control[0]], interior, 
                               chebyshev,N_tray, h);

                interface.Test_Interpolator(bvp, rng, N_tray, c2, stencil);

                
                #ifdef DEBUG 
                /*
                sleep(work_control[0]);
                interface.Print_Interface();
                std::cout << "With stencil composed by interfaces:\n";
                for(int i = 0; i < (int)stencil.size(); i ++){
                    std::cout << "[" << stencil[i][0] << "," << stencil[i][1] <<"]\t"; 
                }
                std::cout << "And positions\n";
                for(int i = 0; i < (int)sten_pos.size(); i ++){
                    std::cout <<"Node " << sten_index[i] << " with position [" <<
                    sten_pos[i](0) << "," << sten_pos[i](1) <<"]\t";
                } 
                */
                #endif 


            }
        }
        std::cout << "Processor " << workerid << " ended its work \n";
        MPI_Comm_free(&workers);
    }
    MPI_Finalize();
}
    


void interface_metadata(std::vector<Eigen::VectorXf> & start, 
                        std::vector<Eigen::VectorXf> & end,
                        std::vector<std::vector<int> > & inter_indexes,
                        std::vector<std::vector<int> > & node_indexes,
                        std::vector<int> n_inter,
                        std::vector<int> n_node,
                        int & N_node,
                        Eigen::VectorXf SW,
                        Eigen::VectorXf NE,
                        std::map<std::vector<int>, int> & interface_index){

    /*Initial values, long and  and shot increments for each direction*/
   std::vector<Eigen::VectorXf> P_cent, lincrement, sincrement;
   Eigen::VectorXf vaux;


   inter_indexes.resize(n_inter[0]*(n_inter[1]-1)+n_inter[1]*(n_inter[0]-1));
   node_indexes.resize(n_inter[0]*(n_inter[1]-1)+n_inter[1]*(n_inter[0]-1));
   vaux.resize(SW.size());

   //Increments are defined
   //ShortIncrement
   for (int i = 0; i < SW.size() ; i++){
        for(int j = 0; j < vaux.size(); j++){
            if (i == j){

                vaux(j) = (NE(j) - SW(j))/(1.0f* (n_node[j]*n_inter[j]+2));

            } else {

                vaux(j) = 0.0f;

            }
        }
        sincrement.push_back(vaux);
    }
    //LongIncrement
    for (int i = 0; i < SW.size() ; i++){

        for(int j = 0; j < vaux.size(); j++){

            if (i == j){

                vaux(j) = (NE(j) - SW(j) - sincrement[i](j))/(1.0f* (n_inter[j]));

            } else {

                vaux(j) = 0.0f;

            }
        }

        lincrement.push_back(vaux);
    }
    //Initial Point
   
    for (int i = 0; i < SW.size(); i++){
            vaux = SW;
            for (int j = 0; j < SW.size(); j++){

                vaux(j) += sincrement[i](j); 

            }

            if(i > 0){

                vaux(i-1) += 0.5* sincrement[i-1](i-1);

            }

            P_cent.push_back(vaux);
    }

    //This section is not valid for + 2D problems

    P_cent[0] += lincrement[1];
    P_cent[1] += lincrement[0];

    int cen = 0, nodei_cen = 0;
    
    for (int i = 0; i < n_inter[1] - 1; i++){

        for (int j = 0; j < n_inter[0]; j++){

            start.push_back(P_cent[0] + lincrement[0]*j);
            end.push_back(P_cent[0] + lincrement[0]*(j+1) - sincrement[0]);

            inter_indexes[cen].push_back(j);
            inter_indexes[cen].push_back(2*i + 1);
            

            for (int l = 0; l < n_node[0]; l++){

                node_indexes[cen].push_back(nodei_cen);
                nodei_cen ++;
            }

            cen++;

            start.push_back(P_cent[1] + lincrement[1]*j);
            end.push_back(P_cent[1] + lincrement[1]*(j+1) - sincrement[1]);

            inter_indexes[cen].push_back(i);
            inter_indexes[cen].push_back(2*j);

            for (int l = 0; l < n_node[0]; l++){

                node_indexes[cen].push_back(nodei_cen);
                nodei_cen ++;

            }

            cen++;
        }


        P_cent[0] += lincrement[1];
        P_cent[1] += lincrement[0];

    }

    N_node = nodei_cen;

    for(int i = 0; i < (int)inter_indexes.size(); i++){
        interface_index[inter_indexes[i]] = i;
    }


}

void interface_metadata(std::vector<Eigen::VectorXf> & node_positions,
                        std::vector<std::vector<int> > & node_indexes,
                        std::vector<std::vector<int> > & inter_indexes,
                        std::map<std::vector<int>, int> & interface_map,
                        std::vector<int> n_inter,
                        std::vector<int> n_node,
                        Eigen::VectorXf SW,
                        Eigen::VectorXf NE){

    /*Initial values, long and  and shot increments for each direction*/
   std::vector<Eigen::VectorXf> P_cent, lincrement, sincrement, start, end;
   Eigen::VectorXf vaux;


   inter_indexes.resize(n_inter[0]*(n_inter[1]-1)+n_inter[1]*(n_inter[0]-1));
   node_indexes.resize(n_inter[0]*(n_inter[1]-1)+n_inter[1]*(n_inter[0]-1));
   vaux.resize(SW.size());

   //Increments are defined
   //ShortIncrement
   for (int i = 0; i < SW.size() ; i++){
        for(int j = 0; j < vaux.size(); j++){
            if (i == j){

                vaux(j) = (NE(j) - SW(j))/(1.0f* (n_node[j]*n_inter[j]+2));

            } else {

                vaux(j) = 0.0f;

            }
        }
        sincrement.push_back(vaux);
    }
    //LongIncrement
    for (int i = 0; i < SW.size() ; i++){

        for(int j = 0; j < vaux.size(); j++){

            if (i == j){

                vaux(j) = (NE(j) - SW(j) - sincrement[i](j))/(1.0f* (n_inter[j]));

            } else {

                vaux(j) = 0.0f;

            }
        }

        lincrement.push_back(vaux);
    }
    //Initial Point
   
    for (int i = 0; i < SW.size(); i++){
            vaux = SW;
            for (int j = 0; j < SW.size(); j++){

                vaux(j) += sincrement[i](j); 

            }

            if(i > 0){

                vaux(i-1) += 0.5* sincrement[i-1](i-1);

            }

            P_cent.push_back(vaux);
    }

    //This section is not valid for + 2D problems

    P_cent[0] += lincrement[1];
    P_cent[1] += lincrement[0];

    int cen = 0, nodei_cen = 0;
    
    for (int i = 0; i < n_inter[1] - 1; i++){

        for (int j = 0; j < n_inter[0]; j++){

            start.push_back(P_cent[0] + lincrement[0]*j);
            end.push_back(P_cent[0] + lincrement[0]*(j+1) - sincrement[0]);

            inter_indexes[cen].push_back(j);
            inter_indexes[cen].push_back(2*i + 1);
            

            for (int l = 0; l < n_node[0]; l++){

                node_indexes[cen].push_back(nodei_cen);
                nodei_cen ++;
            }

            cen++;

            start.push_back(P_cent[1] + lincrement[1]*j);
            end.push_back(P_cent[1] + lincrement[1]*(j+1) - sincrement[1]);

            inter_indexes[cen].push_back(i);
            inter_indexes[cen].push_back(2*j);

            for (int l = 0; l < n_node[0]; l++){

                node_indexes[cen].push_back(nodei_cen);
                nodei_cen ++;

            }

            cen++;
        }


        P_cent[0] += lincrement[1];
        P_cent[1] += lincrement[0];

    }

    for(int i = 0; i < (int)inter_indexes.size(); i++){
        interface_map[inter_indexes[i]] = i;
    }

    fullfill_node_pos(node_positions, inter_indexes, node_indexes, interface_map,
                      start, end);


}

void stencil_metadata(std::vector<std::vector<int> > & stencil,
                      std::vector<Eigen::VectorXf> & stencil_pos,
                      std::vector<int> & stencil_index,
                      std::vector<std::vector<int> > & n_index,
                      int dir,
                      bool interior, 
                      std::vector<int> current_index,
                      std::vector<Eigen::VectorXf> & start, 
                      std::vector<Eigen::VectorXf> & end,
                      std::vector<int> & n_inter,
                      std::vector<int> & n_node,
                      std::map<std::vector<int>, int> & interface_index,  
                      float *params_stencil){

    std::vector< std::vector<int> > aux;
    std::vector< int > vaux;
    vaux.resize(2);
    if(interior == true){
        switch (dir) {
    
            case 0: 
                //Spline is Horizontal

                //Top spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] + 2;
                params_stencil[3] = end[interface_index[vaux]][1];
                aux.push_back(vaux);
                

                //Bottom spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] - 2;
                params_stencil[1] = end[interface_index[vaux]][1];
                aux.push_back(vaux);

                //Letf-down spline Interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1] - 1;
                params_stencil[0] = end[interface_index[vaux]][0];
                aux.push_back(vaux);

                //Left-up spline interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1] + 1;
                aux.push_back(vaux);

                //Right-down spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] - 1;
                aux.push_back(vaux);

                //Right-up spline Interface
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] + 1;
                params_stencil[2] = end[interface_index[vaux]][0];
                aux.push_back(vaux);

            break;
    
            case 1:
                //Vertical

                //Right spline Interface
                vaux[0] = current_index[0] + 1;
                vaux[1] = current_index[1];
                params_stencil[2] = end[interface_index[vaux]][0];
                aux.push_back(vaux);

                //Left spline Interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1];
                params_stencil[0] = end[interface_index[vaux]][0];
                aux.push_back(vaux);

                //Letf-up spline Interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1] + 1;
                aux.push_back(vaux);

                //Right-up spline interface
                vaux[0] = current_index[0] + 1;
                vaux[1] = current_index[1] + 1;
                params_stencil[3] = end[interface_index[vaux]][1];
                aux.push_back(vaux);

                //Letf-down spline Interface
                vaux[0] = current_index[0] - 1;
                vaux[1] = current_index[1] - 1;
                params_stencil[1] = end[interface_index[vaux]][1];
                aux.push_back(vaux);

                //Right-dowm spline interface
                vaux[0] = current_index[0] + 1;
                vaux[1] = current_index[1] - 1;
                aux.push_back(vaux);

    
            break;
    
            default :
            std::cout <<"Wrong direction parameter. Stencils won't be properly built."<< '\n';
            } 

    }else{
        switch (dir) {
    
            case 0: 
                //Horizontal

                //Top spline Interface
                if(current_index[1] + 2 <= 2*n_inter[1] - 3){
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] + 2;
                params_stencil[3] = end[interface_index[vaux]][1];
                aux.push_back(vaux);
                }

                //Bottom spline Interface
                if(current_index[1]-2 >= 0){
                vaux[0] = current_index[0];
                vaux[1] = current_index[1] - 2;
                params_stencil[1] = end[interface_index[vaux]][1];
                aux.push_back(vaux);
                }

                //Left splines
                if(current_index[0] > 0){
                    //Letf-down spline Interface
                    if(current_index[1]-1 >= 0){
                        vaux[0] = current_index[0] - 1;
                        vaux[1] = current_index[1] - 1;
                        params_stencil[0] = end[interface_index[vaux]][0];
                        aux.push_back(vaux);
                    }
                    //Left-up spline interface
                    if(current_index[1]+1 <= 2*n_inter[1] - 2){
                        vaux[0] = current_index[0] - 1;
                        vaux[1] = current_index[1] + 1;
                        aux.push_back(vaux);
                    }
                }

                //Right splines
                if(current_index[0] < n_inter[0] - 1){
                    //Right-down spline Interface
                    if(current_index[1]-1 >= 0){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] - 1;
                        aux.push_back(vaux);
                    }
                    //Right-up spline interface
                    if(current_index[1]+1 <= 2*n_inter[1] - 2){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] + 1;
                        params_stencil[2] = end[interface_index[vaux]][0];
                        aux.push_back(vaux);
                    }
                }

            break;
    
            case 1:
                //Vertical

                //Right spline Interface
                if(current_index[0] + 1 <= n_inter[0] - 2){
                    vaux[0] = current_index[0] + 1;
                    vaux[1] = current_index[1];
                    params_stencil[2] = end[interface_index[vaux]][0];
                    aux.push_back(vaux);
                }

                //Left spline Interface
                if(current_index[0] - 1 >= 0){
                    vaux[0] = current_index[0] - 1;
                    vaux[1] = current_index[1];
                    params_stencil[0] = end[interface_index[vaux]][0];
                    aux.push_back(vaux);
                }

                //Top splines
                if(current_index[1] + 1 <= 2*n_inter[1] - 3 ){
                    //Letf-up spline Interface
                    if(current_index[0] >= 0){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] + 1;
                        aux.push_back(vaux);
                    }
                    //Right-up spline interface
                    if(current_index[0] + 1 <= n_inter[0] - 1){
                        vaux[0] = current_index[0] + 1;
                        vaux[1] = current_index[1] + 1;
                        params_stencil[3] = end[interface_index[vaux]][1];
                        aux.push_back(vaux);
                    }
                }

                //Bottom splines
                if(current_index[1] - 1 >= 0){
                    //Letf-down spline Interface
                    if(current_index[0] >= 0){
                        vaux[0] = current_index[0];
                        vaux[1] = current_index[1] - 1;
                        params_stencil[1] = end[interface_index[vaux]][1];
                        aux.push_back(vaux);
                    }
                    //Right-dowm spline interface
                    if(current_index[0] + 1 <= n_inter[0] - 1){
                        vaux[0] = current_index[0] + 1;
                        vaux[1] = current_index[1] - 1;
                        aux.push_back(vaux);
                    }
                }

    
            break;
    
            default :
            std::cout <<"Wrong direction parameter. Stencils won't be properly built."<< '\n';
            } 
    
    }

    stencil = aux;
    //Points of the stencil are recovered and stored
    Eigen::VectorXf vxfaux, sincrement;
    stencil_index.clear();
    stencil_pos.clear();

    for(std::vector<std::vector<int> >::iterator item = stencil.begin();
        item != stencil.end();
        item++){

        sincrement = (end[interface_index[*item]]-start[interface_index[*item]])/(n_index[interface_index[*item]].size()-1);

        vxfaux = start[interface_index[*item]];

        for(int i = 0; 
            i < (int) n_index[interface_index[*item]].size(); 
            i++){

            stencil_index.push_back(n_index[interface_index[*item]][i]);
            stencil_pos.push_back(vxfaux);
            vxfaux += sincrement;

            }
        }    
}

void set_direction_interior(bool & interior,
                            int & dir,
                            int control,
                            std::vector<std::vector <int> > & i_index,
                            std::vector<int> & i_N){
    interior = true;
    if(i_index[control][0] == 0) interior = false;
    if (i_index[control][0] == i_N[0] -1) interior = false;

    if(i_index[control][1]%2 == 1){

        //In the case of horizontal interfaces, second component of the index is odd
        dir = 0;

        if(i_index[control][0] == 0) interior = false;
        if (i_index[control][0] == i_N[0] -1) interior = false;
        if(i_index[control][1] == 1) interior = false;
        if (i_index[control][1] >= 2*i_N[1]-3) interior = false;

    } else {

        //In the case of vertical ones, its even
        dir = 1;

        if(i_index[control][0] == 0) interior = false;
        if (i_index[control][0] == i_N[0] -2) interior = false;
        if(i_index[control][1] == 0) interior = false;
        if (i_index[control][1] >= 2*(i_N[1]-1)) interior = false;
    }


}

void compute_ipsi(std::vector<std::vector<float> > & ipsi_val,
                  std::vector<Eigen::VectorXf> & sten_position,
                  BVP bvp,
                  float c2){

    //Auxiliary vectors 
    Eigen::VectorXf vaux, sincrement;
    //Psi Matrix, its inverse and Identity are created 
    Eigen::MatrixXf Psi, iPsi, I;
    //Matrix are resized
    Psi.resize(sten_position.size(), sten_position.size());
    I.resize(sten_position.size(), sten_position.size());
    I.setIdentity();
    //And fullfilled
    for(unsigned int i = 0; i < sten_position.size(); i ++){
        for(unsigned int j = 0; j < sten_position.size(); j ++){
            Psi(i,j) = bvp.rbf.Value(sten_position[i], sten_position[j], c2);
        }
    }
    std::ofstream ofile("Output/Debug/Psi.txt");
    ofile << Psi;
    ofile.close();
    float err, cond;
    iPsi= Psi.inverse();
    Eigen::BiCGSTAB<Eigen::MatrixXf > CGS;
    CGS.compute(Psi);
    iPsi = CGS.solveWithGuess(I,iPsi);
    err = CGS.error();
    cond = Psi.norm()*iPsi.norm();
    #ifdef DEBUG
        std::cout << "Inverse matrix computed with error "<< err <<" and condition number" << cond << std::endl;
    #endif
    Psi.resize(0,0);
    ipsi_val.resize(sten_position.size());
    for(unsigned int i = 0; i < sten_position.size(); i ++){
        ipsi_val[i].resize(sten_position.size());
        for(unsigned int j = 0; j < sten_position.size(); j ++){
            ipsi_val[i][j] = iPsi(i,j);
        }
    }
    ofile.open("Output/Debug/iPsi.txt");
    ofile << iPsi;
    ofile.close();
    I.resize(0,0);
    iPsi.resize(0,0);
}

void fullfill_node_pos(std::vector<Eigen::VectorXf> & node_pos,
                       std::vector<std::vector<int> > i_index,
                       std::vector<std::vector<int> > n_index,
                       std::map<std::vector<int>, int> in_index,
                       std::vector<Eigen::VectorXf> start,
                       std::vector<Eigen::VectorXf> end){
    //Auxiliary vectors 
    Eigen::VectorXf vaux, sincrement;
    //File variable
    std::ofstream ofile;
    node_pos.clear();
    //Loop over the points in the system
    for(std::vector<std::vector<int> >::iterator it = i_index.begin();
    it != i_index.end();
    it++){

        sincrement = (end[in_index[*it]]-start[in_index[*it]])/(n_index[in_index[*it]].size()-1);
        vaux = start[in_index[*it]];
        for(int i = 0; i < (int) n_index[in_index[*it]].size(); i++){

            node_pos.push_back(vaux);
            vaux += sincrement;
        }
    }
}