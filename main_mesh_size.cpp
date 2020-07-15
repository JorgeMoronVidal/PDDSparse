#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <map>
#include <stdlib.h> 
#include <fstream>
#include <iterator>
#include <eigen3/Eigen/Dense>
#include <time.h>
#include "stencil.hpp"
#include "BVP.hpp"
#include "node.hpp"
#include "sara_parabolic.hpp"
#include "rectangle.hpp"

#define REQUEST 1
#define REPLY 2
#define TAG_Gi 10
#define TAG_Gj 11
#define TAG_Gval 12
#define TAG_Bi 20
#define TAG_Bval 21
#define DEBUG

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

int ReadFromFile(std::string file,
                 float &h,
                 float &T,
                 float &fac,
                 int &N_tray,
                 int &N_subdo,
                 int &N_node, 
                 Eigen::VectorXf & SW,
                 Eigen::VectorXf & NE);

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
    //Indexes of interfaces, nodes and subdomains
    std::vector<std::vector<int> > interface_indexes, node_indexes, node_interface, subd_index;
    //Number of interfaces and nodes in each direction.
    std::vector<int> i_N, n_N;
    //Initialization of these variables
    std::string file = "configuration.txt";
    float h, fac, T_start;
    int N_tray,N_subdo, N_node; 
    Eigen::VectorXf SW,NE;
    ReadFromFile(file, h, T_start, fac, N_tray,N_subdo,N_node, SW, NE);
    i_N.push_back(N_subdo); i_N.push_back(N_subdo);
    n_N.push_back(N_node); n_N.push_back(N_node);
    if(myid == server){
        printf("*************PDD SPARSE*****************\n");
        printf("***************Ver 0.2******************\n");
        printf("***********Jorge Moron Vidal************\n");
        printf("**********jmoron@math.uc3m.es***********\n");
        printf("h = %e \t fac =%e \t T_start = %0.2f\n",h,fac, T_start);
        printf("Number of subdomains: %d x %d \n", i_N[0], i_N[1]);
        printf("Number of nodes per interface: %d\n",n_N[0]);
        printf("Number of trayectories per node %d\n",N_tray);
        printf("Position of South-West corner of the domain: (%.2f,%.2f)\n",SW(0), SW(1));
        printf("Position of North-East corner of the domain: (%.2f,%.2f)\n",NE(0), NE(1));
    }
    //direction of the interface: 0 if horizontal 1 if vertical
    int dir;

    //Interior is 0 if the interface is not interior 1 if it is
    bool interior;

    //true if the nodes of the interfaces are distirbuted by chebysehv distribution
    bool chebyshev = false;

    //Parameters of the boundaries
    float global_p[4];
    global_p[0] = SW(0);
    global_p[1] = SW(1);
    global_p[2] = NE(0);
    global_p[3] = NE(1);

    //It gives the index where is stored the information of a given interface in a vector
    std::map<std::vector<int>, int> interface_map;
    //Positions of the nodes 
    std::vector<Eigen::VectorXf> node_pos;
    //Initializes the interface metadata variables in function of the specified parameters
    interface_metadata(node_pos, node_indexes, interface_indexes, interface_map, i_N, n_N, SW, NE);

    fullfill_node_interface(node_pos.size(),node_interface, node_indexes);

    #ifdef DEBUG
    FILE * pFile;
    if(myid == server){
        pFile = fopen ("Output/Debug/Node_position.txt","w");
        if (pFile!=NULL)
        {   fprintf(pFile, "index,x,y\n");
            for( unsigned int i = 0; i < node_pos.size(); i++) fprintf(pFile, "%u,%.3f,%.3f \n",i,node_pos[i](0),node_pos[i](1));
        } else {
            std::cout << "Something failed opening Output/Debug/Node_position.txt \n";
        }
        fclose(pFile);
        pFile = fopen("Output/Debug/boundary_global.txt", "w");
        fprintf(pFile,"x,y\n");
        fprintf(pFile,"%.3f,%.3f\n",global_p[0], global_p[1]);
        fprintf(pFile,"%.3f,%.3f\n",global_p[2], global_p[1]);
        fprintf(pFile,"%.3f,%.3f\n",global_p[2], global_p[3]);
        fprintf(pFile,"%.3f,%.3f\n",global_p[0], global_p[3]);
        fprintf(pFile,"%.3f,%.3f\n",global_p[0], global_p[1]);
        fclose(pFile);
        //system("python3 Plot_Node_Position.py");
    }
    #endif
    //c2 is defined
    float c2 = pow(fac*(node_pos[0]-node_pos[1]).norm(),2.0);
    if(myid == server) std::cout <<"c2 is " << c2 << "\n";

    //Subdomain is resized
    subd_index.resize(interface_indexes.size());
    //work_control controls the flow of the main program

    int work_control[2] = {0, 0}, aux = 0;
    if(myid == server){
    //I am the server who distributes task
        printf("Node_id\tProcess\tTime\n");
        do {
            aux = work_control[0];
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
                     REQUEST, world, &status);

            work_control[0] = aux;
            MPI_Send(work_control, 2, MPI_INT, 
                    status.MPI_SOURCE, REPLY, world);

            work_control[0]++;

        }while(work_control[0] < (int)node_pos.size());

        uint32_t N_node = node_pos.size();
        Eigen::SparseMatrix<float> G, I, B;
        I.resize(N_node, N_node);
        I.setIdentity();
        G.resize(N_node, N_node);
        B.resize(N_node, 1);
        
        typedef Eigen::Triplet<float,int> T;
        std::vector<T> T_vec_G, T_vec_B;
        int G_i[N_node*N_node], G_j[N_node*N_node], B_i[N_node];
        float G_val[N_node*N_node], B_val[N_node];
        for(int i = 0; i < server; i++){
            //It ends worker process solving
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
                     REQUEST, world, &status);
            work_control[0] = 0;
            work_control[1] = 1;
            MPI_Send(work_control, 2, MPI_INT, 
                status.MPI_SOURCE, REPLY, world);

            //The server starts accumulating the values of G and B
            MPI_Recv(work_control, 2, MPI_INT, status.MPI_SOURCE, 
                     REQUEST, world, &status);
            //G_i is received
            MPI_Recv(G_i, work_control[0], MPI_INT, status.MPI_SOURCE,
                    TAG_Gi, world, &status);
            //G_j is received
            MPI_Recv(G_j, work_control[0], MPI_INT, status.MPI_SOURCE,
                    TAG_Gj, world, &status);
            //G_val is received
            MPI_Recv(G_val, work_control[0], MPI_FLOAT, status.MPI_SOURCE,
                    TAG_Gval, world, &status);

            for(int j = 0; j < work_control[0]; j++){
                T_vec_G.push_back(T(G_i[j], G_j[j], G_val[j]));
            }

            //B_i is received
            MPI_Recv(B_i, work_control[1], MPI_INT, status.MPI_SOURCE,
                    TAG_Bi, world, &status);
            //B_val is received
            MPI_Recv(B_val, work_control[1], MPI_FLOAT, status.MPI_SOURCE,
                    TAG_Bval, world, &status);

            for(int j = 0; j < work_control[1]; j++){
                T_vec_B.push_back(T(B_i[j], 0, B_val[j]));
            }
        }
        Eigen::VectorXf Bd,ud;
        G.setFromTriplets(T_vec_G.begin(),T_vec_G.end());
        G+=I;
        #ifdef DEBUG 
        std::ofstream ofile("Output/Debug/G_m.txt");
        /*for (int k=0; k<G.outerSize(); ++k)
            for (Eigen::SparseMatrix<float>::InnerIterator it(G,k); it; ++it)
           {
                ofile << it.value() << "\t"
                << it.row() << "\t"   // row index
                << it.col() << "\n";   // col index (here it is equal to k)
            }
        */
        ofile << G;
        ofile.close();
        #endif
        B.setFromTriplets(T_vec_B.begin(), T_vec_B.end());
        Bd = Eigen::VectorXf(B);
        B.resize(0,0);
        T_vec_G.resize(0);
        T_vec_B.resize(0);
        //Firs Step in the solution's computation
        Eigen::SparseLU<Eigen::SparseMatrix<float> > solver_LU;
        solver_LU.compute(G);
        ud = solver_LU.solve(Bd);
        //Second Step in the solution's computation
        Eigen::BiCGSTAB<Eigen::SparseMatrix<float> > solver_IT;
        solver_IT.compute(G);
        ud = solver_IT.solveWithGuess(Bd,ud);

        ofile.open("Output/Debug/B_m.txt");
        for(int k = 0; k < (int)Bd.size(); k++){
            ofile << Bd(k) << std::endl;
        }
        ofile.close();
        std::cout <<" u was computed with error "<< solver_IT.error() << "\n";
        ofile.open("Output/solution.txt");
        ofile << "Node_i, Node_x, Node_y,sol_analytic,sol_PDDS,err\n";
        for(int i = 0; i< (int)node_pos.size(); i++){
           ofile << i <<"," <<node_pos[i](0) << "," << node_pos[i](1) << ","
           << Equation_u(node_pos[i], T_start) << "," << ud(i) <<","<<
           fabs(Equation_u(node_pos[i], T_start)-ud(i))/Equation_u(node_pos[i], T_start)<<"\n";
        }
        ofile.close();
        for(int i = 0; i< (int)node_pos.size(); i++){
           std::cout<< "Node " << i << "("<< node_pos[i](0) << "," << node_pos[i](1) << ") PDDS = " << ud(i) <<
           " Analytic = "<< Equation_u(node_pos[i], T_start) << 
           " r.err = " << fabs(Equation_u(node_pos[i], T_start)-
            ud(i))/Equation_u(node_pos[i], T_start)<<"\n";
        }



    } else {
        //Node object
        Node node_s;
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
        clock_t tic, toc;
        //Boundary Value problem
        BVP bvp;
        std::map<std::string, pfscalar> scalar_init;
        std::map<std::string, pfvector> vector_init;
        std::map<std::string, pfmatrix> matrix_init;
        std::map<std::string, pfscalarN> scalarN_init;
        std::map<std::string, std::string> string_init;
        //G and B storage vector
        std::vector<float> G, B;
        std::vector<int>  G_j, G_i, B_i;
        //G and B storage vectors for each node
        float B_temp;
        std::vector<float> G_temp;
        std::vector<int>  G_j_temp;
        //Stencil
        Stencil stencil;
        //BVP initialization 
        scalar_init["f"] = Equation_f;
        scalar_init["c"] = Equation_c;
        scalar_init["u"] = Equation_u;
        scalar_init["g"] = Equation_g;
        scalar_init["p"] = Equation_p;
        vector_init["b"] = Equation_b;
        matrix_init["sigma"] = Equation_sigma;
        bvp.Boundary_init(Rectangle2D, Stopping);
        bvp.BVP_init(2,scalar_init,scalarN_init,vector_init, matrix_init,string_init, Equation_RBF);
        while(work_control[1] == 0){

            MPI_Send(work_control, 2, MPI_INT, server, REQUEST, world);
            MPI_Comm_rank(workers, &workerid);

            MPI_Recv(work_control, 2, MPI_INT, server, REPLY, world, 
            MPI_STATUS_IGNORE);

            if(work_control[1] == 0){
                std::cout << "Processor " << myid << " is solving " <<work_control[0] << std::endl; 
                tic = clock();
                G_temp.clear();
                B_temp = 0.0f;
                G_j_temp.clear();
                G_temp.clear();
                stencil.Reset();
                if(node_interface[work_control[0]].size() == 1){
                    set_direction_interior(interior, dir, node_interface[work_control[0]][0], interface_indexes, i_N);
                    stencil.Init(dir, interior, interface_indexes[node_interface[work_control[0]][0]], 
                                 node_pos, node_indexes, interface_map, i_N, global_p);
                    
                    node_s.init(node_pos[work_control[0]],
                                h,
                                N_tray,
                                work_control[0],
                                interface_indexes[node_interface[work_control[0]][0]],
                                interface_indexes[node_interface[work_control[0]][0]]);
                } else {
                    stencil.Init(interface_indexes[node_interface[work_control[0]][0]],
                                 interface_indexes[node_interface[work_control[0]][1]],
                                 interface_indexes[node_interface[work_control[0]][2]],
                                 interface_indexes[node_interface[work_control[0]][3]],
                                 node_pos, node_indexes, interface_map, i_N, global_p);
                    node_s.init(node_pos[work_control[0]],
                                h,
                                N_tray,
                                work_control[0],
                                interface_indexes[node_interface[work_control[0]][1]],
                                interface_indexes[node_interface[work_control[0]][1]]);
                }
                stencil.Compute_ipsi(bvp, c2);
                stencil.Print(work_control[0]);
                Print_interface(work_control[0], node_interface, node_indexes, node_pos);
                node_s.Solve_PDDSparse(T_start, bvp, rng, N_tray, c2, stencil, G_j_temp, G_temp, B_temp);
                B.push_back(B_temp);
                B_i.push_back(work_control[0]);
                for(unsigned int k = 0;k < G_temp.size(); k ++){
                    G_i.push_back(work_control[0]);
                    G_j.push_back(G_j_temp[k]);
                    G.push_back(G_temp[k]);
                }
                toc = clock();
                printf("%d\t%d\t%0.3f\n", work_control[0], myid, ((float)(toc - tic) / CLOCKS_PER_SEC));
            }
        }
        work_control[0] = (int) G.size();
        work_control[1] = (int) B.size();
        MPI_Send(work_control, 2, MPI_INT, server, REQUEST, world);
        //G_i is sent
        MPI_Send(&G_i[0], work_control[0], MPI_INT, server, TAG_Gi, world);
        //G_j is sent
        MPI_Send(&G_j[0], work_control[0], MPI_INT, server, TAG_Gj, world);
        //G_val is sent
        MPI_Send(&G[0], work_control[0], MPI_FLOAT, server,TAG_Gval, world);
        //B_i is sent
        MPI_Send(&B_i[0], work_control[1], MPI_INT, server, TAG_Bi, world);
        //B_val is sent
        MPI_Send(&B[0], work_control[1], MPI_FLOAT, server, TAG_Bval, world);
        std::cout << "Processor " << workerid << " ended its work \n";
        MPI_Comm_free(&workers);
    }
    MPI_Finalize();
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

                vaux(j) = (NE(j) - SW(j))/(1.0f* ((n_node[j])*n_inter[j]));

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

                vaux(j) = (NE(j) - SW(j))/(1.0f* (n_inter[j]));

            } else {

                vaux(j) = 0.0f;

            }
        }

        lincrement.push_back(vaux);
    }
    //Initial Point
   
    for (int i = 0; i < SW.size(); i++){
            vaux = SW;
            P_cent.push_back(vaux);
    }

    //This section is not valid for + 2D problems

    P_cent[0] += lincrement[1];
    P_cent[1] += lincrement[0];

    int cen = 0;
    for (int i = 0; i < n_inter[1] - 1; i++){

        for (int j = 0; j < n_inter[0]; j++){

            start.push_back(P_cent[0] + lincrement[0]*j);
            end.push_back(P_cent[0] + lincrement[0]*(j+1));

            inter_indexes[cen].push_back(j);
            inter_indexes[cen].push_back(2*i + 1);
        
            cen++;

            start.push_back(P_cent[1] + lincrement[1]*j);
            end.push_back(P_cent[1] + lincrement[1]*(j+1));

            inter_indexes[cen].push_back(i);
            inter_indexes[cen].push_back(2*j);

            cen++;
        }


        P_cent[0] += lincrement[1];
        P_cent[1] += lincrement[0];

    }

    for(int i = 0; i < (int)inter_indexes.size(); i++){
        interface_map[inter_indexes[i]] = i;
    }
    unsigned int N_nodes = fullfill_node_index(node_indexes, n_inter, n_node, interface_map);
    fullfill_node_pos(node_positions, inter_indexes, node_indexes, interface_map,start, end, N_nodes);

}
unsigned int fullfill_node_index(std::vector< std::vector<int> > & node_indexes, 
                         std::vector<int> n_inter, 
                         std::vector<int> n_node, 
                         std::map<std::vector<int>, int> interface_map){
    //Vaux stores the current interface index.
    //cross stores the index of the nodes that constitutes 
    std::vector<int> vaux, cross;
    unsigned int index, cent = 0, cross_cent = 0;
    cross.resize((n_inter[0]-1)*(n_inter[1]-1));
    //As It is done in Paco's original code we start by fullfilling the vertical interfaces
    cross.clear();
    for(int i = 0; i < n_inter[0] - 1; i++){
        for(int j = 0; j <= n_inter[1]*2 - 2; j += 2){

            vaux.clear();
            vaux.push_back(i);
            vaux.push_back(j);
            index = interface_map[vaux];
            node_indexes[index].push_back(cent);
            if(j != 0) cross[(j/2-1)*(n_inter[1]-1) + i] = cent;
            for(int k = 0; k < n_node[1]; k ++){
                cent ++;
                node_indexes[index].push_back(cent);
            }
            if(j == n_inter[1]*2 -2) {
                cent ++;
            }
        }
    }
    //Horizontal interfaces
    for(int j = 1; j <= 2*n_inter[1]-3; j +=2 ){
        for(int i = 0; i < n_inter[0]; i ++){

            vaux.clear();
            vaux.push_back(i);
            vaux.push_back(j);
            index = interface_map[vaux];
            if(i != 0){
                node_indexes[index].push_back(cross[cross_cent]);
                cross_cent ++;
            } else {
                node_indexes[index].push_back(cent);
                cent ++;
            }
            for(int k = 0; k < n_node[0] - 1; k ++){
                node_indexes[index].push_back(cent);
                cent ++;
            }
            if(i != n_inter[0] -1 ) {
                node_indexes[index].push_back(cross[cross_cent]);
            } else {
                node_indexes[index].push_back(cent);
                cent ++;
            }
        }
    }
    return cent;
}
void fullfill_node_pos(std::vector<Eigen::VectorXf> & node_pos,
                       std::vector<std::vector<int> > i_index,
                       std::vector<std::vector<int> > n_index,
                       std::map<std::vector<int>, int> interface_map,
                       std::vector<Eigen::VectorXf> start,
                       std::vector<Eigen::VectorXf> end,
                       unsigned int N_nodes){
    
    //Auxiliary vectors 
    Eigen::VectorXf vaux, sincrement;
    int index;
    node_pos.clear();
    //We resize the node positions vector
    node_pos.resize(N_nodes);
    //Loop over the points in the system
    for(std::vector<std::vector<int> >::iterator it = i_index.begin();
    it != i_index.end();
    it++){
        index = interface_map[*it];
        sincrement = (end[index]-start[index])/(n_index[index].size()-1);
        vaux = start[index];
        for(int i = 0; i < (int) n_index[index].size(); i++){
            if(node_pos[n_index[index][i]].size() == 0){
                node_pos[n_index[index][i]] = vaux;
            }
            vaux += sincrement;
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

void fullfill_node_interface(unsigned int N_node,
                             std::vector <std::vector<int> >& node_interface, 
                             std::vector <std::vector<int> > node_indexes){
    node_interface.clear();
    node_interface.resize(N_node);
    
    for(int i = 0; i < (int)node_indexes.size(); i++){
        for(int j = 0; j < (int)node_indexes[i].size(); j++){
            node_interface[node_indexes[i][j]].push_back(i);
        }
    }
    FILE *pf;
    pf = fopen("Output/Debug/Node_interface.txt", "w");
    for(int i = 0; i < (int)node_interface.size(); i++){
        for(int j = 0; j < (int)node_interface[i].size(); j++){
            fprintf(pf,"%d %d ",i,node_interface[i][j]);
        }
        fprintf(pf,"\n");
    }
    fclose(pf);

}

void Print_interface(int node,
                    std::vector< std::vector<int> > node_interface,
                    std::vector< std::vector<int> > node_index,
                    std::vector< Eigen::VectorXf > node_pos){
    char filename[100];
    FILE *pFile;
    sprintf(filename,"Output/Debug/Interface_%d.txt", node);
    pFile = fopen(filename, "w");
    int index;
    fprintf(pFile,"index,x,y\n");
    for(int i = 0; i < (int) node_interface[node].size(); i++){
        for(int j = 0; j < (int) node_index[node_interface[node][i]].size(); j++){
            index = node_index[node_interface[node][i]][j];
            fprintf(pFile,"%d,%f,%f\n",index, node_pos[index](0),node_pos[index](1));
        }
    }
    fclose(pFile);
}

int ReadFromFile(std::string file,
                 float &h,
                 float &T,
                 float &fac,
                 int &N_tray,
                 int &N_subdo,
                 int &N_node, 
                 Eigen::VectorXf & SW,
                 Eigen::VectorXf & NE)
  {
    std::ifstream myfile(file.c_str());
    std::string line, cent;
    //If T is not specified then the problem is suposed to be elliptic
    T = INFINITY;
  if (myfile)  // same as: if (myfile.good())
    {   std::cout << "READING\n";
        while (getline( myfile, line )){  // same as: while (getline( myfile, line ).good())

            if (line.c_str()[0]!='%') //First of all we check if the line is a comment
            { 
                cent.clear(); //Centinel needs to be empty

                for(std::string::iterator it = line.begin(); it != line.end(); it++){
                //Loop over all the characters of the line
                    
                    if(*it != ' ' && *it != '\t'){
                    //If it is not a space nor a tabulation then we concatenate the term to cent
                        cent += *it;
                    }
                    else{
                    //If we get a space or a tabulation then we check the string
                        if(cent == "h="){
                        //TypeOfPDE
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            h = atof(cent.c_str());
                            break;
                        }
                        if(cent == "T="){
                        //TypeOfPDE
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            T = atof(cent.c_str());
                            break;
                        }
                        if(cent == "fac="){
                        //TypeOfPDE
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            fac = atof(cent.c_str());
                            break;
                        }

                        if(cent == "SW="){
                        //South West point of the domain 
                            SW.resize(2);
                            while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            for(int i = 0; i < 2; i++){
                                while((*it != ' ') && (*it != '\t') && (it != line.end())){

                                cent += *it;
                                it++;

                                }
                                SW(i) = atof(cent.c_str());
                                cent.clear();
                                while( *it == ' ' or *it == '\t'){
                                        it++;
                                }
                    
                            }

                            break;
                        }

                        if(cent == "NE="){
                        //South West point of the domain 
                            NE.resize(2);
                            while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            for(int i = 0; i < 2; i++){
                            while((*it != ' ') && (*it != '\t') && (it != line.end())){

                                cent += *it;
                                it++;

                            }
                            NE(i) = atof(cent.c_str());
                            cent.clear();
                            while( *it == ' ' or *it == '\t'){
                                        it++;
                            }
                    
                            }

                            break;
                        }

                        if(cent == "N_tray="){
                        //Number of trayectories per node

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            N_tray = atoi(cent.c_str());
                            break;
                        }

                        if(cent == "N_subdo="){
                        //Number of domains in each direction

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            N_subdo = atoi(cent.c_str());
                            break;
                        }

                        if(cent == "N_node="){
                        //Number nodes per interface

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            N_node = atoi(cent.c_str());
                            break;
                        }
                        
                    }
                }
            }
        }
        myfile.close();
        return 1;
    }

    else {
    
        std::cout << "WARNING: Cofiguration file was not found\n";
        std::cout << "Make sure it exist a configuration.txt file in the program's directory.\n";
        return 0;
    }    
 }


