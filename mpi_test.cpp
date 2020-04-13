#include <mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <vector>

#define REQUEST 1
#define REPLY 2
/*
void interface_metadata(std::vector<Eigen::VectorXf> & start, 
                        std::vector<Eigen::VectorXf> & end,
                        std::vector<std::vector<int>> & inter_indexes,
                        std::vector<std::vector<int>>  & node_indexes,
                        std::vector<int> n_inter,
                        std::vector<int> n_node,
                        Eigen::VectorXf  SW,
                        Eigen::VectorXf NE);
*/
int main(int argc, char *argv[]) {
    
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
    MPI_Group_excl(world_group, 1, ranks, &worker_group);
    MPI_Comm_create(world, worker_group, &workers);
    MPI_Group_free(&worker_group);
    MPI_Group_free(&world_group);
    std::vector<int> i_index;
    i_index.resize(100);

    int work_control[2] = {0, 0}, aux = 0;

    if(myid == server){
    /*I am the server who distributes task and stores the solution*/
        do {
            aux = work_control[0];
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
                     REQUEST, world, &status);
            work_control[0] = aux;
            MPI_Send(work_control, 2, MPI_INT, 
                    status.MPI_SOURCE, REPLY, world);

            work_control[0]++;

        }while(work_control[0] < (int)i_index.size());

        for(int i = 0; i < server; i++){
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
                     REQUEST, world, &status);
            work_control[0] = 0;
            work_control[1] = 1;
            MPI_Send(work_control, 2, MPI_INT, 
                status.MPI_SOURCE, REPLY, world);
        }

    } else {

        while(work_control[1] == 0){

            MPI_Send(work_control, 2, MPI_INT, server, REQUEST, world);
            MPI_Comm_rank(workers, &workerid);

            MPI_Recv(work_control, 2, MPI_INT, server, REPLY, world, 
            MPI_STATUS_IGNORE);

            if(work_control[1] == 0){

                std::cout << "Process " << workerid << " is solving itfc "
                << work_control[0] << '\n';
                sleep(1);
            }
        }

        std::cout << "Processor " << workerid << " ended its work \n";
        MPI_Comm_free(&workers);
    }
    MPI_Finalize();
}
    

/*
void interface_metadata(std::vector<Eigen::VectorXf> & start, 
                        std::vector<Eigen::VectorXf> & end,
                        std::vector<std::vector<int>> & inter_indexes,              
                        std::vector<std::vector<int>> & node_indexes,
                        std::vector<int> n_inter,
                        std::vector<int> n_node,
                        Eigen::VectorXf  SW,
                        Eigen::VectorXf NE){

     Initial values, long and  and shot increments for each direction
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

            for (int j = 0; j < SW.size(); j++){

                vaux(j) = sincrement[i](j); 

            }

            if(i > 0 ){

                vaux(i-1) += 0.5* sincrement[i-1](i-1);

            }

            P_cent.push_back(vaux);
    }

    //This section is not valid for + 2D problems

    P_cent[0] += lincrement[1];

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
        }


        P_cent[0] += lincrement[1];

    }

    P_cent[1] += lincrement[0];

    for (int i = 0; i < n_inter[0] - 1; i++){

        for (int j = 0; j < n_inter[1]; j++){

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


        P_cent[1] += lincrement[0];

    }


    for(int i = 0; i < 2; i++){

        std::cout <<"sincrement " << i;

        for(int j = 0; j < 2; j++){

            std::cout << " " << sincrement[i](j);
        }
        std::cout <<'\n';
    }
    for(int i = 0; i < 2; i++){

        std::cout <<"lincrement " << i;

        for(int j = 0; j < 2; j++){

            std::cout << " " << lincrement[i](j);
        }

        std::cout <<'\n';
    }
    
    cen = 0;
    std::cout << "Horizontal" << '\n';
    for(int i = 0; i < n_inter[0]*(n_inter[1]-1); i++){
        std::cout << "Interface" << '\n';
        for(int j = 0; j < 2; j++){

            std::cout << " " << inter_indexes[i][j];
        }
        std::cout << '\n';
        std::cout <<"start " ;

        for(int j = 0; j < 2; j++){

            std::cout << " " << start[i](j);
        }

        std::cout <<" \n end " ;

        for(int j = 0; j < 2; j++){

            std::cout << " " << end[i](j);
        }
        std::cout << '\n' << "Node index:";
        for (int l = 0; l < n_node[0]; l++){

                std::cout << " " << node_indexes[i][l];
                cen ++;

        }

        std::cout <<'\n';
    }

    std::cout << "Vertical" << '\n';
    for(int i = n_inter[0]*(n_inter[1]-1); i < 2*n_inter[0]*(n_inter[1]-1); i++){
        for(int j = 0; j < 2; j++){

            std::cout << " " << inter_indexes[i][j];
        }
        std::cout << '\n';
        std::cout <<"start " ;

        for(int j = 0; j < 2; j++){

            std::cout << " " << start[i](j);
        }
        std::cout <<"\n end " ;

        for(int j = 0; j < 2; j++){

            std::cout << " " << end[i](j);
        }
        std::cout << '\n' << "Node index:";
        for (int l = 0; l < n_node[0]; l++){

                std::cout << " " << node_indexes[i][l];
                cen ++;

        }

        std::cout <<'\n';
    }
    
}
*/