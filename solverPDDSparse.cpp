#include"solverPDDSparse.hpp"

SolverPDDS::SolverPDDS(int argc, char *argv[]){
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
}

void SolverPDDS::interface_metadata(std::vector<Eigen::VectorXf> & node_positions,
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
unsigned int SolverPDDS::fullfill_node_index(std::vector< std::vector<int> > & node_indexes, 
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
void SolverPDDS::fullfill_node_pos(std::vector<Eigen::VectorXf> & node_pos,
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

void SolverPDDS::set_direction_interior(bool & interior,
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

void SolverPDDS::fullfill_node_interface(unsigned int N_node,
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

void SolverPDDS::Print_interface(int node,
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


