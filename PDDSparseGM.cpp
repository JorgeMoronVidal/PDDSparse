#include "PDDSparseGM.hpp"
PDDSJob::PDDSJob(void){
    index[0] = -1;
    N[0] = 0;
    N[1] = 0;
}
void PDDSJob::Send_To_Worker(MPI_Status & status, MPI_Comm & world){
    MPI_Send(N, 2, MPI_UNSIGNED, status.MPI_SOURCE, TODO_JOB_UINT, world);
    printf("N sent from server to %d\n", status.MPI_SOURCE);
    MPI_Send(index, 1, MPI_INT, status.MPI_SOURCE, TODO_JOB_INT, world);
    printf("Index sent from server to %d\n", status.MPI_SOURCE);
    MPI_Send(h, 1, MPI_DOUBLE, status.MPI_SOURCE, TODO_JOB_DOUBLE, world);
    printf("h sent from server to %d\n", status.MPI_SOURCE);
}
void PDDSJob::Recieve_From_Server( int server, MPI_Comm & world){
    int myid;
    MPI_Comm_rank(world, &myid);
    MPI_Recv(N, 2, MPI_UNSIGNED, server, TODO_JOB_UINT, world, MPI_STATUS_IGNORE);
    printf("N received from server in %d\n", myid);
    MPI_Recv(index, 1, MPI_INT, server, TODO_JOB_INT, world, MPI_STATUS_IGNORE);
    printf("Index received from server in  %d\n", myid);
    MPI_Recv(h, 1, MPI_DOUBLE, server, TODO_JOB_DOUBLE, world, MPI_STATUS_IGNORE);
    printf("h received from server in  %d\n", myid);
}
/*void PDDSJob::Send_To_Server(int server, MPI_Comm & world){
    int myid;
    MPI_Comm_rank(world, &myid);
    MPI_Send(N, 2, MPI_UNSIGNED, server, DONE_JOB_UINT, world);
    printf("N sent to server from %d\n", myid);
    MPI_Send(index, 1, MPI_INT, server, DONE_JOB_INT, world);
    printf("Index sent to server from %d\n", myid);
}
void PDDSJob::Recieve_From_Worker(MPI_Status status, MPI_Comm & world){
    MPI_Recv(N, 2, MPI_UNSIGNED, status.MPI_SOURCE, DONE_JOB_UINT, world, MPI_STATUS_IGNORE);
    printf("N received in server from %d\n", status.MPI_SOURCE);
    MPI_Recv(index, 1, MPI_INT, status.MPI_SOURCE, DONE_JOB_INT, world, MPI_STATUS_IGNORE);
    printf("Index received in server from %d\n", status.MPI_SOURCE);
}*/
void PDDSparseGM::MPI_Configuration(int argc, char *argv[]){
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
    sprintf(debug_fname,"Output/Debug/Debug_output_p%d.txt",myid);
    FILE *dfile;
    dfile = fopen(debug_fname,"w");
    fclose(dfile);
}
PDDSparseGM::PDDSparseGM(int argc, char *argv[]){
  iN.resize(2);
  nN.resize(2);
  h0 = 0.0;
  T_start =INFINITY;
  eps = 0.5;
  SW.resize(2);
  NE.resize(2);
  N = 0;
  fac = 1;
  MPI_Configuration(argc, argv);
}
PDDSparseGM::PDDSparseGM(int argc, char *argv[],std::string file){ 
  MPI_Configuration(argc, argv);
  ReadFromFile(file);
  parameters.resize(4);
  parameters[0] = SW[0];
  parameters[1] = SW[1];
  parameters[2] = NE[0];
  parameters[3] = NE[1];
  if(myid == server){
    Fullfill_interfaces();
  }
}
void PDDSparseGM::ReadFromFile(std::string file){
  std::ifstream myfile(file.c_str());
    std::string line, cent;
    //If T is not specified then the problem is suposed to be elliptic
    T_start = INFINITY;
    if (myfile)  // same as: if (myfile.good())
    {   while (getline( myfile, line )){  // same as: while (getline( myfile, line ).good())

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
                        if(cent == "h0="){
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            h0 = atof(cent.c_str());
                            break;
                        }
                        if(cent == "eps="){
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            eps = atof(cent.c_str());
                            break;
                        }
                        if(cent == "fac="){
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
                        if(cent == "T="){
                         
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            T_start = atof(cent.c_str());
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

                        if(cent == "N="){
                        //Number of trayectories per node

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                                N = atoi(cent.c_str()); 
                            break;
                        }
                        if(cent == "N_job="){
                        //Number of trayectories per node

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            while( it != line.end()){

                                cent += *it;
                                it++;

                            }
                            
                            N_job = atoi(cent.c_str());
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
                            iN.resize(2);
                            iN[0] = atoi(cent.c_str());
                            iN[1] = atoi(cent.c_str());
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
                            
                            nN.resize(2);
                            nN[0] = atoi(cent.c_str());
                            nN[1] = atoi(cent.c_str());
                            break;
                        }
                        
                    }
                }
            }
        }
        myfile.close();

    }

    else {
    
        std::cout << "WARNING: Cofiguration file was not found\n";
        std::cout << "Make sure it exist a configuration.txt file in the program's directory.\n";
    }    
}
void PDDSparseGM::Fullfill_interfaces(void){
   /*Initial values, long and  and shot increments for each direction*/
   std::vector<Eigen::VectorXd> P_cent, lincrement, sincrement, start, end;
   Eigen::VectorXd vaux;
   std::vector<std::vector<int>> node_index;
   interfaces.resize(iN[0]*(iN[1]-1)+ iN[1]*( iN[0]-1));
   subdomains.resize(iN[0]*(iN[1]-1)+ iN[1]*( iN[0]-1));
   node_index.resize(iN[0]*(iN[1]-1)+ iN[1]*( iN[0]-1));
   vaux.resize(SW.size());

   //Increments are defined
   //ShortIncrement
   for (int i = 0; i < SW.size() ; i++){
        for(int j = 0; j < vaux.size(); j++){
            if (i == j){

                vaux(j) = (NE(j) - SW(j))/(1.0* ((nN[j])*iN[j]));
            } else {

                vaux(j) = 0.0;

            }
        }
        sincrement.push_back(vaux);
    }
    //LongIncrement
    for (int i = 0; i < SW.size() ; i++){

        for(int j = 0; j < vaux.size(); j++){

            if (i == j){

                vaux(j) = (NE(j) - SW(j))/(1.0* (iN[j]));

            } else {

                vaux(j) = 0.0;

            }
        }

        lincrement.push_back(vaux);
    }
    //Initial Point
   
    for (int i = 0; i < SW.size(); i++){
            vaux = SW;
            P_cent.push_back(vaux);
    }
    //c2 is computed as the square of fac*internodal space.
    c2 = pow(fac*sincrement[0].norm(),2.0);
    //This section is not valid for + 2D problems

    P_cent[0] += lincrement[1];
    P_cent[1] += lincrement[0];

    int cen = 0;

    for (int i = 0; i < iN[1] - 1; i++){

        for (int j = 0; j < iN[0]; j++){

            start.push_back(P_cent[0] + lincrement[0]*j);
            end.push_back(P_cent[0] + lincrement[0]*(j+1));

            interfaces[cen].label.push_back(j);
            interfaces[cen].label.push_back(2*i + 1);
            subdomains[cen].label.push_back(j);
            subdomains[cen].label.push_back(2*i + 1);
            cen++;

            start.push_back(P_cent[1] + lincrement[1]*j);
            end.push_back(P_cent[1] + lincrement[1]*(j+1));

            interfaces[cen].label.push_back(i);
            interfaces[cen].label.push_back(2*j);
            subdomains[cen].label.push_back(i);
            subdomains[cen].label.push_back(2*j);

            cen++;
        }


        P_cent[0] += lincrement[1];
        P_cent[1] += lincrement[0];

    }

    for(int i = 0; i < (int)interfaces.size(); i++){
        interface_map[interfaces[i].label] = i;
    }
    #ifdef INTERSECTIONS_YES
    //stvaux stores the current interface index.
    //cross stores the index of the nodes that constitutes 
    std::vector<int> stvaux, cross;
    unsigned int index, cent = 0, cross_cent = 0;
    cross.resize((iN[0]-1)*(iN[1]-1));
    //As It is done in Paco's original code we start by fullfilling the vertical interfaces
    cross.clear();
    for(int i = 0; i < iN[0] - 1; i++){
        for(int j = 0; j <= iN[1]*2 - 2; j += 2){
            stvaux.clear();
            stvaux.push_back(i);
            stvaux.push_back(j);
            index = interface_map[stvaux];
            node_index[index].push_back(cent);
            if(j != 0) cross[(j/2-1)*(iN[1]-1) + i] = cent;
            for(int k = 0; k < nN[1]; k ++){
                cent ++;
                node_index[index].push_back(cent);
            }
            if(j == iN[1]*2 -2) {
                cent ++;
            }
        }
    }
    //Horizontal interfaces
    for(int j = 1; j <= 2*iN[1]-3; j +=2 ){
        for(int i = 0; i < iN[0]; i ++){
            stvaux.clear();
            stvaux.push_back(i);
            stvaux.push_back(j);
            index = interface_map[stvaux];
            if(i != 0){
                node_index[index].push_back(cross[cross_cent]);
                cross_cent ++;
            } else {
                node_index[index].push_back(cent);
                cent ++;
            }
            for(int k = 0; k < nN[0] - 1; k ++){
                node_index[index].push_back(cent);
                cent ++;
            }
            if(i != iN[0] -1 ) {
                node_index[index].push_back(cross[cross_cent]);
            } else {
                node_index[index].push_back(cent);
                cent ++;
            }
        }
    }
    #else 
    //stvaux stores the current interface index.
    std::vector<int> stvaux;
    unsigned int index, cent = 0;
    for(int i = 0; i < iN[0] - 1; i++){
        for(int j = 0; j <= iN[1]*2 - 2; j += 2){
            stvaux.clear();
            stvaux.push_back(i);
            stvaux.push_back(j);
            index = interface_map[stvaux];
            node_index[index].push_back(cent);
            cent ++;
            //printf("[%d % d] nN = %d \n", i,j, nN[1]);
            for(int k = 0; k < nN[1] -1; k ++){
                node_index[index].push_back(cent);
                cent ++;
            }
        }
    }
    //Horizontal interfaces
    for(int j = 1; j <= 2*iN[1]-3; j +=2 ){
        for(int i = 0; i < iN[0]; i ++){
            stvaux.clear();
            stvaux.push_back(i);
            stvaux.push_back(j);
            index = interface_map[stvaux];
            node_index[index].push_back(cent);
            cent ++;
            //printf("[%d % d] nN = %d \n", i,j, nN[0]);
            for(int k = 0; k < nN[0] - 1; k ++){
                node_index[index].push_back(cent);
                cent ++;
            }
        }
    }
    #endif
    nNodes = (int)cent;
    FILE *dfile;
    dfile = fopen(debug_fname,"a");
    fprintf(dfile,"nNodes is %d\n",nNodes);
    fclose(dfile);
    //Loop over the points in the system
    double disp;
    disp = 0.5*(end[0]-start[0]).norm()/(node_index[0].size());
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        index = interface_map[(*it).label];
        (*it).index.resize(node_index[index].size());
        (*it).solution.resize(node_index[index].size());
        (*it).position.resize(node_index[index].size());
        sincrement[0] = (end[index]-start[index])/(node_index[index].size()-1);
        #ifdef INTERSECTIONS_YES
        vaux = start[index];
        #else
        if((*it).label[1]%2 == 0){
            //Vertical
            if((*it).label[1] == 0){
                end[index][1] = end[index][1] - disp;
            } else {
                if((*it).label[1] == 2*(iN[1] - 1)){
                    start[index][1] = start[index][1] + disp;
                } else {
                    start[index][1] = start[index][1] + disp;
                    end[index][1] = end[index][1] - disp;
                }
            }
        } else {
            //Horizontal
            if((*it).label[0] == 0){
                end[index][0] = end[index][0] - disp;
            } else {
                if((*it).label[0] == iN[0] - 1){
                start[index][0] = start[index][0] + disp;
                } else {
                    start[index][0] = start[index][0] + disp;
                    end[index][0] = end[index][0] - disp;
                }
            }
        }
        sincrement[0] = (end[index]-start[index])/(node_index[index].size()-1);
        vaux = start[index];
        #endif
        for(int i = 0; i < (int) node_index[index].size(); i++){
            (*it).index[i] = node_index[index][i];
            (*it).position[i] = vaux;
            //printf("%d %f %f\n", node_index[index][i],(*it).position[i][0],(*it).position[i][1]);
            vaux += sincrement[0];
        }
    }
    //getchar();
}
void PDDSparseGM::Solve(BVP bvp){
    bool done=true;
    double start = MPI_Wtime();
    FILE *pFile;
    if(myid==server) Print_Problem();
    //Compute the PDDSparse Matrix
    if(myid==server){

        PDDSJob auxjob;
        //Job list is fullfilled
        for(int knot_index = 0; knot_index < nNodes; knot_index++){
            auxjob.index[0] = knot_index;
            auxjob.N[0] = N; auxjob.N[1] = N_job;
            for(int jobs_per_knot = 0; jobs_per_knot < N/N_job; jobs_per_knot ++){
                work_control[1] = 1;
                Send_Node(auxjob);
            }
        }
        //G and B are received from the workers
        for(int process = 0; process < server; process++){
             Receive_G_B_2();
        }
        Compute_Solution_2(bvp);
        pFile = fopen(debug_fname,"a");
        double end = MPI_Wtime();
        fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
        fprintf(pFile,"Process %d ended its work.\n",myid);
        fclose(pFile);
    } else {
        //Start and end time for the node
        double knot_start, knot_end;
        //G and B storage vectors for each node
        double B_temp;
        PDDSJob auxjob;
        std::vector<double> G_temp;
        std::vector<int>  G_j_temp;
        //Stencil
        Stencil stencil;
        GMSolver solver(bvp ,parameters, h0, (unsigned int) myid +1);
        c2 = pow(fac*(1.0/nN[0]),2.0);
        G_i.clear();G_j.clear(); G.clear(); G_CT.clear(); G_var.clear();
        B.clear();B_CT.clear();B_i.clear();B_var.clear();
        do{
            done = Receive_Node(auxjob);
            if(!done){
                knot_start = MPI_Wtime();
                G_temp.clear();
                B_temp = 0.0;
                G_j_temp.clear();
                G_temp.clear();
                stencil = Recieve_Stencil_Data();
                stencil.Compute_ipsi(bvp,c2,debug_fname);
                if(stencil.Is_Interior()){ 
                    solver.Solve(position,c2,stencil,G_j_temp,G_temp,B_temp,auxjob.N[0], auxjob.N[1]);
                }else{
                    solver.Solve_mix(position,c2,1.0,stencil,G_j_temp,G_temp,B_temp,auxjob.N[0], auxjob.N[1]);
                }
                B.push_back(B_temp);
                B_CT.push_back(solver.B_CT);
                B_i.push_back(auxjob.index[0]);
                B_var.push_back(solver.var_B);
                for(unsigned int k = 0;k < G_temp.size(); k ++){
                    G_i.push_back(auxjob.index[0]);
                    G_j.push_back(G_j_temp[k]);
                    G.push_back(G_temp[k]);
                    G_CT.push_back(solver.G_CT[k]);
                    G_var.push_back(solver.var_G[k]);
                }
                knot_end = MPI_Wtime();
                
            }
        }while(!done);
        pFile = fopen(debug_fname,"a");
        Send_G_B_2();
        double end = MPI_Wtime();
    }
    //MPI_Finalize();
}
void PDDSparseGM::Solve(BVP bvp, std::string file){

    ReadFromFile(file);
    if(myid == server){
          Fullfill_interfaces();
    }
    Solve(bvp);
}
void PDDSparseGM::Send_Node(PDDSJob job){
    int work_control_aux[2];
    MPI_Recv(work_control_aux, 2, MPI_INT, MPI_ANY_SOURCE, REQUEST_NODE, world, &status);
    MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_NODE, world);
    job.Send_To_Worker(status, world);
    double pos[2];
    pos[0] = Node_Position(job.index[0])[0];
    pos[1] = Node_Position(job.index[0])[1];
    MPI_Send(pos, 2, MPI_DOUBLE, status.MPI_SOURCE, NODE_POSITION, world);
    if(work_control[1] == 1){
        FILE *pFile;
        pFile = fopen ("Output/Debug/Node_position.txt","a");
        fprintf(pFile, "%d,%.4f,%.4f \n",job.index[0],pos[0],pos[1]);
        fclose(pFile);
        pFile = fopen(debug_fname, "a");
        fprintf(pFile,"Sending node %d.\n",job.index[0]);
        fclose(pFile);
        Send_Stencil_Data(job.index[0]);
    } else {
        printf("Ending process.\n");
    }
}
void PDDSparseGM::Send_Node_Loop(PDDSJob job){
    int work_control_aux[2];
    MPI_Recv(work_control_aux, 2, MPI_INT, MPI_ANY_SOURCE, REQUEST_NODE, world, &status);
    MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_NODE, world);
    job.Send_To_Worker(status, world);
    double pos[2];
    pos[0] = Node_Position(job.index[0])[0];
    pos[1] = Node_Position(job.index[0])[1];
    MPI_Send(pos, 2, MPI_DOUBLE, status.MPI_SOURCE, NODE_POSITION, world);
    if(work_control[1] == 1){
        FILE *pFile;
        pFile = fopen ("Output/Debug/Node_position.txt","a");
        fprintf(pFile, "%d,%.4f,%.4f \n",job.index[0],pos[0],pos[1]);
        fclose(pFile);
        pFile = fopen(debug_fname, "a");
        fprintf(pFile,"Sending node %d.\n",job.index[0]);
        fclose(pFile);
        Send_Stencil_Data_Loop(job.index[0]);
    } else {
        printf("Ending process.\n");
    }
}
std::vector<std::vector<int> > PDDSparseGM::Get_Interfaces(int index){
    std::vector<std::vector<int> > output;
    //printf("Node %d is in interfaces: ",index);
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        if((*it).In(index) == true) {
            output.push_back((*it).label);
            //printf("[%d %d]", (*it).label[0],(*it).label[1]);
        }
    }
    //printf("\n");
    return output;
}
std::map<direction, std::vector<std::vector <int> > > PDDSparseGM::Labels_Stencil(int index){
    std::map<direction, std::vector<std::vector <int> > > output;
    std::vector< std::vector<int> > node_interfaces;
    Interface inter;
    std::vector< int > interface_east, interface_west,interface_north, interface_south;
    int max_v = 0, max_h = 0;
    node_interfaces = Get_Interfaces(index);
    FILE *dfile;
    switch(node_interfaces.size()){
        case 1:
            inter = interfaces[interface_map[node_interfaces[0]]];
            for(int i = 0; i < (int) inter.N().size(); i ++){
                if(interface_map.count(inter.N()[i]) > 0) output[North].push_back(inter.N()[i]);
            }
            for(int i = 0; i < (int) inter.S().size(); i ++){
                if(interface_map.count(inter.S()[i]) > 0) output[South].push_back(inter.S()[i]);
            }
            for(int i = 0; i < (int) inter.E().size(); i ++){
                if(interface_map.count(inter.E()[i]) > 0) output[East].push_back(inter.E()[i]);
            }
            for(int i = 0; i < (int) inter.W().size(); i ++){
                if(interface_map.count(inter.W()[i]) > 0) output[West].push_back(inter.W()[i]);
            }
            break;
        case 4:
            max_h = std::max(node_interfaces[0][0], node_interfaces[1][0]);
            max_h = std::max(max_h, node_interfaces[2][0]);
            max_h = std::max(max_h, node_interfaces[3][0]);
            max_v = std::max(node_interfaces[0][1], node_interfaces[1][1]);
            max_v = std::max(max_v, node_interfaces[2][1]);
            max_v = std::max(max_v, node_interfaces[3][1]);
            if(node_interfaces[0][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[0][0] == max_h) interface_east = node_interfaces[0];
                else interface_west = node_interfaces[0];
            } else {
                //Vertical interface
                if(node_interfaces[0][1] == max_v) interface_north = node_interfaces[0];
                else interface_south = node_interfaces[0];
            }
            if(node_interfaces[1][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[1][0] == max_h) interface_east = node_interfaces[1];
                else interface_west = node_interfaces[1];
            } else {
                //Vertical interface
                if(node_interfaces[1][1] == max_v) interface_north = node_interfaces[1];
                else interface_south = node_interfaces[1];
            }
            if(node_interfaces[2][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[2][0] == max_h) interface_east = node_interfaces[2];
                else interface_west = node_interfaces[2];
            } else {
                //Vertical interface
                if(node_interfaces[2][1] == max_v) interface_north = node_interfaces[2];
                else interface_south = node_interfaces[2];
            }
            if(node_interfaces[3][1]%2 == 1){
                //Horizontal interface
                if(node_interfaces[3][0] == max_h) interface_east = node_interfaces[3];
                else interface_west = node_interfaces[3];
            } else {
                //Vertical interface
                if(node_interfaces[3][1] == max_v) interface_north = node_interfaces[3];
                else interface_south = node_interfaces[3];
            }
            inter = interfaces[interface_map[interface_west]];
            if(interface_map.count(inter.N()[0]) > 0) output[North].push_back(inter.N()[0]);
            inter = interfaces[interface_map[interface_east]];
            if(interface_map.count(inter.N()[0]) > 0) output[North].push_back(inter.N()[0]);
            inter = interfaces[interface_map[interface_west]];
            if(interface_map.count(inter.S()[0]) > 0) output[South].push_back(inter.S()[0]);
            inter = interfaces[interface_map[interface_east]];
            if(interface_map.count(inter.S()[0]) > 0) output[South].push_back(inter.S()[0]);
            inter = interfaces[interface_map[interface_south]];
            if(interface_map.count(inter.E()[0]) > 0) output[East].push_back(inter.E()[0]);
            inter = interfaces[interface_map[interface_north]];
            if(interface_map.count(inter.E()[0]) > 0) output[East].push_back(inter.E()[0]);
            inter = interfaces[interface_map[interface_south]];
            if(interface_map.count(inter.W()[0]) > 0) output[West].push_back(inter.W()[0]);
            inter = interfaces[interface_map[interface_north]];
            if(interface_map.count(inter.W()[0]) > 0) output[West].push_back(inter.W()[0]);
            break;
        default:
            dfile = fopen(debug_fname,"a");
            fprintf(dfile,"Something went wrong in Labels_Stencil.\n");
            fclose(dfile);
    }
    /*
    printf("South: ");
    for(int i = 0; i < (int) output[South].size(); i++) printf("[%d %d] ",output[South][i][0],output[South][i][1]);
    printf("\n");
    printf("North: ");
    for(int i = 0; i < (int) output[North].size(); i++) printf("[%d %d] ",output[North][i][0],output[North][i][1]);
    printf("\n");
    printf("East: ");
    for(int i = 0; i < (int) output[East].size(); i++) printf("[%d %d] ",output[East][i][0],output[East][i][1]);
    printf("\n");
    printf("West: ");
    for(int i = 0; i < (int) output[West].size(); i++) printf("[%d %d] ",output[West][i][0],output[West][i][1]);
    printf("\n"); */
    return output;
}
void PDDSparseGM::Send_Stencil_Data(int index){
    std::map<direction, std::vector<std::vector <int> > > labels = Labels_Stencil(index);
    //The parameters of the stencil are computed
    double s_params[4];
    for(int i = 0; i < 4; i++) s_params[i] = parameters[i];
    if(labels[South].size()>0){
        //s_params[0] = interfaces[interface_map[labels[South][0]]].position[0][0];
        s_params[1] = interfaces[interface_map[labels[South][0]]].position[0][1];
    }
    if(labels[West].size()>0){
        s_params[0] = interfaces[interface_map[labels[West][0]]].position[0][0];
        //s_params[1] = interfaces[interface_map[labels[West][0]]].position[0][1];
    }
    int aux_size;
    if(labels[North].size()>0){
        aux_size = (int) interfaces[interface_map[labels[North][labels[North].size()-1]]].position.size();
        //s_params[2] = interfaces[interface_map[labels[North][labels[North].size()-1]]].position[aux_size - 1][0];
        s_params[3] = interfaces[interface_map[labels[North][labels[North].size()-1]]].position[aux_size - 1][1];
    }
    if(labels[East].size()>0){
        aux_size = (int) interfaces[interface_map[labels[East][labels[East].size()-1]]].position.size();
        s_params[2] = interfaces[interface_map[labels[East][labels[East].size()-1]]].position[aux_size - 1][0];
        //s_params[3] = interfaces[interface_map[labels[East][labels[East].size()-1]]].position[aux_size - 1][1];
    }
    //The indexes and positions of the stencils are computed
    std::map<direction, std::vector<int> > s_index;
    std::map<direction, std::vector<double> > s_x, s_y;
    std::vector<direction> s_dir;
    Interface s_inter, aux_inter;
    std::vector<int> aux_label;
    FILE *dfile;
    s_dir.resize(4);
    s_dir[0] = North; s_dir[1] = South; s_dir[2] = East; s_dir[3] = West;
    for(std::vector<direction>::iterator it = s_dir.begin(); it != s_dir.end(); it++){
    if(labels[*it].size() > 0){
        s_inter = interfaces[interface_map[labels[*it][0]]];
        aux_label.resize(2);
        aux_label[0] = s_inter.label[0]; aux_label[1] = s_inter.label[1];
        switch(*it){
            case North:
                aux_label[0] --;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case South:
                aux_label[0] --;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
                break;
            case East:
                aux_label[1] += -2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case West:
                aux_label[1] += -2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = (int)((1-STEN_ELONG)*(aux_inter.index.size()-1));
                    i < (int) aux_inter.index.size(); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            default:
                dfile = fopen(debug_fname,"a");
                fprintf(dfile,"Something went wrong during stencil formation. 1\n");
                fclose(dfile);
        }
        for(int i = 0; i < (int) labels[*it].size(); i++){
            s_inter = interfaces[interface_map[labels[*it][i]]];
            for(int j = 0;
            j < (int) s_inter.index.size(); j ++){
                s_index[*it].push_back(s_inter.index[j]);
                s_x[*it].push_back(s_inter.position[j][0]);
                s_y[*it].push_back(s_inter.position[j][1]);
            }
        }
        aux_label.resize(2);
        aux_label[0] = s_inter.label[0]; aux_label[1] = s_inter.label[1];
        switch(*it){
            case North:
                aux_label[0] ++;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case South:
                aux_label[0] ++;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
                break;
            case East:
                aux_label[1] += 2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            case West:
                aux_label[1] += 2;
                if(interface_map.count(aux_label)>0){
                aux_inter = interfaces[interface_map[aux_label]];
                for(int i = 0;
                    i < (int) ((STEN_ELONG)*(aux_inter.index.size())); i ++){
                        s_index[*it].push_back(aux_inter.index[i]);
                        s_x[*it].push_back(aux_inter.position[i][0]);
                        s_y[*it].push_back(aux_inter.position[i][1]);
                    }
                }
            break;
            default:
                dfile = fopen(debug_fname,"w");
                fprintf(dfile,"Something went wrong during stencil formation. 2\n");
                fclose(dfile);
        }
    }
    }
    
    //Parameters of the stencils are sent to the worker
    MPI_Send(s_params, 4, MPI_DOUBLE, status.MPI_SOURCE, PARAMETERS, world);
    int sizes[4];
    sizes[0] = (int) s_index[North].size();sizes[1] = (int) s_index[South].size();
    sizes[2] = (int) s_index[East].size();sizes[3] = (int) s_index[West].size();
    //Sizes of the stencil's vectors are sent to the worker
    MPI_Send(sizes, 4, MPI_INT, status.MPI_SOURCE, SIZES, world);
    //Stencil indexes and positions are sent to the worker 
    int *i_north, *i_south, *i_west, *i_east;
    double *x_north, *x_south, *x_west, *x_east, 
           *y_north, *y_south, *y_west, *y_east;
    //North stencil
    i_north = new int[sizes[0]];
    x_north = new double[sizes[0]];
    y_north = new double[sizes[0]];
    for(int i = 0; i < sizes[0]; i++){
        i_north[i] = s_index[North][i];
        x_north[i] = s_x[North][i];
        y_north[i] = s_y[North][i];
    }
    MPI_Send(i_north, sizes[0], MPI_INT, status.MPI_SOURCE, IND_NORTH, world);
    MPI_Send(x_north, sizes[0], MPI_DOUBLE, status.MPI_SOURCE, X_NORTH, world);
    MPI_Send(y_north, sizes[0], MPI_DOUBLE, status.MPI_SOURCE, Y_NORTH, world);
    delete i_north; delete x_north; delete y_north;
    //South stencil
    i_south = new int[sizes[1]];
    x_south = new double[sizes[1]];
    y_south = new double[sizes[1]];
    for(int i = 0; i < sizes[1]; i++){
        i_south[i] = s_index[South][i];
        x_south[i] = s_x[South][i];
        y_south[i] = s_y[South][i];
    }
    MPI_Send(i_south, sizes[1], MPI_INT, status.MPI_SOURCE, IND_SOUTH, world);
    MPI_Send(x_south, sizes[1], MPI_DOUBLE, status.MPI_SOURCE, X_SOUTH, world);
    MPI_Send(y_south, sizes[1], MPI_DOUBLE, status.MPI_SOURCE, Y_SOUTH, world);
    delete i_south; delete x_south; delete y_south;
    //East stencil
    i_east = new int[sizes[2]];
    x_east = new double[sizes[2]];
    y_east = new double[sizes[2]];
    for(int i = 0; i < sizes[2]; i++){
        i_east[i] = s_index[East][i];
        x_east[i] = s_x[East][i];
        y_east[i] = s_y[East][i];
    }
    MPI_Send(i_east, sizes[2], MPI_INT, status.MPI_SOURCE, IND_EAST, world);
    MPI_Send(x_east, sizes[2], MPI_DOUBLE, status.MPI_SOURCE, X_EAST, world);
    MPI_Send(y_east, sizes[2], MPI_DOUBLE, status.MPI_SOURCE, Y_EAST, world);
    delete i_east; delete x_east; delete y_east;
    //West stencil
    i_west = new int[sizes[3]];
    x_west = new double[sizes[3]];
    y_west = new double[sizes[3]];
    for(int i = 0; i < sizes[3]; i++){
        i_west[i] = s_index[West][i];
        x_west[i] = s_x[West][i];
        y_west[i] = s_y[West][i];
    }
    MPI_Send(i_west, sizes[3], MPI_INT, status.MPI_SOURCE, IND_WEST, world);
    MPI_Send(x_west, sizes[3], MPI_DOUBLE, status.MPI_SOURCE, X_WEST, world);
    MPI_Send(y_west, sizes[3], MPI_DOUBLE, status.MPI_SOURCE, Y_WEST, world);
    delete i_west; delete x_west; delete y_west;
}
void PDDSparseGM::Send_Stencil_Data_Loop(int index){
    std::vector< std::vector<int> > node_interfaces;
    node_interfaces = Get_Interfaces(index);
    int aux_size[1];
    if(node_interfaces.size() == 1){
        int aux_int[2] = {node_interfaces[0][0], node_interfaces[0][1]};
        aux_size[0]=2;
        MPI_Send(aux_size,1,MPI_INT,status.MPI_SOURCE,INTER_NUMBER,world);
        MPI_Send(aux_int,2,MPI_INT,status.MPI_SOURCE,INTER_LABELS,world);
    }else{
        int aux_int[8] = {node_interfaces[0][0], node_interfaces[0][1],
                          node_interfaces[1][0], node_interfaces[1][1],
                          node_interfaces[2][0], node_interfaces[2][1],
                          node_interfaces[3][0], node_interfaces[3][1],
                        };
        aux_size[0]=8;
        MPI_Send(aux_size,1,MPI_INT,status.MPI_SOURCE,INTER_NUMBER,world);
        MPI_Send(aux_int,8,MPI_INT,status.MPI_SOURCE,INTER_LABELS,world);
    }
    Send_Stencil_Data(index);
}
bool PDDSparseGM::Receive_Node(PDDSJob & job){
    FILE *dfile;
    MPI_Comm_rank(workers, &workerid);
    MPI_Send(work_control, 2, MPI_INT, server, REQUEST_NODE, world);
    MPI_Recv(work_control, 2, MPI_INT, server, REPLY_NODE, world, MPI_STATUS_IGNORE);
    if(work_control[1] == 1){
        job.Recieve_From_Server(server, world);
        dfile = fopen(debug_fname,"a");
        fprintf(dfile,"Process %d is solving node %d .\n",workerid,job.index[0]);
        fclose(dfile);
        double pos[2];
        MPI_Recv(pos, 2, MPI_DOUBLE, server, NODE_POSITION, world, MPI_STATUS_IGNORE);
        position.resize(2);
        position[0] = pos[0]; position[1] = pos[1];
        return false;
    } else {
        FILE *dfile;
        dfile = fopen(debug_fname,"a");
        fprintf(dfile,"Process %d ended its receiving nodes.\n",workerid);
        fclose(dfile);
        return true;
    }
}
Stencil PDDSparseGM::Recieve_Stencil_Data(void){
    double s_params[4];
    MPI_Recv(s_params, 4, MPI_DOUBLE, server, PARAMETERS, world, MPI_STATUS_IGNORE);
    int sizes[4];
    MPI_Recv(sizes, 4, MPI_INT, server, SIZES, world, MPI_STATUS_IGNORE);
    int *i_north, *i_south, *i_west, *i_east;
    double *x_north, *x_south, *x_west, *x_east, 
           *y_north, *y_south, *y_west, *y_east;
    std::map<direction, std::vector<int> > s_index;
    std::map<direction, std::vector<double> > s_x, s_y;
    i_north = new int[sizes[0]]; x_north = new double[sizes[0]]; y_north = new double[sizes[0]];
    MPI_Recv(i_north, sizes[0], MPI_INT, server, IND_NORTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_north, sizes[0], MPI_DOUBLE, server, X_NORTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_north, sizes[0], MPI_DOUBLE, server, Y_NORTH, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[0]; i++){
        s_index[North].push_back(i_north[i]);
        s_x[North].push_back(x_north[i]);
        s_y[North].push_back(y_north[i]);
    }
    delete i_north; delete x_north; delete y_north;
    i_south = new int[sizes[1]]; x_south = new double[sizes[1]]; y_south = new double[sizes[1]]; 
    MPI_Recv(i_south, sizes[1], MPI_INT, server, IND_SOUTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_south, sizes[1], MPI_DOUBLE, server, X_SOUTH, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_south, sizes[1], MPI_DOUBLE, server, Y_SOUTH, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[1]; i++){
        s_index[South].push_back(i_south[i]);
        s_x[South].push_back(x_south[i]);
        s_y[South].push_back(y_south[i]);
    }
    delete i_south; delete x_south; delete y_south;
    i_east = new int[sizes[2]]; x_east = new double[sizes[2]]; y_east = new double[sizes[2]];
    MPI_Recv(i_east, sizes[2], MPI_INT, server, IND_EAST, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_east, sizes[2], MPI_DOUBLE, server, X_EAST, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_east, sizes[2], MPI_DOUBLE, server, Y_EAST, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[2]; i++){
        s_index[East].push_back(i_east[i]);
        s_x[East].push_back(x_east[i]);
        s_y[East].push_back(y_east[i]);
    }
    delete i_east; delete x_east; delete y_east;
    i_west = new int[sizes[3]]; x_west = new double[sizes[3]]; y_west = new double[sizes[3]];
    MPI_Recv(i_west, sizes[3], MPI_INT, server, IND_WEST, world, MPI_STATUS_IGNORE);
    MPI_Recv(x_west, sizes[3], MPI_DOUBLE, server, X_WEST, world, MPI_STATUS_IGNORE);
    MPI_Recv(y_west, sizes[3], MPI_DOUBLE, server, Y_WEST, world, MPI_STATUS_IGNORE);
    for(int i = 0; i < sizes[3]; i++){
        s_index[West].push_back(i_west[i]);
        s_x[West].push_back(x_west[i]);
        s_y[West].push_back(y_west[i]);
    }
    delete i_west; delete x_west; delete y_west;
    Stencil output;
    output.Init(s_index,s_x,s_y,s_params);
    return output;
}
Stencil PDDSparseGM::Recieve_Stencil_Data_Loop(void){
    int size[1];
    MPI_Recv(size,1,MPI_INT, server, INTER_NUMBER, world, MPI_STATUS_IGNORE);
    int *Iindexes = new int[size[0]];
    MPI_Recv(Iindexes,size[0],MPI_INT,server, INTER_LABELS,world,MPI_STATUS_IGNORE);
    std::vector<std::vector<int>> subd_index;
    FILE *ox_file, *oy_file, *osol_file;
    char aux_fname[256];
    sprintf(aux_fname,"Input/Patch_%d/x_0.txt",myid);
    if((ox_file = fopen(aux_fname,"w"))==NULL){
        sprintf(aux_fname,"mkdir -p Input/Patch_%d",myid);
        system(aux_fname);
        sprintf(aux_fname,"Input/Patch_%d/x_0.txt",myid);
        ox_file = fopen(aux_fname,"w");
    }
    sprintf(aux_fname,"Input/Patch_%d/x_1.txt",myid);
    oy_file = fopen(aux_fname,"w");
    sprintf(aux_fname,"Input/Patch_%d/value.txt",myid);
    osol_file = fopen(aux_fname,"w");
    if(size[0] == 2){
        subd_index.resize(2);
        subd_index[0].resize(2);
        subd_index[1].resize(2);
        if(Iindexes[1]%2== 0){
            //horizontal patch
            std::cout << "Horizontal patch\n";
            subd_index[0][0] = Iindexes[0]; subd_index[0][1] = Iindexes[1]/2;
            subd_index[1][0] = Iindexes[0]+1; subd_index[1][1] = Iindexes[1]/2;
            std::vector<double> x_W, x_E;
            std::vector<double> y_W, y_E;
            std::vector<double> z_W, z_E;
            unsigned int i,j;
            sprintf(aux_fname,"Output/Subdomains/X_%d%d.txt",subd_index[0][0],subd_index[0][1]);
            x_W = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Y_%d%d.txt",subd_index[0][0],subd_index[0][1]);
            y_W = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Sol_%d%d.txt",subd_index[0][0],subd_index[0][1]);
            z_W = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/X_%d%d.txt",subd_index[1][0],subd_index[1][1]);
            x_E = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Y_%d%d.txt",subd_index[1][0],subd_index[1][1]);
            y_E = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Sol_%d%d.txt",subd_index[1][0],subd_index[1][1]);
            z_E = Read_File(aux_fname);
            for(j = 0; j < y_W.size(); j++){
                fprintf(oy_file,"%e\n",y_W[j]);
                for(i = 0; i < x_W.size(); i++){
                    fprintf(osol_file," %e",z_W[x_W.size()*j + i]);
                }
                for(i = 1; i < x_E.size(); i++){
                    fprintf(osol_file," %e",z_E[x_W.size()*j + i]);
                }
                fprintf(osol_file,"\n");
            }
            fclose(osol_file);
            fclose(oy_file);
            for(i = 0; i < x_W.size(); i++){
                    fprintf(ox_file,"%e\n",x_W[i]);
            }
            for(i = 1; i < x_E.size(); i++){
                    fprintf(ox_file,"%e\n",x_E[i]);
            }
            fclose(ox_file);
        }else{
            //vertical patch
            std::cout << "Vertical patch\n";
            subd_index[0][0] = Iindexes[0]; subd_index[0][1] = (Iindexes[1]-1)/2;
            subd_index[1][0] = Iindexes[0]; subd_index[1][1] = 1+((Iindexes[1]-1)/2);
            std::vector<double> x_N, x_S;
            std::vector<double> y_N, y_S;
            std::vector<double> z_N, z_S;
            unsigned int i,j;
            sprintf(aux_fname,"Output/Subdomains/X_%d%d.txt",subd_index[0][0],subd_index[0][1]);
            x_S = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Y_%d%d.txt",subd_index[0][0],subd_index[0][1]);
            y_S = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Sol_%d%d.txt",subd_index[0][0],subd_index[0][1]);
            z_S = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/X_%d%d.txt",subd_index[1][0],subd_index[1][1]);
            x_N = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Y_%d%d.txt",subd_index[1][0],subd_index[1][1]);
            y_N = Read_File(aux_fname);
            sprintf(aux_fname,"Output/Subdomains/Sol_%d%d.txt",subd_index[1][0],subd_index[1][1]);
            z_N = Read_File(aux_fname);
            for(j = 0; j < y_S.size(); j++){
                fprintf(oy_file,"%e\n",y_S[j]);
                for(i = 0; i < x_S.size(); i++){
                    fprintf(osol_file," %e",z_S[x_S.size()*j + i]);
                }
                fprintf(osol_file,"\n");
            }
            for(j = 1; j < y_N.size(); j++){
                fprintf(oy_file,"%e\n",y_N[j]);
                for(i = 0; i < x_S.size(); i++){
                    fprintf(osol_file," %e",z_N[x_S.size()*j + i]);
                }
                fprintf(osol_file,"\n");
            }
            fclose(osol_file);
            fclose(oy_file);
            for(i = 0; i < x_S.size(); i++){
                    fprintf(ox_file," %e\n",x_S[i]);
            }
            fclose(ox_file);
        }
    }else{
        //Knot belongs to 4 interfaces (It is on an intersection)
        std::cout << "Squared patch\n";
        std::vector<int> aux_index;
        aux_index.resize(2);
        subd_index.resize(2);
        std::vector<double> x_W, x_E;
        std::vector<double> y_W, y_E;
        std::vector<double> z_W, z_E;
        unsigned int i,j,cent = 0;
        for(int l = 0; l < 8; l = l +2){
            if(Iindexes[l+1]%2==0){
                aux_index[0] = Iindexes[l];
                aux_index[1] = Iindexes[l+1]/2;
                subd_index[0] = aux_index;
                aux_index[0] = Iindexes[l] + 1;
                aux_index[1] = Iindexes[l+1]/2;
                subd_index[1] = aux_index;
                i = 0; j = 0;
                sprintf(aux_fname,"Output/Subdomains/X_%d%d.txt",subd_index[0][0],subd_index[0][1]);
                x_W = Read_File(aux_fname);
                sprintf(aux_fname,"Output/Subdomains/Y_%d%d.txt",subd_index[0][0],subd_index[0][1]);
                y_W = Read_File(aux_fname);
                sprintf(aux_fname,"Output/Subdomains/Sol_%d%d.txt",subd_index[0][0],subd_index[0][1]);
                z_W = Read_File(aux_fname);
                sprintf(aux_fname,"Output/Subdomains/X_%d%d.txt",subd_index[1][0],subd_index[1][1]);
                x_E = Read_File(aux_fname);
                sprintf(aux_fname,"Output/Subdomains/Y_%d%d.txt",subd_index[1][0],subd_index[1][1]);
                y_E = Read_File(aux_fname);
                sprintf(aux_fname,"Output/Subdomains/Sol_%d%d.txt",subd_index[1][0],subd_index[1][1]);
                z_E = Read_File(aux_fname);
                for(j = cent; j < y_W.size(); j++){
                    fprintf(oy_file,"%e\n",y_W[j]);
                    for(i = 0; i < x_W.size(); i++){
                        fprintf(osol_file," %e",z_W[x_W.size()*j + i]);
                    }
                    for(i = 1; i < x_E.size(); i++){
                        fprintf(osol_file," %e",z_E[x_W.size()*j + i]);
                    }
                    fprintf(osol_file,"\n");
                }
                cent ++;
            }
        }
        fclose(osol_file);
        fclose(oy_file);
        for(i = 0; i < x_W.size(); i++){
                fprintf(ox_file,"%e\n",x_W[i]);
        }
        for(i = 1; i < x_E.size(); i++){
                fprintf(ox_file,"%e\n",x_E[i]);
        }
        fclose(ox_file);
    }
    return Recieve_Stencil_Data();
}
void PDDSparseGM::Send_G_B(void){
    work_control[0] = (int) G.size();
    work_control[1] = (int) B.size();
    MPI_Send(work_control, 2, MPI_INT, server, REQUEST_NODE, world);
    //G_i is sent
    int *G_i_aux;
    G_i_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_i_aux[i] = G_i[i];
    G_i.resize(0);
    MPI_Send(G_i_aux, work_control[0], MPI_INT, server, TAG_Gi, world);
    delete G_i_aux;
    //G_j is sent
    int *G_j_aux;
    G_j_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_j_aux[i] = G_j[i];
    G_j.resize(0);
    MPI_Send(G_j_aux, work_control[0], MPI_INT, server, TAG_Gj, world);
    delete G_j_aux;
    //G_val is sent
    double *G_aux;
    G_aux = new double[work_control[0]];
    for(int i = 0; i < work_control[0]; i++){
         G_aux[i] = G[i];
    }
    G.resize(0);
    MPI_Send(G_aux, work_control[0], MPI_DOUBLE, server,TAG_Gval, world);
    for(int i = 0; i < work_control[0]; i++){
         G_aux[i] = G_var[i];
    }
    MPI_Send(G_aux, work_control[0], MPI_DOUBLE, server,TAG_Gvar, world);
    delete G_aux;
    //B_i is sent
    int *B_i_aux;
    B_i_aux = new int[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_i_aux[i] = B_i[i];
    B_i.resize(0);
    MPI_Send(B_i_aux, work_control[1], MPI_INT, server, TAG_Bi, world);
    delete B_i_aux;
    //B is sent
    double *B_aux;
    B_aux = new double[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_aux[i] = B[i];
    B.resize(0);
    MPI_Send(B_aux, work_control[1], MPI_DOUBLE, server, TAG_Bval, world);
    for(int i = 0; i < work_control[1]; i++) B_aux[i] = B_var[i];
    B.resize(0);
    MPI_Send(B_aux, work_control[1], MPI_DOUBLE, server, TAG_Bvar, world);
    delete B_aux;
    FILE *dfile;
    dfile = fopen(debug_fname,"a");
    fprintf(dfile,"Processor %d ended sending its contribution to G and B\n", workerid);
    fclose(dfile);
    //MPI_Comm_free(&workers);
}
void PDDSparseGM::Send_G_B_2(void){
    work_control[0] = (int) G.size();
    work_control[1] = (int) B.size();
    MPI_Send(work_control, 2, MPI_INT, server, REQUEST_NODE, world);
    //G_i is sent
    int *G_i_aux;
    G_i_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_i_aux[i] = G_i[i];
    G_i.resize(0);
    MPI_Send(G_i_aux, work_control[0], MPI_INT, server, TAG_Gi, world);
    delete G_i_aux;
    //G_j is sent
    int *G_j_aux;
    G_j_aux = new int[work_control[0]];
    for(int i = 0; i < work_control[0]; i++) G_j_aux[i] = G_j[i];
    G_j.resize(0);
    MPI_Send(G_j_aux, work_control[0], MPI_INT, server, TAG_Gj, world);
    delete G_j_aux;
    //G_val is sent
    double *G_aux;
    G_aux = new double[work_control[0]];
    for(int i = 0; i < work_control[0]; i++){
         G_aux[i] = G[i];
    }
    G.resize(0);
    MPI_Send(G_aux, work_control[0], MPI_DOUBLE, server,TAG_Gval, world);
    for(int i = 0; i < work_control[0]; i++){
         G_aux[i] = G_var[i];
    }
    G_var.resize(0);
    MPI_Send(G_aux, work_control[0], MPI_DOUBLE, server,TAG_Gvar, world);
    for(int i = 0; i < work_control[0]; i++){
         G_aux[i] = G_CT[i];
    }
    G_CT.resize(0);
    MPI_Send(G_aux, work_control[0], MPI_DOUBLE, server,TAG_GCT, world);
    delete G_aux;
    //B_i is sent
    int *B_i_aux;
    B_i_aux = new int[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_i_aux[i] = B_i[i];
    B_i.resize(0);
    MPI_Send(B_i_aux, work_control[1], MPI_INT, server, TAG_Bi, world);
    delete B_i_aux;
    //B is sent
    double *B_aux;
    B_aux = new double[work_control[1]];
    for(int i = 0; i < work_control[1]; i++) B_aux[i] = B[i];
    B.resize(0);
    MPI_Send(B_aux, work_control[1], MPI_DOUBLE, server, TAG_Bval, world);
    for(int i = 0; i < work_control[1]; i++) B_aux[i] = B_var[i];
    B_var.resize(0);
    MPI_Send(B_aux, work_control[1], MPI_DOUBLE, server, TAG_Bvar, world);
    for(int i = 0; i < work_control[1]; i++) B_aux[i] = B_CT[i];
    B_CT.resize(0);
    MPI_Send(B_aux, work_control[1], MPI_DOUBLE, server, TAG_BCT, world);
    delete B_aux;
    FILE *dfile;
    dfile = fopen(debug_fname,"a");
    fprintf(dfile,"Processor %d ended sending its contribution to G and B\n", workerid);
    fclose(dfile);
    //MPI_Comm_free(&workers);
}
void PDDSparseGM::Receive_G_B(void){
    //It ends worker process solving
    MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
    REQUEST_NODE, world, &status);
    work_control[0] = 0;
    work_control[1] = 0;
    MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_NODE, world);
    //The server accumulates the values of G and B
    MPI_Recv(work_control, 2, MPI_INT, status.MPI_SOURCE, REQUEST_NODE, world, &status);
    G.resize(work_control[0]);
    G_i.resize(work_control[0]);
    G_j.resize(work_control[0]);
    //G_i is received
    int *G_i_aux, *G_j_aux;
    double *G_aux;
    G_i_aux = new int[work_control[0]];
    MPI_Recv(G_i_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gi, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_i[i] = G_i_aux[i];
    delete G_i_aux;
    //G_j is received
    G_j_aux = new int[work_control[0]];
    MPI_Recv(G_j_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gj, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_j[i] = G_j_aux[i];
    delete G_j_aux;
    //G is received
    G_aux = new double[work_control[0]];
    MPI_Recv(G_aux, work_control[0], MPI_DOUBLE, status.MPI_SOURCE, TAG_Gval, world, &status);
    for(int i = 0; i < work_control[0]; i++) G[i] = G_aux[i];
    //Triplets vector is fullfilled
    for(int j = 0; j < work_control[0]; j++){
        T_vec_G.push_back(T(G_i[j], G_j[j], G[j]));
    }
    //G_var is received 
    MPI_Recv(G_aux, work_control[0], MPI_DOUBLE, status.MPI_SOURCE, TAG_Gvar, world, &status);
    for(int i = 0; i < work_control[0]; i++) G[i] = G_aux[i];
    for(int j = 0; j < work_control[0]; j++){
        T_vec_Gvar.push_back(T(G_i[j], G_j[j], G[j]));
    }
    delete G_aux;
    G.resize(0);
    G_i.resize(0);
    G_j.resize(0);
    B.resize(work_control[1]);
    B_i.resize(work_control[1]);
    double *B_aux;
    int *B_i_aux;
    //B_i is received
    B_i_aux = new int[work_control[1]];
    MPI_Recv(B_i_aux, work_control[1], MPI_INT, status.MPI_SOURCE, TAG_Bi, world, &status);
    for(int i = 0; i< work_control[1]; i++) B_i[i] = B_i_aux[i];
    delete B_i_aux;
    //B is received
    B_aux = new double[work_control[1]];
    MPI_Recv(B_aux, work_control[1], MPI_DOUBLE, status.MPI_SOURCE, TAG_Bval, world, &status);
    for(int i = 0; i < work_control[1]; i++) B[i] = B_aux[i];
    for(int j = 0; j < work_control[1]; j++){
        T_vec_B.push_back(T(B_i[j], 0, B[j]));
    }
    MPI_Recv(B_aux, work_control[1], MPI_DOUBLE, status.MPI_SOURCE, TAG_Bvar, world, &status);
    for(int i = 0; i < work_control[1]; i++) B[i] = B_aux[i];
    for(int j = 0; j < work_control[1]; j++){
        T_vec_Bvar.push_back(T(B_i[j], 0, B[j]));
    }
    delete B_aux;
}
void PDDSparseGM::Receive_G_B_2(void){
    //It ends worker process solving
    MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, 
    REQUEST_NODE, world, &status);
    work_control[0] = 0;
    work_control[1] = 0;
    MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_NODE, world);
    //The server accumulates the values of G and B
    MPI_Recv(work_control, 2, MPI_INT, status.MPI_SOURCE, REQUEST_NODE, world, &status);
    G.resize(work_control[0]);
    G_i.resize(work_control[0]);
    G_j.resize(work_control[0]);
    //G_i is received
    int *G_i_aux, *G_j_aux;
    double *G_aux;
    G_i_aux = new int[work_control[0]];
    MPI_Recv(G_i_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gi, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_i[i] = G_i_aux[i];
    delete G_i_aux;
    //G_j is received
    G_j_aux = new int[work_control[0]];
    MPI_Recv(G_j_aux, work_control[0], MPI_INT, status.MPI_SOURCE, TAG_Gj, world, &status);
    for(int i = 0; i < work_control[0]; i++) G_j[i] = G_j_aux[i];
    delete G_j_aux;
    //G is received
    G_aux = new double[work_control[0]];
    MPI_Recv(G_aux, work_control[0], MPI_DOUBLE, status.MPI_SOURCE, TAG_Gval, world, &status);
    for(int i = 0; i < work_control[0]; i++) G[i] = G_aux[i];
    //Triplets vector is fullfilled
    for(int j = 0; j < work_control[0]; j++){
        T_vec_G.push_back(T(G_i[j], G_j[j], G[j]));
    }
    //G_var is received 
    MPI_Recv(G_aux, work_control[0], MPI_DOUBLE, status.MPI_SOURCE, TAG_Gvar, world, &status);
    for(int i = 0; i < work_control[0]; i++) G[i] = G_aux[i];
    for(int j = 0; j < work_control[0]; j++){
        T_vec_Gvar.push_back(T(G_i[j], G_j[j], G[j]));
    }
    //G_CT is received 
     MPI_Recv(G_aux, work_control[0], MPI_DOUBLE, status.MPI_SOURCE, TAG_GCT, world, &status);
    for(int i = 0; i < work_control[0]; i++) G[i] = G_aux[i];
    for(int j = 0; j < work_control[0]; j++){
        T_vec_GCT.push_back(T(G_i[j], G_j[j], G[j]));
    }
    delete G_aux;
    G.resize(0);
    G_i.resize(0);
    G_j.resize(0);
    B.resize(work_control[1]);
    B_i.resize(work_control[1]);
    double *B_aux;
    int *B_i_aux;
    //B_i is received
    B_i_aux = new int[work_control[1]];
    MPI_Recv(B_i_aux, work_control[1], MPI_INT, status.MPI_SOURCE, TAG_Bi, world, &status);
    for(int i = 0; i< work_control[1]; i++) B_i[i] = B_i_aux[i];
    delete B_i_aux;
    //B is received
    B_aux = new double[work_control[1]];
    MPI_Recv(B_aux, work_control[1], MPI_DOUBLE, status.MPI_SOURCE, TAG_Bval, world, &status);
    for(int i = 0; i < work_control[1]; i++) B[i] = B_aux[i];
    for(int j = 0; j < work_control[1]; j++){
        T_vec_B.push_back(T(B_i[j], 0, B[j]));
    }
    //B_var is received
    MPI_Recv(B_aux, work_control[1], MPI_DOUBLE, status.MPI_SOURCE, TAG_Bvar, world, &status);
    for(int i = 0; i < work_control[1]; i++) B[i] = B_aux[i];
    for(int j = 0; j < work_control[1]; j++){
        T_vec_Bvar.push_back(T(B_i[j], 0, B[j]));
    }
    //B_CT is received
    MPI_Recv(B_aux, work_control[1], MPI_DOUBLE, status.MPI_SOURCE, TAG_BCT, world, &status);
    for(int i = 0; i < work_control[1]; i++) B[i] = B_aux[i];
    for(int j = 0; j < work_control[1]; j++){
        T_vec_BCT.push_back(T(B_i[j], 0, B[j]));
    }
    delete B_aux;
}
Eigen::VectorXd PDDSparseGM::Node_Position(int node_index){
    Eigen::VectorXd output;
    output.resize(2);
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        if((*it).In(node_index)){
            output = (*it).Node_Position(node_index);
        }
    }
    return output;
}
void PDDSparseGM::Compute_Solution(BVP bvp){
    Eigen::SparseMatrix<double> G_sparse, I_sparse, B_sparse;
    unsigned int size = (unsigned int) nNodes;
    
    FILE *ofile;
    G_sparse.resize(size, size);
    G_sparse.setFromTriplets(T_vec_Gvar.begin(),T_vec_Gvar.end());
    T_vec_Gvar.resize(0);
    ofile = fopen("Output/Debug/G_var.csv","w");
    fprintf(ofile, "Gvar_i,Gvar_j,Gvar_ij\n");
    for (int k=0; k<G_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    G_sparse.resize(size, size);
    G_sparse.setFromTriplets(T_vec_G.begin(),T_vec_G.end());
    T_vec_G.resize(0);
    I_sparse.resize(size, size);
    I_sparse.setIdentity();
    G_sparse+=I_sparse;
    I_sparse.resize(0,0);
    ofile = fopen("Output/Debug/G.csv","w");
    fprintf(ofile, "G_i,G_j,G_ij\n");
    for (int k=0; k<G_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    B_sparse.resize(size, 1);
    B_sparse.setFromTriplets(T_vec_Bvar.begin(), T_vec_Bvar.end());
    T_vec_Bvar.resize(0);
    ofile = fopen("Output/Debug/B_var.csv","w");
    fprintf(ofile, "B_i,B_j,B_ij\n");
    for (int k=0; k<B_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    B_sparse.resize(size, 1);
    B_sparse.setFromTriplets(T_vec_B.begin(), T_vec_B.end());
    T_vec_B.resize(0);
    ofile = fopen("Output/Debug/B.csv","w");
    fprintf(ofile, "B_i,B_j,B_ij\n");
    for (int k=0; k<B_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    Eigen::VectorXd ud,Bd,guess;
    Bd = Eigen::VectorXd(B_sparse);
    guess = Bd*0;
    B_sparse.resize(0,0);
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_LU;
    //solver_LU.compute(G_sparse);
    //solver_LU.factorize(G_sparse);
    //ud = solver_LU.solve(Bd);
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_IT;
    solver_IT.compute(G_sparse);
    ud = solver_IT.solveWithGuess(Bd,guess);
    ofile = fopen("Output/solution.csv","w");
    fprintf(ofile,"Knot_index,x,y,sol_analytic,sol_PDDS,err,rerr\n");
    double sol, err, rerr;
    for(int i = 0; i< (int)ud.size(); i++){
        sol = bvp.u.Value(Node_Position(i),T_start);
        err = ud(i) - sol;
        rerr = fabs(err)/sol;
        fprintf(ofile,"%d,%e,%e,%e,%e,%e,%e\n",i,Node_Position(i)[0],Node_Position(i)[1],
        sol,ud(i),err,rerr);
    }
    fclose(ofile);
}
void PDDSparseGM::Compute_Solution_2(BVP bvp){
    Eigen::SparseMatrix<double> G_sparse, I_sparse, B_sparse;
    unsigned int size = (unsigned int) nNodes;
    FILE *ofile;
    G_sparse.resize(size, size);
    G_sparse.setFromTriplets(T_vec_Gvar.begin(),T_vec_Gvar.end());
    T_vec_Gvar.resize(0);
    G_sparse_var = G_sparse;
    ofile = fopen("Output/Debug/G_var.csv","w");
    fprintf(ofile, "Gvar_i,Gvar_j,Gvar_ij\n");
    for (int k=0; k<G_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    G_sparse.resize(size, size);
    G_sparse.setFromTriplets(T_vec_G.begin(),T_vec_G.end());
    T_vec_G.resize(0);
    I_sparse.resize(size, size);
    I_sparse.setIdentity();
    G_sparse+=I_sparse;
    ofile = fopen("Output/Debug/G.csv","w");
    fprintf(ofile, "G_i,G_j,G_ij\n");
    for (int k=0; k<G_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    B_sparse.resize(size, 1);
    B_sparse.setFromTriplets(T_vec_Bvar.begin(), T_vec_Bvar.end());
    T_vec_Bvar.resize(0);
    ofile = fopen("Output/Debug/B_var.csv","w");
    fprintf(ofile, "B_i,B_j,B_ij\n");
    for (int k=0; k<B_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    B_sparse.resize(size, 1);
    B_sparse.setFromTriplets(T_vec_B.begin(), T_vec_B.end());
    T_vec_B.resize(0);
    ofile = fopen("Output/Debug/B.csv","w");
    fprintf(ofile, "B_i,B_j,B_ij\n");
    for (int k=0; k<B_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    Eigen::VectorXd ud,ud_CT,Bd,guess;
    Bd = Eigen::VectorXd(B_sparse);
    guess = Bd*0;
    B_sparse.resize(0,0);
    /*Eigen::SparseLU<Eigen::SparseMatrix<double> > solver_LU;
    solver_LU.compute(G_sparse);
    solver_LU.factorize(G_sparse);
    ud = solver_LU.solve(Bd);*/
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_IT;
    solver_IT.compute(G_sparse);
    ud = solver_IT.solveWithGuess(Bd,guess);
    knot_solutions = ud;
    G_sparse.resize(size, size);
    G_sparse.setFromTriplets(T_vec_GCT.begin(),T_vec_GCT.end());
    G_sparse+=I_sparse;
    I_sparse.resize(0,0);
    ofile = fopen("Output/Debug/G_CT.csv","w");
    fprintf(ofile, "G_i,G_j,G_ij\n");
    for (int k=0; k<G_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    T_vec_GCT.resize(0);
    B_sparse.resize(size, 1);
    B_sparse.setFromTriplets(T_vec_BCT.begin(), T_vec_BCT.end());
    T_vec_BCT.resize(0);
    ofile = fopen("Output/Debug/B_CT.csv","w");
    fprintf(ofile, "B_i,B_j,B_ij\n");
    for (int k=0; k<B_sparse.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(B_sparse,k); it; ++it)
        {
            fprintf(ofile,"%ld,%ld,%e\n",it.row(),it.col(),it.value());
        }
    fclose(ofile);
    solver_IT.compute(G_sparse);
    Bd = Eigen::VectorXd(B_sparse);
    ud_CT = solver_IT.solveWithGuess(Bd,guess);
    bias = (ud-ud_CT).cwiseAbs();
    ofile = fopen("Output/solution.csv","w");
    fprintf(ofile,"Knot_index,x,y,sol_analytic,sol_PDDS,err,rerr,sbias\n");
    double sol, err, rerr;
    for(int i = 0; i< (int)ud.size(); i++){
        sol = bvp.u.Value(Node_Position(i),T_start);
        err = ud(i) - sol;
        rerr = fabs(err)/sol;
        fprintf(ofile,"%d,%e,%e,%e,%e,%e,%e,%e\n",i,Node_Position(i)[0],Node_Position(i)[1],
        sol,ud(i),err,rerr,bias(i));
    }
    fclose(ofile);
}
void PDDSparseGM::Print_Interface(void){
    if(myid==server){
    int cent = 0;
    for(std::vector<Interface>::iterator it = interfaces.begin();
    it != interfaces.end();
    it++){
        (*it).Print(cent);
        cent ++;
    }
    }
}
void PDDSparseGM::Print_Problem(void){
    FILE *dfile;
    dfile = fopen(debug_fname,"a");
    fprintf(dfile,"**********PDDSparse Solver*************\n");
    fprintf(dfile,"***************Ver 1.0******************\n");
    fprintf(dfile,"***********Jorge Moron Vidal************\n");
    fprintf(dfile,"**********jmoron@math.uc3m.es***********\n");
    fprintf(dfile,"h0 = %e \t T_start = %0.2f\n",h0, T_start);
    fprintf(dfile,"Number of trayectories per knot %d\n",N);
    fprintf(dfile,"Number of trayectories per job %d\n", N_job);
    fprintf(dfile,"Number of subdomains/interfaces: %d x %d \n", iN[0], iN[1]);
    fprintf(dfile,"Number of knots per interface: %d\n",nN[0]);
    fprintf(dfile,"Position of South-West corner of the domain: (%.2f,%.2f)\n",SW(0), SW(1));
    fprintf(dfile,"Position of North-East corner of the domain: (%.2f,%.2f)\n",NE(0), NE(1));
    fclose(dfile);
}
void PDDSparseGM::Read_Solution(void){
   std::ifstream solution_file("Output/solution.csv");
   std::string line, aux_string;
   std::vector<double> csv_row;
   uint16_t aux_index;
   if(solution_file.is_open()){
       std::cout << "Reading solution file" << std::endl;
       csv_row.clear();
       getline(solution_file, line);
       csv_row.resize(std::count(line.begin(),line.end(),',') +1);
       aux_string.clear(); 
       while(getline(solution_file, line)){
           aux_index = 0;
           for(std::string::iterator character = line.begin();
           character != line.end();
           character ++){
               if (*character == ','){
                   csv_row[aux_index] = atof(aux_string.c_str());
                   aux_string.clear();
                   aux_index ++;
               } else {
                   aux_string += *character;
               }
           }
           csv_row[aux_index] = atof(aux_string.c_str());
           aux_string.clear();
           for(std::vector<Interface>::iterator interface = interfaces.begin();
               interface != interfaces.end();
               interface ++){
                if( (*interface).In((int) csv_row[0])){
                    for(uint32_t interface_index = 0;
                        interface_index < (*interface).index.size();
                        interface_index ++){
                        if((*interface).index[interface_index] == (int) csv_row[0]){
                            (*interface).solution[interface_index] = csv_row[4];
                            /* printf("Node %d is in interface [%d %d] with solution %f \n",
                            (int) csv_row[0],(*interface).label[0],(*interface).label[1],
                            csv_row[4]);*/
                        }
                    }
                }
            }
        }

   } else {
       std::cout << "Output/solution.csv couldn't be opened" << std::endl;
   }
   solution_file.close();
}
void PDDSparseGM::Fullfill_subdomains(BVP bvp){
    subdomains.clear();
    Subdomain subdomain_aux;
    std::vector<int> label_aux;
    std::map<direction, Interface> interfaces_map_aux;
    label_aux.resize(2);
    for(int subd_y = 0; subd_y < iN[1]; subd_y ++){
        for(int subd_x = 0; subd_x < iN[0]; subd_x++){
            label_aux[0] = subd_x;
            interfaces_map_aux.clear();
            label_aux[1] = subd_y*2;
            if(interface_map.count(label_aux) > 0 ){
                interfaces_map_aux[East] = interfaces[interface_map[label_aux]];
            }
            label_aux[0] = subd_x-1;
            label_aux[1] = subd_y*2;
            if(interface_map.count(label_aux) > 0 ){
                interfaces_map_aux[West] = interfaces[interface_map[label_aux]];
            }
            label_aux[0] = subd_x;
            label_aux[1] = subd_y*2 +1;
            if(interface_map.count(label_aux) > 0 ){
                interfaces_map_aux[North] = interfaces[interface_map[label_aux]];
            }
            label_aux[0] = subd_x;
            label_aux[1] = subd_y*2 -1;
            if(interface_map.count(label_aux) > 0 ){
                interfaces_map_aux[South] = interfaces[interface_map[label_aux]];
            }
            label_aux[0] = subd_x;
            label_aux[1] = subd_y;
            subdomain_aux.Init(label_aux, interfaces_map_aux, parameters, bvp);
            subdomains.push_back(subdomain_aux);
        }
    }
}
void PDDSparseGM::Solve_Subdomains(BVP bvp){
    if(myid == server){
        Read_Solution();
        Fullfill_subdomains(bvp);
        for(std::vector<Subdomain>::iterator it_subdomain = subdomains.begin();
            it_subdomain != subdomains.end();
            it_subdomain ++){
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, ASK_SERVER, world, &status);
            work_control[0] = SEND_WORK;
            MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
            (*it_subdomain).Send_To_Worker(status,world);
        }
        for(int i = 0; i<server; i++){
            MPI_Recv(work_control, 2, MPI_INT, MPI_ANY_SOURCE, ASK_SERVER, world, &status);
            work_control[0] = END_WORKER;
            MPI_Send(work_control, 2, MPI_INT, status.MPI_SOURCE, REPLY_WORKER, world);
        }
    } else {
        Subdomain subdomain;
        do{
            MPI_Send(work_control, 2, MPI_INT, server, ASK_SERVER, world);
            MPI_Recv(work_control, 2, MPI_INT, server, REPLY_WORKER, world, &status);
            if(work_control[0] == SEND_WORK){
                subdomain.Recieve_From_Server(server, world);
                subdomain.Solve(world);
            }
        }while(work_control[0] == SEND_WORK);
    }
    //MPI_Finalize();
}
void PDDSparseGM::Solve_NumVR(BVP bvp, std::vector<double> h_vec, std::vector<int> N_vec){
    bool done=true;
    double start = MPI_Wtime();
    FILE *pFile;
    //if(myid==server) Print_Problem();
    //Compute the PDDSparse Matrix
    if(myid==server){
        system("mv Output/solution.csv Output/solution_nvarred.csv");
        system("mv Output/Debug/B.csv Output/Debug/B_nvarred.csv");
        system("mv Output/Debug/G.csv Output/Debug/G_nvarred.csv");
        system("mv Output/Debug/B_CT.csv Output/Debug/B_CT_nvarred.csv");
        system("mv Output/Debug/G_CT.csv Output/Debug/G_CT_nvarred.csv");
        system("mv Output/Debug/B_var.csv Output/Debug/B_var_nvarred.csv");
        system("mv Output/Debug/G_var.csv Output/Debug/G_var_nvarred.csv");  
        PDDSJob auxjob;
        //Job list is fullfilled
        for(int knot_index = 0; knot_index < nNodes; knot_index++){
            auxjob.index[0] = knot_index;
            //auxjob.N[0] = N_vec[knot_index]; auxjob.N[1] = 10000;
            auxjob.N[0] = N_vec[knot_index]; auxjob.N[1] = N_job; auxjob.h[0] = h_vec[knot_index];
            for(int jobs_per_knot = 0; jobs_per_knot < N_vec[knot_index]/N_job; jobs_per_knot ++){
                work_control[1] = 1;
                Send_Node_Loop(auxjob);
            }
        }
        //G and B are received from the workers
        for(int process = 0; process < server; process++){
             Receive_G_B();
        }
        Compute_Solution_2(bvp);
        double end = MPI_Wtime();
    } else {
        //Start and end time for the node
        double knot_start, knot_end;
        //G and B storage vectors for each node
        double B_temp;
        PDDSJob auxjob;
        std::vector<double> G_temp;
        std::vector<int>  G_j_temp;
        //Stencil
        Stencil stencil;
        VectorFunction grad;
        std::stringstream dirname;
        dirname << "Input/Patch_"<<myid;
        GMSolver solver(bvp ,parameters, h0, (unsigned int) myid +1);
        c2 = pow(fac*(1.0/nN[0]),2.0);
        G_i.clear();G_j.clear(); G.clear(); G_CT.clear(); G_var.clear();
        B.clear();B_CT.clear();B_i.clear();B_var.clear();
        do{
            done = Receive_Node(auxjob);
            if(!done){
                knot_start = MPI_Wtime();
                B_temp = 0.0;
                G_j_temp.clear();
                G_temp.clear();
                stencil = Recieve_Stencil_Data_Loop();
                stencil.Compute_ipsi(bvp,c2,debug_fname);
                grad.Init(2,dirname.str());
                solver.h = auxjob.h[0];
                if(stencil.Is_Interior()){ 
                    printf("Knot %d is interior\n",auxjob.index[0]);
                    solver.Solve(position,c2,stencil,G_j_temp,G_temp,B_temp,auxjob.N[0], auxjob.N[1],grad);
                }else{
                    printf("Knot %d is exterior\n",auxjob.index[0]);
                    solver.Solve_mix(position,c2,1.0,stencil,G_j_temp,G_temp,B_temp,auxjob.N[0], auxjob.N[1],grad);
                }
                B.push_back(B_temp);
                B_i.push_back(auxjob.index[0]);
                B_var.push_back(solver.var_B);
                for(unsigned int k = 0;k < G_temp.size(); k ++){
                    G_i.push_back(auxjob.index[0]);
                    G_j.push_back(G_j_temp[k]);
                    G.push_back(G_temp[k]);
                    G_var.push_back(solver.var_G[k]);
                }
                knot_end = MPI_Wtime();
                pFile = fopen("Output/Debug/Node_debug.csv","a");
                fprintf(pFile,"%d,%f,%f,%f\n",auxjob.index[0],
                solver.var_B,solver.APL,(knot_end-knot_start)/60);
                fclose(pFile);
                pFile = fopen("Output/Debug/G_var.csv","a");
                for(unsigned int i = 0; i < solver.var_G.size(); i++){
                    fprintf(pFile,"%d,%d,%f\n",auxjob.index[0],
                    G_j_temp[i],solver.var_G[i]);
                }
                fclose(pFile);
            }
        }while(!done);
        pFile = fopen(debug_fname,"a");
        fprintf(pFile,"Process %d is done solving nodes \n", myid);
        Send_G_B();
        fprintf(pFile,"Process %d is done sending G \n", myid);
        double end = MPI_Wtime();
        fprintf(pFile,"Process %d used %f  hours \n", myid, (end-start)/3600);
        fclose(pFile);
        MPI_Comm_free(&workers);
    }
    //MPI_Finalize();
}
void PDDSparseGM::Compute_h_N(BVP bvp, double eps, std::vector<double> & h_vec, std::vector<int> & N_vec){
    if(myid==server){
        std::map<direction, std::vector<std::vector <int> > > labels;
        double s_params[4],h_aux,h_max;
        int aux_size;
        std::vector<double> knot_est_variance;
        Eigen::VectorXd X_aux, N_aux,E_P_aux;
        X_aux.resize(2); N_aux.resize(2); E_P_aux.resize(2);
        knot_est_variance.resize(nNodes);
        h_vec.resize(nNodes);
        N_vec.resize(nNodes);
        //var_G*u vector
        FILE *pfile;
        pfile = fopen("Output/Debug/h_and_N.csv","w");
        fprintf(pfile,"Knot,h,N\n");
        for(int i = 0; i < nNodes; i++) knot_est_variance[i] = 0.0;
        for(int k=0; k<G_sparse_var.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(G_sparse_var,k); it; ++it)
            {
                knot_est_variance[it.row()] += knot_solutions[it.col()]*knot_solutions[it.col()]*it.value();
            }
        }
        for(int i = 0; i < nNodes; i++){
            Labels_Stencil(i);
            //Stencil parameters
            for(int i = 0; i < 4; i++) s_params[i] = parameters[i];
            if(labels[South].size()>0){
                //s_params[0] = interfaces[interface_map[labels[South][0]]].position[0][0];
                s_params[1] = interfaces[interface_map[labels[South][0]]].position[0][1];
            }
            if(labels[West].size()>0){
                s_params[0] = interfaces[interface_map[labels[West][0]]].position[0][0];
                //s_params[1] = interfaces[interface_map[labels[West][0]]].position[0][1];
            }
            if(labels[North].size()>0){
                aux_size = (int) interfaces[interface_map[labels[North][labels[North].size()-1]]].position.size();
                //s_params[2] = interfaces[interface_map[labels[North][labels[North].size()-1]]].position[aux_size - 1][0];
                s_params[3] = interfaces[interface_map[labels[North][labels[North].size()-1]]].position[aux_size - 1][1];
            }
            if(labels[East].size()>0){
                aux_size = (int) interfaces[interface_map[labels[East][labels[East].size()-1]]].position.size();
                s_params[2] = interfaces[interface_map[labels[East][labels[East].size()-1]]].position[aux_size - 1][0];
                //s_params[3] = interfaces[interface_map[labels[East][labels[East].size()-1]]].position[aux_size - 1][1];
            }
            h_max = pow(bvp.boundary.Dist(s_params, X_aux, E_P_aux, N_aux)/(
            (N_aux.transpose()*bvp.sigma.Value(X_aux,0.0)).norm()*2*0.5826),2.0);
            h_aux = 0.3333*eps*h0/bias(i);
            h_vec[i] = std::min(h_max,h_aux);
            N_vec[i] = std::max((int)(9.0*knot_est_variance[i]/(eps*eps)),10000);
            fprintf(pfile,"%d,%e,%d\n",i,h_vec[i],N_vec[i]);
        }
    }
}
std::vector<double> PDDSparseGM::Read_File(char fname[256]){
    std::vector<double> output;
    char *line_buf = NULL;
    size_t line_buf_size = 0;
    int line_count = 0;
    ssize_t line_size;
    FILE *fp = fopen(fname, "r");
    if (!fp)
    {
        fprintf(stderr, "Error opening file '%s'\n", fname);
        return output;
    }
    /* Get the first line of the file. */
    line_size = getline(&line_buf, &line_buf_size, fp);
    /* Loop through until we are done with the file. */
    while (line_size > 0)
    {
        /* Increment our line count */
        line_count++;
        /* Show the line details */
        output.push_back(atof(line_buf));
        /* Get the next line */
        line_size = getline(&line_buf, &line_buf_size, fp);
    }
   /* Free the allocated line buffer */
   free(line_buf);
   line_buf = NULL;
   fclose(fp);
   return output;
}