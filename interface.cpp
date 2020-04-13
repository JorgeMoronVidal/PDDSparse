#include "interface.hpp"

Interface::Interface(void){
	interior = false;
	direction = -1;
}

void Interface::Init (Eigen::VectorXf start, 
           Eigen::VectorXf end,
           std::vector<int>  i_index,
           std::vector<int>  n_index,
           int direc,
           std::vector< int > subd_index,
           bool inter,
           bool chebyshev,
           int N_Steps,
           float discretization){

	std::srand(n_index[1]);
	Eigen::VectorXf X, increment;
	X = start;
	interface_index = i_index;
	increment = (end-start)/(n_index.size()-1);
	interior = inter;
	direction = direc;
	node_array.resize(n_index.size());
	subdomain_index = subd_index;

	if(chebyshev){

		std::cout <<"Chebyshev distribution of nodes not defined yet."  << '\n';

	} else {

		for (unsigned int i = 0; i < n_index.size(); i++){

			node_array[i].init(X,
							   discretization,
							   N_Steps,
							   n_index[i],
							   i_index,
							   subd_index);
			X += increment;
		}
	}

	Interface_it = node_array.begin();

}

void Interface::Solve(BVP bvp,
				  gsl_rng *rng,
                  float * parameters_stencil,
                  float * parameters_surface,
                  int N_tray,
                  float c2,
                  std::vector< std::vector<float> > & iPsi,
                  std::vector< Eigen::VectorXf > & stencil_position,
                  std::vector< int > & stencil_index,
                  std::vector<float> & G,
                  std::vector<float> & B,
                  std::vector<int> & G_j,
                  std::vector<int> & G_i,
                  std::vector<int> & B_i){

//Each entry of solution is the solution in a point
float B_node;
std::vector<float> G_node;
//Loop over all the nodes in the interface
for(Interface_it = node_array.begin();
	Interface_it != node_array.end();
	Interface_it ++){
	(*Interface_it).Solve_PDDSparse(bvp, rng, parameters_stencil, parameters_surface,
									N_tray, c2, iPsi, stencil_position, stencil_index, G_node,
									B_node);
	#ifdef DEBUG
		std::cout << "Node " << (*Interface_it).i_node << " was solved \n";
	#endif
		
	for(int i = 0; i < (int)stencil_position.size(); i++){
		G.push_back(G_node[i]);
		G_i.push_back((*Interface_it).i_node);
		G_j.push_back(stencil_index[i]);
	}
	B_i.push_back((*Interface_it).i_node);
	B.push_back(B_node);
	}
}
void Interface::Print_Interface(void){
	std::cout << "\nInterface "; 
	for(unsigned int i = 0; i < interface_index.size(); i++) std::cout << interface_index[i] << ' ';
	std::cout << '\n';
	std::cout << "Interior = " << interior << '\n';
	std::cout << "Direction = " << direction << '\n';
	std::cout <<  "Is part of subdomains: ";
	for(unsigned int i = 0; i < subdomain_index.size(); i++) std::cout << subdomain_index[i] << ' ';
	std::cout << '\n';
	std::cout << "And stores the nodes" << '\n';
	for(unsigned int i = 0; i < node_array.size(); i++) std::cout << "Node " << node_array[i].i_node 
	<< " : " << node_array[i].x0(0) << " " << node_array[i].x0(1) << " ";
	std::cout << '\n';
}