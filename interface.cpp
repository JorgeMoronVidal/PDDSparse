#include "interface.hpp"

Interface::Interface(void){
	interior = false;
	direction = -1;
}


void Interface::Init (std::vector <Eigen::VectorXf> node_pos,
           std::vector<int>  i_index,
           std::vector<int>  n_index,
           int direc,
           std::vector< int > subd_index,
           bool inter,
           bool chebyshev,
           int N_Steps,
           float discretization){

	interface_index = i_index;
	interior = inter;
	direction = direc;
	node_array.resize(n_index.size());
	subdomain_index = subd_index;

	if(chebyshev){

		std::cout <<"Chebyshev distribution of nodes not defined yet."  << '\n';

	} else {

		for (unsigned int i = 0; i < n_index.size(); i++){

			node_array[i].init(node_pos[n_index[i]],
							   discretization,
							   N_Steps,
							   n_index[i],
							   i_index,
							   subd_index);
		}
	}

	Interface_it = node_array.begin();

}

void Interface::Solve(BVP bvp,
                  gsl_rng *rng,
                  int N_tray,
                  float c2,
                  Stencil stencil,
                  std::vector<float> & G,
                  std::vector<float> & B,
                  std::vector<int> & G_j,
                  std::vector<int> & G_i,
                  std::vector<int> & B_i){

	//Each entry of solution is the solution in a point
	float B_node;
	std::vector<float> G_node;
	std::vector<int> stencil_index;
	stencil.Compute_ipsi(bvp, c2);
	//Loop over all the nodes in the interface
	for(Interface_it = node_array.begin();
		Interface_it != node_array.end();
		Interface_it ++){
		(*Interface_it).Solve_PDDSparse(bvp, rng, stencil.stencil_parameters, stencil.global_parameters,
									N_tray, c2, stencil, stencil_index, G_node,
									B_node);
		#ifdef DEBUG
		std::cout << "Node " << (*Interface_it).i_node << " was solved \n";
		#endif

		for(unsigned int i = 0; i < stencil_index.size(); i++){
			G.push_back(G_node[i]);
			G_i.push_back((*Interface_it).i_node);
			G_j.push_back(stencil_index[i]);
		}
		B_i.push_back((*Interface_it).i_node);
		B.push_back(B_node);
	}
}

void Interface::Test_Interpolator(BVP bvp,
                  gsl_rng *rng,
                  int N_tray,
                  float c2,
                  Stencil stencil){
	//Each entry of solution is the solution in a point
	float B_node;
	std::vector<float> G_node;
	std::vector<int> stencil_index;
	stencil.Compute_ipsi(bvp, c2);
	//Loop over all the nodes in the interface
	for(Interface_it = node_array.begin();
		Interface_it != node_array.end();
		Interface_it ++){
		(*Interface_it).Test_Interpolator(bvp, rng, stencil.stencil_parameters, stencil.global_parameters,N_tray, c2, stencil);
		#ifdef DEBUG
		std::cout << "Node " << (*Interface_it).i_node << " was tested \n";
		#endif

	}
}

void Interface::Test_G(BVP bvp,
                  gsl_rng *rng,
                  int N_tray,
                  float c2,
                  Stencil stencil,
                  std::vector<float> & G,
                  std::vector<float> & B,
                  std::vector<int> & G_j,
                  std::vector<int> & G_i,
                  std::vector<int> & B_i){

	//Each entry of solution is the solution in a point
	float B_node;
	std::vector<float> G_node;
	std::vector<int> stencil_index;
	stencil.Compute_ipsi(bvp, c2);
	//Loop over all the nodes in the interface
	for(Interface_it = node_array.begin();
		Interface_it != node_array.end();
		Interface_it ++){
		(*Interface_it).Test_G(bvp, rng, stencil.stencil_parameters, stencil.global_parameters,
									N_tray, c2, stencil, stencil_index, G_node,
									B_node);
		#ifdef DEBUG
		std::cout << "Node " << (*Interface_it).i_node << " was G_tested \n";
		#endif

		for(unsigned int i = 0; i < stencil_index.size(); i++){
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