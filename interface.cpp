#include "interface.hpp"

Interface::Interface(void){
	interior = false;
	direction = -1;
}

void Interface::Init(Eigen::VectorXf start, 
                   Eigen::VectorXf end,
                   std::vector<int> i_index,
                   std::vector<int> subd_index,
                   std::vector <int>n_index,
                   int direc,
                   Eigen::SparseMatrix<float> Psim,
                   bool inter,
                   bool chebyshev,
                   float tol,
                   float discretization){

	std::srand(n_index[0]);
	Eigen::VectorXf X, increment;
	Eigen::SparseMatrix<float> H;
	H.resize(Psim.rows(),Psim.cols());
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
							   tol,
							   discretization,
							   (unsigned int) std::rand(),
							   n_index[i],
							   i_index,
							   subd_index,
							   H,
							   Psim.row(n_index[i]));
			X += increment;
		}
	}

	Interface_it = node_array.begin();
}

void Interface::Set_Stencil(std::vector<Interface> Sten, float *parameters){

	for(std::vector<Interface>::iterator Stencil_it = Sten.begin(); 
		Stencil_it != Sten.end(); 
		Stencil_it ++){
		
		for(std::vector<Node>::iterator Node_it = (*Stencil_it).node_array.begin();
			Node_it != (*Stencil_it).node_array.end();
			Node_it ++){

			stencil_position.push_back((*Node_it).x0);
			stencil_j.push_back((*Node_it).i_node);

		}

	}
}

void Interface::Print_Interface(void){
	std::cout << "\nInterface "; 
	for(unsigned int i = 0; i < interface_index.size(); i++) std::cout << interface_index[i] << ' ';
	std::cout << '\n';
	std::cout << "Interior = " << interior << '\n';
	std::cout <<  "Is part of subdomains: ";
	for(unsigned int i = 0; i < subdomain_index.size(); i++) std::cout << subdomain_index[i] << ' ';
	std::cout << '\n';
	std::cout << "And stores the nodes" << '\n';
	for(unsigned int i = 0; i < node_array.size(); i++) std::cout << "Node " << node_array[i].i_node 
	<< " : " << node_array[i].x0(0) << " " << node_array[i].x0(1) << " ";
	std::cout << '\n';
}