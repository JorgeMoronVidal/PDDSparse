#include <iostream>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseLU>
#include <vector>
#include <iterator>

int main(){
	typedef Eigen::Triplet<float,int> T;
	Eigen::SparseMatrix<float> H, b, I, H_inv, Psi;
	std::vector<float> values;
	std::vector<int> i,j;
	std::vector<T> T_vec_H, T_vec_b;
	H.resize(4,4);
	I.resize(4,4);
	Psi.resize(4,4);
	I.setIdentity();
	b.resize(4,1);

	i.push_back(0); j.push_back(0); values.push_back(2.0f);
	i.push_back(2); j.push_back(0); values.push_back(1.0f);
	i.push_back(1); j.push_back(1); values.push_back(3.0f);
	i.push_back(1); j.push_back(3); values.push_back(4.5f);
	i.push_back(2); j.push_back(2); values.push_back(5.5f);
	i.push_back(3); j.push_back(1); values.push_back(6.5f);

	for(int k = 0; k < i.size(); k ++){
		T_vec_H.push_back(T(i[k],j[k],values[k]));
	}

	for(int k = 0; k < 4; k++){
		T_vec_b.push_back(T(k,0,1.0f));
	}
	Eigen::MatrixXf Matriz;
	Matriz.resize(3,3);
	Matriz(1,2) = 2;
	Matriz(2,1) = 1;
	H.setFromTriplets(T_vec_H.begin(),T_vec_H.end());
	b.setFromTriplets(T_vec_b.begin(),T_vec_b.end());
	Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;
	solver.compute(H);
	H_inv = solver.solve(I); 

}