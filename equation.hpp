#include <eigen3/Eigen/Core>

float Equation_c(Eigen::VectorXf X, Eigen::VectorXf N);
float Equation_f(Eigen::VectorXf X, Eigen::VectorXf N);
float Equation_g(Eigen::VectorXf X, Eigen::VectorXf N);
Eigen::VectorXf Equation_b(Eigen::VectorXf X, Eigen::VectorXf N);
Eigen::VectorXf Equation_F(Eigen::VectorXf X, Eigen::VectorXf N);
Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, Eigen::VectorXf N);
float Equation_u(Eigen::VectorXf X, Eigen::VectorXf N);
float Equation_RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);
float Equation_d2dx2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);
float Equation_d2dy2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);
float Equation_uRBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2);