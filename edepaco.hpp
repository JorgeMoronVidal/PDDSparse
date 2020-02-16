#include <eigen3/Eigen/Core>


float cpaco(Eigen::VectorXf X, Eigen::VectorXf N);
float fpaco(Eigen::VectorXf X, Eigen::VectorXf N);
float ppaco(Eigen::VectorXf X, Eigen::VectorXf N);
float gpaco(Eigen::VectorXf X, Eigen::VectorXf N);
float psipaco(Eigen::VectorXf X, Eigen::VectorXf N);
float varpsipaco(Eigen::VectorXf X, Eigen::VectorXf N);
bool Mixpaco(Eigen::VectorXf X, Eigen::VectorXf N);
bool Dirichpaco(Eigen::VectorXf X, Eigen::VectorXf N);
bool Robinpaco(Eigen::VectorXf X, Eigen::VectorXf N);
Eigen::VectorXf bpaco(Eigen::VectorXf X, Eigen::VectorXf N);
Eigen::VectorXf Fpaco(Eigen::VectorXf X, Eigen::VectorXf N);
Eigen::MatrixXf sigmapaco(Eigen::VectorXf X, Eigen::VectorXf N);
float upaco(Eigen::VectorXf X, Eigen::VectorXf N);