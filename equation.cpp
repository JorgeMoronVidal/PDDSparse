#include "equation.hpp"
#include <math.h>
#include <iostream>

float   Equation_c(Eigen::VectorXf X, Eigen::VectorXf N){
    return -(pow(X(0),2.0)+pow(X(1),2.0))*exp(-X(0));
}
float   Equation_f(Eigen::VectorXf X, Eigen::VectorXf N){
    float aux = (pow(X(0),2.0)+pow(X(1),2.0))*exp(-X(0));
    //In order to enforce the solution to be positive
    //aa = 1.5f, bb =2.0f, C = 2.1f
    //res  = res= (aa^2+aux).*cos(aa*x) + (bb^2+aux).*sin(bb*y) + aa*sin(aa*x) + C*aux;
    return (2.25f + aux)*cos(1.5f*X(0)) + (4.0f + aux)*
    sin(2.0f*X(1)) + 1.5f*sin(1.5f*X(0)) + 2.1f*aux;
}
  float Equation_g(Eigen::VectorXf X, Eigen::VectorXf N){
    return Equation_u(X,N);
}
  Eigen::VectorXf Equation_b(Eigen::VectorXf X, Eigen::VectorXf N){
    Eigen::VectorXf b;
    b = Eigen::VectorXf::Constant(2,1.0f);
    return b;
}
  Eigen::VectorXf Equation_F(Eigen::VectorXf X, Eigen::VectorXf N){
    Eigen::VectorXf F(2);
    F = Eigen::VectorXf::Zero(2);
    return F;
}
  Eigen::MatrixXf Equation_sigma(Eigen::VectorXf X, Eigen::VectorXf N){
    Eigen::MatrixXf sigma(2,2);
    sigma << 1.41421356237f ,0.0f,
            0.0f, 1.41421356237f;
    return sigma;
}
  float Equation_u(Eigen::VectorXf X, Eigen::VectorXf N){
    // res= cos(aa*x) + sin(bb*y) + C;
    return cos(1.5f * X(0)) + sin(2.0f*X(1)) + 2.1f;
}

  float Equation_RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    float r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return sqrt(r2 + c2);
}

float Equation_d2dx2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    float r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return (pow(x(1)-xj(1),2)+c2)/pow(sqrt(r2+c2), 3);
}

float Equation_d2dy2RBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    float r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return (pow(x(0)-xj(0),2)+c2)/pow(sqrt(r2+c2), 3);
}

float Equation_uRBF(Eigen::VectorXf x , Eigen::VectorXf xj, float c2){
    return Equation_d2dx2RBF(x,xj,c2) + Equation_d2dy2RBF(x,xj,c2) + 2.0;
}