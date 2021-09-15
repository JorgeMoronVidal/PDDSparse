#include <eigen3/Eigen/Core>
#include <math.h>
#include <iostream>
#define L 10.0
inline double Equation_c(Eigen::VectorXd X, double t){
    return 0.0;
}
inline double Equation_f(Eigen::VectorXd X, double t){
    return 2.0;
}
inline double Equation_g(Eigen::VectorXd X, double t){
    return 0.0;
}
inline Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t){
    return X*0.0;
}
inline Eigen::MatrixXd Equation_sigma(Eigen::VectorXd, double){
    return Eigen::MatrixXd::Identity(2,2) * 1.41421356237;
}
Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t){
    Eigen::VectorXd F(2);
    double aux = 0.0;
    for(int i = 0.0; i < 10.0; i ++){
        for(int j = 0.0; j < 10.0; j++){
            aux += -pow(-1.0,(double)i+j)*M_PI*sin((2.0*i+1.0)*M_PI*X(0)/L)*cos((2.0*j+1.0)*M_PI*X(1)/L)*
            (1.0/((2.0*j+1.0)*((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1))));
        }
    }

    F(0) = aux;
    aux = 0.0;
    for(int i = 0.0; i < 10.0; i ++){
        for(int j = 0.0; j < 10.0; j++){
            aux += -pow(-1.0,(double)i+j)*cos((2.0*i+1.0)*M_PI*X(0)/L)*M_PI*sin((2.0*j+1.0)*M_PI*X(1)/L)*
            (1.0/((2.0*i+1.0)*((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1))));
        }
    }
    F(1) = aux;
    F = -Equation_sigma(X,t).transpose()*F*L*32.0/pow(M_PI,4.0);
    return F;
}
double Equation_u(Eigen::VectorXd X, double t){
    double aux = 0.0;
    for(int i = 0.0; i < 200.0; i ++){
        for(int j = 0.0; j < 200.0; j++){
            aux += pow(-1.0,(double)i+j)*cos((2.0*i+1.0)*M_PI*X(0)/L)*cos((2.0*j+1.0)*M_PI*X(1)/L)*
            (1.0/((2.0*i+1.0)*(2.0*j+1.0)*((2.0*i+1)*(2.0*i+1)+(2.0*j+1)*(2.0*j+1))));
        }
    }

    return aux*pow(L,2.0)*32.0/pow((double)M_PI,4.0);

}

double Equation_RBF(Eigen::VectorXd x , Eigen::VectorXd xj, double c2){
    double r2 = pow(x(0)-xj(0),2) + pow(x(1)-xj(1),2);
    return sqrt(r2 + c2);
}