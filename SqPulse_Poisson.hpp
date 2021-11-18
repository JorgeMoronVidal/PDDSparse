#include <eigen3/Eigen/Core>
#include <math.h>
#include <iostream>
#define N_KS 3
//A = height of the pulse 
//w = width of the pulse 
#define A_PROB  2.0
#define W_PROB  3.0
inline double A_k(int k){
    if(k == 0) return A_PROB*W_PROB;
    return 2.0*A_PROB*sin(0.5*k*W_PROB)/(1.0*k);
}
inline double Equation_u(Eigen::VectorXd X, double t){
    double out  = 0.0;
    for(int kx = -N_KS; kx < N_KS+1; kx++){
        for(int ky = -N_KS; ky < N_KS+1; ky++){
            out += A_k(kx)*A_k(ky)*cos(X(0)*kx)*cos(X(1)*ky);
        }
    }
    return out*pow(1.0/(2*N_KS+1),2);
}
inline double Equation_dudx(Eigen::VectorXd & X){
    double out  = 0.0;
    for(int kx = -N_KS; kx < N_KS+1; kx++){
        for(int ky = -N_KS; ky < N_KS+1; ky++){
            out += A_k(kx)*A_k(ky)*kx*sin(X(0)*kx)*cos(X(1)*ky);
        }
    }
    return out*(-1.0)*pow(1.0/(2*N_KS+1),2);
}
inline double Equation_dudy(Eigen::VectorXd & X){
    double out  = 0.0;
for(int kx = -N_KS; kx < N_KS+1; kx++){
        for(int ky = -N_KS; ky < N_KS+1; ky++){
            out += A_k(kx)*A_k(ky)*ky*cos(X(0)*kx)*sin(X(1)*ky);
        }
    }
    return out*(-1.0)*pow(1.0/(2*N_KS+1),2);
}
inline double Equation_d2udxdy(Eigen::VectorXd & X){
    double out  = 0.0;
for(int kx = -N_KS; kx < N_KS+1; kx++){
        for(int ky = -N_KS; ky < N_KS+1; ky++){
            out += A_k(kx)*A_k(ky)*ky*kx*sin(X(0)*kx)*sin(X(1)*ky);
        }
    }
    return out*1.0*pow(1.0/(2*N_KS+1),2);
}
inline double Equation_d2udx2(Eigen::VectorXd & X){
    double out  = 0.0;
for(int kx = -N_KS; kx < N_KS+1; kx++){
        for(int ky = -N_KS; ky < N_KS+1; ky++){
            out += A_k(kx)*A_k(ky)*kx*kx*cos(X(0)*kx)*cos(X(1)*ky);
        }
    }
    return out*(-1.0)*pow(1.0/(2*N_KS+1),2);
}
inline double Equation_d2udy2(Eigen::VectorXd & X){
    double out  = 0.0;
for(int kx = -N_KS; kx < N_KS+1; kx++){
        for(int ky = -N_KS; ky < N_KS+1; ky++){
            out += A_k(kx)*A_k(ky)*ky*ky*cos(X(0)*kx)*cos(X(1)*ky);
        }
    }
    return out*(-1.0)*pow(1.0/(2*N_KS+1),2);
}
inline double Equation_c(Eigen::VectorXd X, double t){
    return -0.0;
}
inline double Equation_g(Eigen::VectorXd X, double t){
    return Equation_u(X,t);
}
inline double Equation_p(Eigen::VectorXd X, double t){
    return Equation_u(X,t);
}
inline Eigen::VectorXd Equation_b(Eigen::VectorXd X, double t){
    Eigen::VectorXd Output(2);
    return Output* 0.0;
}
inline Eigen::MatrixXd Equation_sigma(Eigen::VectorXd X, double t){
    return Eigen::MatrixXd::Identity(2,2) * 1.41421356237;
}
inline double Equation_f(Eigen::VectorXd X, double t){
    return -Equation_d2udx2(X)-Equation_d2udy2(X);
}
inline double Equation_Varphi(Eigen::VectorXd X,Eigen::VectorXd normal, double t){
    return 0.0;
}
inline double Equation_Psi(Eigen::VectorXd X, Eigen::VectorXd normal, double t){
    return 0.0;
}
inline Eigen::VectorXd Equation_F(Eigen::VectorXd X, double t){
    Eigen::VectorXd F(2);
    F(0) = Equation_dudx(X);
    F(1) = Equation_dudy(X);
    F = -Equation_sigma(X,t).transpose()*F;
    return F;
}
inline bool Stopping_mix(Eigen::VectorXd X){
    //if(fabs(X(0) + 1.0) < 1E-8) return false;
    return true;
}
inline double Equation_RBF(Eigen::VectorXd x , Eigen::VectorXd xj, double c2){
    return 1/(pow((x-xj).norm(),2) + c2);
}
