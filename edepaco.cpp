#include "edepaco.hpp"

float cpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return -pow(X.norm(),2)* exp(-X(0));
}
float fpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return 2.0f * sin(7.0f*X(0) + 9.0f*X(1)) + 130.0f*
    cos(7.0f*X(0) + 9.0f*X(1)) -cpaco(X, N) * upaco(X, N);
}
float ppaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return 2.1f + cos(7.0f * X(0) + 9.0f * X(1));
}
float gpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return 2.1f + cos(7.0f * X(0) + 9.0f * X(1));
}
float psipaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return -(N(0) * 7.0f + N(1) * 9.0f) * sin(7.0f * X(0) + 9.0f * X(1)) 
            + upaco(X,N);
}
float varpsipaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return -1.0f;
}

bool Mixpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    if( X(1) < 0.0f){
        return false;
    }
    return true;
}

bool Dirichpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return true;
}

bool Robinpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return false;
}

Eigen::VectorXf bpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    Eigen::VectorXf b(X.size());
    b << 1.0f, 0.0f;
    return b;
}
Eigen::VectorXf Fpaco(Eigen::VectorXf X, Eigen::VectorXf N){
    Eigen::VectorXf F(X.size());
    F << 7.0f, 9.0f;
    F = F*sin(7.0f * X(0) + 9.0f * X(1) );
    return F;
}
Eigen::MatrixXf sigmapaco(Eigen::VectorXf X, Eigen::VectorXf N){
    Eigen::MatrixXf sigma(X.size(), X.size());
    sigma << 1.41421356237f, 0.0f,
             0.0f, 1.41421356237f;
    return sigma;
}
float upaco(Eigen::VectorXf X, Eigen::VectorXf N){
    return 2.1f + cos(7.0f * X(0) + 9.0f * X(1));
}