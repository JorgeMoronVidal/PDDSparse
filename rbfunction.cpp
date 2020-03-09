#include "rbfunction.hpp"

RBFunction::RBFunction(void){
    function = Default_RBF;
}

void RBFunction::Init(pRBF input){
    function = input;
}


float RBFunction::Value (Eigen::VectorXf x, 
                          Eigen::VectorXf x_i,
                          float c){
    return function(x, x_i, c);
}

float Default_RBF(Eigen::VectorXf x, 
                  Eigen::VectorXf x_i,
                  float c){
    return 0.0f;
}