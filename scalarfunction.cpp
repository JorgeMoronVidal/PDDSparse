#include"scalarfunction.hpp"

ScalarFunction::ScalarFunction(void){
    function = Default_Scalar;
    analytic = true;
}

void ScalarFunction::Init(pfscalar input){
    function = input;
    analytic = true;
}

void ScalarFunction::Init(unsigned int dim, 
                          std::string input){
    lookuptable.Init(dim, input);
}

inline float ScalarFunction::Value(Eigen::VectorXf position, 
                            Eigen::VectorXf normal){
    if(analytic){
        return function(position, normal);
    }

    return lookuptable.Eval(position);
}

