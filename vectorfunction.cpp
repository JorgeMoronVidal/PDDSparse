#include"vectorfunction.hpp"

VectorFunction::VectorFunction(void){
    function = Default_Vector;
    analytic = true;
}

void VectorFunction::Init(pfvector input){
    function = input;
    analytic = true;
}

void VectorFunction::Init(int dim, 
                          std::string input){
    std::string aux;
    lookuptable.resize(dim);
    for(int i = 0; i < dim; i++){
        aux = input + "/" + std::to_string(i) +"_component";
        lookuptable[i].Init(dim, aux);
    }
    analytic = false;
}

Eigen::VectorXf VectorFunction::Value(Eigen::VectorXf position, 
                            Eigen::VectorXf normal){
    if(analytic){
        return function(position, normal);
    }

    Eigen::VectorXf out;
    out.resize(position.size());
    for(int i = 0; i < position.size(); i++ ){

        out(i) = lookuptable[i].Eval(position);
    }
    return out;
}

Eigen::VectorXf Default_Vector(Eigen::VectorXf position, Eigen::VectorXf normal){
    Eigen::VectorXf out;
    out.resize(position.size());
    return out;
}
