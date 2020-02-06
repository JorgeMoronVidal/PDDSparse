#include"BVP.hpp"

float Default_Scalar(Eigen::VectorXf X)
{   
    return 0.0f;
}

Eigen::VectorXf Default_Vector(Eigen::VectorXf X)
{   
    Eigen::VectorXf out;
    out.resize(X.size());
    return out;
}

Eigen::MatrixXf Default_Matrix(Eigen::VectorXf X)
{
    Eigen::MatrixXf out;
    out.resize(X.size());
    return out;
}

BVP::BVP(void){
    f =  Default_Scalar;
    c = Default_Scalar;
    u = Default_Scalar;
    varphi = Default_Scalar;

    F = Default_Vector;
    mu = Default_Vector;
    b = Default_Vector;
    psi = Default_Vector;

    sigma = Default_Matrix;
}
void BVP::BVP_init(std::map<std::string, pfscalar> map_fscalar,
              std::map<std::string, pfvector> map_fvector,
              std::map<std::string, pfmatrix> sigma,
              std::map<std::string, std::string> map_lut)
{       
    /* We read the maps and initialize the object components*/
    std::string aux;

    for(std::map<std::string, pfscalar>::iterator it = map_fscalar.begin();
        it != map_fscalar.end(); 
        it ++){

        aux = it->first;
        if(aux == "f"){
            f = it->second;
        } else if (aux == "c"){
            c = it->second;
        }else if (aux == "u"){
            u = it->second;
        }else if (aux == "varphi"){
            varphi = it->second;
        }
    } 

    
    for(std::map<std::string, pfvector>::iterator it = map_fvector.begin();
        it != map_fvector.end(); 
        it ++){

        aux = it->first;
        if(aux == "F"){
            F = it->second;
        } else if (aux == "mu"){
            mu = it->second;
        }else if (aux == "b"){
            b = it->second;
        }else if (aux == "psi"){
            psi = it->second;
        }
    }
 

}

void BVP::Surf_init(pfbound boundary, pfstop stopf)
{

}

void BVP::Surf_init(std::string boundary , pfstop stopf)
{

}