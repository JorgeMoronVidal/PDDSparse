#include"BVP.hpp"

BVP::BVP(void){

    //initialization of the control map
    control["f"] = false;
    control["c"] = false;
    control["u"] = false;
    control["g"] = false;
    control["varphi"] = false;
    control["F"] = false;
    control["mu"] = false;
    control["b"] = false;
    control["psi"] = false;
    control["sigma"] = false;
    control["domain"] = false;
    control["stop"] = false;
}

void BVP::BVP_init(int dim,
              std::map<std::string, pfscalar> map_fscalar,
              std::map<std::string, pfvector> map_fvector,
              std::map<std::string, pfmatrix> map_fmatrix,
              std::map<std::string, std::string> map_lut)
{       
    /* We read the maps and initialize the object components*/

    for(std::map<std::string, pfscalar>::iterator it = map_fscalar.begin();
        it != map_fscalar.end(); 
        it ++){

        if(it->first == "f"){

            f.Init(it->second);

        } else if (it->first == "c"){

            c.Init(it->second);

        }else if (it->first == "u"){

            u.Init(it->second);

        }else if (it->first == "varphi"){

            varphi.Init(it->second);

        }
        control[it->first] = true;
        std::cout << it->first << " was defined as an analytic function.\n";
    } 

    
    for(std::map<std::string, pfvector>::iterator it = map_fvector.begin();
        it != map_fvector.end(); 
        it ++){

        if(it->first == "F"){

            F.Init(it->second);

        } else if (it->first == "mu"){

            mu.Init(it->second);

        }else if (it->first == "b"){

            b.Init(it->second);

        }else if (it->first == "psi"){

            psi.Init(it->second);

        }

        control[it->first] = true;
        std::cout << it->first << " was defined as an analytic function. \n";


    }
     for(std::map<std::string, pfmatrix>::iterator it = map_fmatrix.begin();
        it != map_fmatrix.end(); 
        it ++){

        if(it->first == "sigma"){
            sigma.Init(it->second);
        }

        std::cout << it->first << " was defined as an analytic function. \n";
    }

    for(std::map<std::string,std::string>::iterator it = map_lut.begin();
        it != map_lut.end(); 
        it ++){

            if(!control[it->first]){

                if(it->first == "f"){

                    f.Init(dim, it->second);

                }else if (it->first == "c"){

                    c.Init(dim, it->second);

                }else if (it->first == "u"){

                    u.Init(dim, it->second);

                }else if (it->first == "varphi"){

                    varphi.Init(dim, it->second);

                }else if(it->first == "F"){

                    F.Init(dim, it->second);

                } else if (it->first == "mu"){

                    mu.Init(dim, it->second);

                }else if (it->first == "b"){

                    b.Init(dim, it->second);

                }else if (it->first == "psi"){

                    psi.Init(dim, it->second);

                }else if (it->first == "sigma"){

                    sigma.Init(dim, it->second);

                }
        }

        std::cout << it->first << " was defined as an interpolated look up table.\n";
    }

    for(std::map<std::string,bool>::iterator it = control.begin();
        it != control.end(); 
        it ++){

        if(!it->second){
            std::cout<< "Function " << it->first << "of the BVP wasn't defined.\n";
        }
    }

}

void BVP::Boundary_init(pfbound bound, pfstop stopf)
{
    boundary._init_(bound, stopf);
}

void BVP::Boundary_init(int dim, std::string bound , pfstop stopf)
{
    boundary._init_(dim, bound, stopf);
}