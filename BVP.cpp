#include"BVP.hpp"

BVP::BVP(void){

    //initialization of the control map
    control["f"] = false;
    control["c"] = false;
    control["u"] = false;
    control["varphi"] = false;
    control["F"] = false;
    control["mu"] = false;
    control["b"] = false;
    control["psi"] = false;
    control["sigma"] = false;
    control["domain"] = false;
    control["stop"] = false;
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

            f.Init(it->second);

        } else if (aux == "c"){

            c.Init(it->second);

        }else if (aux == "u"){

            u.Init(it->second);

        }else if (aux == "varphi"){

            varphi.Init(it->second);

        }
        control[aux] = true;
    } 

    
    for(std::map<std::string, pfvector>::iterator it = map_fvector.begin();
        it != map_fvector.end(); 
        it ++){

        aux = it->first;
        if(aux == "F"){

            F.Init(it->second);

        } else if (aux == "mu"){

            mu.Init(it->second);

        }else if (aux == "b"){

            b.Init(it->second);

        }else if (aux == "psi"){

            psi.Init(it->second);

        }

        control[aux] = true;
        
    }
 

}

void BVP::Surf_init(pfbound boundary, pfstop stopf)
{

}

void BVP::Surf_init(std::string boundary , pfstop stopf)
{

}