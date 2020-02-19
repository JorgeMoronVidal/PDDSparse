#include <iostream>
#include <fstream>
#include <iterator>
#include "BVP.hpp"
#include "equation.hpp"
#include "rectangle.hpp"
#define N 101

int main() {
    BVP bvp;

    //We initialize the variables to charge
    std::map<std::string, pfscalar> scalar_init;
    std::map<std::string, pfvector> vector_init;
    std::map<std::string, pfmatrix> matrix_init;
    std::map<std::string, std::string> string_init;

    scalar_init["f"] = Equation_f;
    scalar_init["c"] = Equation_c;
    scalar_init["u"] = Equation_u;
    scalar_init["g"] = Equation_g;
    vector_init["b"] = Equation_b;
    vector_init["F"] = Equation_F;
    matrix_init["sigma"] = Equation_sigma;
    bvp.Boundary_init(Rectangle2D, Stopping);
    bvp.BVP_init(2,scalar_init, vector_init, matrix_init,string_init);
    
    float params[4] = {-1.0f, -1.0f, 1.0f, 1.0f};
    Eigen::VectorXf Pos, Norm, E_P;
    std::string chain;
    Pos.resize(2);
    Norm = Pos;
    E_P = Pos;
    
    float X[N] , Y[N];
    std::ofstream file_x, file_y, file_val;

    for(int i = 0; i < N; i++){

        X[i] = -3.0f + 0.06f * i;
        Y[i] = -3.0f + 0.06f * i;

    }


    for(std::map<std::string, pfscalar>::iterator it = scalar_init.begin();
        it != scalar_init.end(); 
        it ++){
            
            chain = "lookuptable/"+it->first;
            std::cout << chain << "\n";
            string_init[it->first] = chain;
            file_x.open(chain + "/x_0.txt");
            file_x.close();
            file_y.open(chain+ "/x_1.txt");
            file_y.close();
            file_val.open(chain + "/value.txt");
            file_val.close();
            file_x.open("lookuptable/"+it->first + "/x_0.txt",std::ios::app);
            file_y.open("lookuptable/"+it->first + "/x_1.txt",std::ios::app);
            file_val.open("lookuptable/"+it->first + "/value.txt",std::ios::app);

            for(int i = 0; i < N; i++){
                file_x << X[i] << "\n";
                file_y << Y[i] << "\n";
                Pos(0) = Y[i];
                for(int j = 0; j < N; j++){
                    
                    Pos(1) = X[i];
                    file_val << it->second(Pos,Norm) << " ";
                }
                file_val << "\n";
            }
            file_x.close();
            file_y.close();
            file_val.close();

        }
        
        Eigen::VectorXf test;

        test.resize(2);
    
        for(std::map<std::string, pfvector>::iterator it = vector_init.begin();
        it != vector_init.end(); 
        it ++){
            
            for(int i = 0; i < 2; i++){

                chain = "lookuptable/"+it->first + "/" + std::to_string(i) +"_component";
                std::cout << chain << "\n";
                string_init[it->first] = chain;
                file_x.open(chain + "/x_0.txt");
                file_x.close();
                file_y.open(chain+ "/x_1.txt");
                file_y.close();
                file_val.open(chain + "/value.txt");
                file_val.close();
                file_x.open(chain + "/x_0.txt",std::ios::app);
                file_y.open(chain + "/x_1.txt",std::ios::app);
                file_val.open(chain + "/value.txt",std::ios::app);

                for(int j = 0; j < N; j++){

                    file_x << X[j] << "\n";
                    file_y << Y[j] << "\n";
                    Pos(0) = Y[j];

                    for(int l = 0; l < N; l++){
                    
                        Pos(1) = X[l];
                        file_val << it->second(Pos,Norm)(i) << " ";
                    }
                    file_val << "\n";
                }
                file_x.close();
                file_y.close();
                file_val.close();
            }
            
        }


}