#include<iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <eigen3/Eigen/Core>
#include <fstream>
#include <unordered_map>

int ReadFromFile(std::string file,
                 float &h,
                 int &N_tray,
                 int &N_subdo,
                 int &N_node, 
                 Eigen::VectorXf & SW,
                 Eigen::VectorXf & NE)
  {
    std::ifstream myfile(file);
    std::string line, cent;
  
  if (myfile)  // same as: if (myfile.good())
    {
        while (getline( myfile, line )){  // same as: while (getline( myfile, line ).good())

            if (line.front()!='%') //First of all we check if the line is a comment
            { 
                cent.clear(); //Centinel needs to be empty

                for(std::string::iterator it = line.begin(); it != line.end(); it++){
                //Loop over all the characters of the line
                    
                    if(*it != ' ' && *it != '\t'){
                    //If it is not a space nor a tabulation then we concatenate the term to cent
                        cent += *it;
                    }
                    else{
                    //If we get a space or a tabulation then we check the string
                        if(cent == "h="){
                        //TypeOfPDE
                            while( *it == ' ' or *it == '\t'){

                                        it++;

                            }

                            cent.clear();

                            for(it; it != line.end(); it++){

                                cent += *it;

                            }
                            
                            h = atof(cent.c_str());
                            break;
                        }

                        if(cent == "SW="){
                        //South West point of the domain 
                            SW.resize(2);
                            while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            for(int i = 0; i < 2; i++){
                                for(it; *it != ' ' && *it != '\t' && it != line.end()
                                ; it++){

                                    cent += *it;

                                }
                                SW(i) = atof(cent.c_str());
                                cent.clear();
                                while( *it == ' ' or *it == '\t'){
                                        it++;
                                }
                    
                            }

                            break;
                        }

                        if(cent == "NE="){
                        //South West point of the domain 
                            NE.resize(2);
                            while( *it == ' ' or *it == '\t'){
                                        it++;

                            }
                            cent.clear();
                            for(int i = 0; i < 2; i++){
                                for(it; *it != ' ' && *it != '\t' && it != line.end()
                                ; it++){

                                    cent += *it;

                                }
                                NE(i) = atof(cent.c_str());
                                cent.clear();
                                while( *it == ' ' or *it == '\t'){
                                        it++;
                                }
                    
                            }

                            break;
                        }

                        if(cent == "N_tray="){
                        //Number of trayectories per node

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            for(it; it != line.end(); it++){

                                cent += *it;

                            }
                            
                            N_tray = atoi(cent.c_str());
                            break;
                        }

                        if(cent == "N_subdo="){
                        //Number of domains in each direction

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            for(it; it != line.end(); it++){

                                cent += *it;

                            }
                            
                            N_subdo = atoi(cent.c_str());
                            break;
                        }

                        if(cent == "N_node="){
                        //Number nodes per interface

                            while( *it == ' ' or *it == '\t'){

                                        it++;
                            }

                            cent.clear();

                            for(it; it != line.end(); it++){

                                cent += *it;

                            }
                            
                            N_node = atoi(cent.c_str());
                            break;
                        }
                        
                    }
                }
            }
        }
        myfile.close();
        return 1;
    }

    else {
    
        std::cout << "WARNING: Cofiguration file was not found\n";
        std::cout << "Make sure it exist a configuration.txt file in the program's directory.\n";
        return 0;
    }    
 }

 int main(void){
    std::string file = "configuration.txt";
    float h;
    int N_tray,N_subdo, N_node; 
    Eigen::VectorXf SW,NE;
    ReadFromFile(file, h,N_tray,N_subdo,N_node, SW, NE);

    printf("*************PDD SPARSE*****************\n");
    printf("***************Ver 0.2******************\n");
    printf("***********Jorge Moron Vidal************\n");
    printf("**********jmoron@math.uc3m.es***********\n");
    return 1;
 }