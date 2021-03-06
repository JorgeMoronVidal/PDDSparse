#include "lookuptable.hpp"
#include <iostream>
#include <stdlib.h> 
#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <string>

LookUpTable::LookUpTable(void){
}
void LookUpTable::Init(int dim, std::string file){

    //x and len are dinamically allocated
    x = new double*[dim];
    len = new int[dim];
    int  count = 0;
    std::ifstream infile;
    std::ostringstream  aux;
    std::string line;
    
    for(int i = 0; i < dim; i++){
        //We charge the grid coordinates
        aux << file + "/x_";
        aux << i;
        aux <<".txt";
        //We check if the file actually exist
        infile.open(aux.str().c_str(),std::ios::in);
        if(! infile){
            std::cout << aux.str() << " couldn't be opened.\n Make sure "<< 
            "x_"<<i<<".txt" << " is available in " <<
            file << std::flush;
            std::terminate();
        }
        //The file is ridden
        while (getline( infile, line )){
            count++;
        }
        len[i] = count;
        count = 0;
        x[i] = new double[len[i]];
        infile.clear();
        infile.seekg(0, std::ios::beg);
        while (getline( infile, line )){
            x[i][count] = atof(line.c_str());
            count++;
        }
        count = 0;
        infile.close();
    }
    
    //The file which stores the values of the function in the grid is open
    infile.open((file + "/value.txt").c_str());
    if(! infile){
        std::cout << file + "/value.txt" << " couldn't be opened.\n Make sure "<< 
        "value.txt" << " is available in " <<file << std::flush;
        std::terminate();
    }
    std::string cent;
    std::string::iterator it;
    
    int counter[dim];
    int mesh_size = 1;
    for (int i = 0; i < dim; i++){
        mesh_size *= len[i];
        counter[dim] = 0;
    }
    //It's values are stored in z
    z = new double[mesh_size];
    while (getline( infile, line )){
        
        it = line.begin();
        do{
            while(*it == ' '){

                it++;

            }
            while((*it != ' ') && (it != line.end())){

                cent += *it;
                it++;

            }
            z[counter[1] * len[0] + counter[0]] = atof(cent.c_str());
            counter[0]++;
            cent.clear();
            } while(it != line.end());

        counter[1]++;//row number
        counter[0] = 0;//column number
        
    }
    //Then the spline is initialized
    infile.close();
    T = gsl_interp2d_bicubic;
    spline = gsl_spline2d_alloc(T, len[0], len[1]);
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();
    gsl_spline2d_init(spline, x[0], x[1], z, len[0], len[1]);
}

//Evaluation of the interpolated function when only a position is asked
float LookUpTable::Eval(Eigen::VectorXf position){
    
    return gsl_spline2d_eval(spline, position(0), position(1), xacc, yacc);

}
//Evaluation of the interpolated function when a position and a normal vector have to be computed
float LookUpTable::Eval(float *params, 
                        Eigen::VectorXf & position, 
                        Eigen::VectorXf & exitpoint,
                        Eigen::VectorXf & normal)
{
    float dist = gsl_spline2d_eval(spline, position(0), position(1), xacc, yacc);

    normal(0) = gsl_spline2d_eval_deriv_x(spline, position(0), position(1), xacc, yacc);
    normal(1) = gsl_spline2d_eval_deriv_y(spline, position(0), position(1), xacc, yacc);

    exitpoint = position - normal* dist;

    if(dist > 0.0){

        position = exitpoint;

    }

    return dist;
}