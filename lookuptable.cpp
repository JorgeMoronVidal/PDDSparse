#include "lookuptable.hpp"
#include <iostream>
#include <fstream>
#include <string>

LookUpTable::LookUpTable(int dim, std::string label){

    x = new double*[dim];
    len = new int[dim];
    int  count = 0;
    std::ifstream infile;
    std::string aux,line;

    for(int i = 0; i < dim; i++){
        aux = label + "/x_" + std::to_string(i) + ".txt";
        infile.open(aux,std::ios::in);
        while (getline( infile, line )){
            count++;
        }
        len[i] = count;
        count = 0;
        x[i] = new double[len[i]];
        infile.clear();
        infile.seekg(0, std::ios::beg);
        while (getline( infile, line )){
            x[i][count] = stod(line);
            count++;
        }
        count = 0;
        infile.close();
    }
    infile.open(label + "/value.txt");

    std::string cent;
    std::string::iterator it;

    int counter[dim];
    int mesh_size = 1;
    for (int i = 0; i < dim; i++){
        mesh_size *= len[i];
        counter[dim] = 0;
    }
    z = new double[mesh_size];
    while (getline( infile, line )){
        it = line.begin();
        do{
            while(*it == ' '){

                it++;

            }
            for(it; *it != ' ' && it != line.end()
                ; it++){

                cent += *it;

            }
            z[counter[1] * len[0] + counter[0]] = stod(cent);
            cent.clear();
            counter[0]++;
            } while(it != line.end());

        counter[1]++;//row number
        counter[0] = 0;//column number
    }
    infile.close();
    T = gsl_interp2d_bicubic;
    spline = gsl_spline2d_alloc(T, len[0], len[1]);
    xacc = gsl_interp_accel_alloc();
    yacc = gsl_interp_accel_alloc();
    gsl_spline2d_init(spline, x[0], x[1], z, len[0], len[1]);
}

float LookUpTable::Eval(Eigen::VectorXf & position){
    float dist = gsl_spline2d_eval(spline, position(0), position(1), xacc, yacc);

    if(dist > 0.0){
        position(0) = position(0) - gsl_spline2d_eval_deriv_x(spline, position(0), position(1), xacc, yacc)*dist;
        position(1) = position(1) - gsl_spline2d_eval_deriv_y(spline, position(0), position(1), xacc, yacc)*dist;
    }
    
    return dist;
}

float LookUpTable::Eval(Eigen::VectorXf & position, Eigen::VectorXf & normal){
    float dist = gsl_spline2d_eval(spline, position(0), position(1), xacc, yacc);

    normal(0) = gsl_spline2d_eval_deriv_x(spline, position(0), position(1), xacc, yacc);
    normal(1) = gsl_spline2d_eval_deriv_y(spline, position(0), position(1), xacc, yacc);

    if(dist > 0.0){

        position = position - normal* dist;

    }

    return dist;
}