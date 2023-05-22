#include "Moron_Poisson.hpp"
#include<iostream>
#include<stdio.h>
int main(){
    FILE *pf;
    Eigen::VectorXd aux;
    aux.resize(2);
    pf = fopen("Output/Debug/test_MM.csv","w");
    for(double x = -10.0; x < 10.0; x = x +0.5)
        for(double y = -10.0; y < 10.0; y = y +0.5){
            aux(0) = x; aux(1) = y;
            fprintf(pf,"%e,%e,%e,%e,%e,%e,%e,%e\n",x,y,Equation_u(aux,0.0),
            Equation_dudx(aux),Equation_dudy(aux),Equation_d2udx2(aux),
            Equation_d2udy2(aux), Equation_d2udxdy(aux));
        }
    fclose(pf);
    return 0;
}