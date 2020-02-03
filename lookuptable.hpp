#ifndef BORDERHEADERDEF
#define BORDERHEADERDEF
#include <eigen3/Eigen/Core>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

class LookUpTable{

    private:

    double** x;
    double* z;
    int *len;

    const gsl_interp2d_type *T;
    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    
    public:
    
    LookUpTable(int dim, std::string label);
    
    float Eval(Eigen::VectorXf & position);
    float Eval(Eigen::VectorXf & position, Eigen::VectorXf & normal);
};
#endif