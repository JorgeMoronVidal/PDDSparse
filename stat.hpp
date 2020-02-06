#ifndef STAT
#define STAT

#include <eigen3/Eigen/Core>

/* Stat takes care of the variables and procedures involved
in the statistical characterization of the procedure and the
solution */
struct Stat {
    /*
    -sum is the sum of the solutions computed 
     for each trajectory [with control variates]
    -sq_summ is the sum of the square of the 
     solutions computed for each trajectory [with
     control variates]
    -var is the variance of the result
    -std is the standard deviation of the result
    -summ_0 is the sum of the solutions computed 
     for ecah trajectory [without control variates]
     -sq_summ is the sum of the square of the 
     solutions computed for each trajectory [without
     control variates]
    -summ_xi is the sum of the control variable
    -sqsumm_xi is the sum of the control variable's
     square value
    -varxi is the variance of the control variable.
    -crossum is the sum of xi*result
    */

    float summ, sq_summ, var, std,
    summ_0, sq_summ0, summ_xi, sqsumm_xi,
    var_xi, crosssumm, covar;

    /*
    -counter stores the total number of trayectories computed.
    -counterN stores the number of steps per path
    -summN stores the summ of counterN for various trayectories
    */
    int  counter, counterN, summN;
    Stat(void);
    
    /*Returns the average of the solutions for the trayectories computed
      and updates all the Stat quantities*/
    float update(Eigen::VectorXf X, float Y, float Z, float xi, int N);
    /*Fast version of update where some quantities of Stat are not upgraded*/
    float fast_update(Eigen::VectorXf X, float Y, float Z, float xi);

};  
#endif