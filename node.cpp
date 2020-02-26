#include "node.hpp"
Node::Node(void){

    Y = 1;
    Z = 0;
    

}
void Node::init (Eigen::VectorXf X_init, 
                 float tol, 
                 float discretization, 
                 unsigned int random_seed){

    x0 = X_init;
    X.resize(x0.size());
    N.resize(x0.size());
    tolerance = tol;
    h = discretization;
    seed = random_seed;

}


void Node::Solve_FKAK(BVP bvp){

    /*
    -sum is the sum of the solutions computed 
    for each trajectory [with control variates]
    -sq_summ is the sum of the square of the 
    solutions computed for each trajectory [with
    control variates]
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

    float sol = 0.0f, summ = 0.0f, sq_summ = 0.0f, summ_0 = 0.0f, sq_summ0 = 0.0f,
          xi = 0.0f, summ_xi = 0.0f, sqsumm_xi = 0.0f, var_xi = 0.0f, crossumm = 0.0f;
        
    /*
    -counter stores the total number of trayectories computed.
    -counterN stores the number of steps per path
    -summN stores the summ of counterN for various trayectories
    */

    int  counter, counterN, summN;

    //random number generator
    thread_local std::mt19937 i_random(seed);
    Eigen::VectorXf increment;
    increment.resize(X.size());

    for(int i = 0; i < 1000; i++){
            for(int j = 0; j < increment.size(); j++) increment(j) = Random_Normal(i_random());
            xi+= Y*(bvp.sigma(X,t).transpose()*F(X,t)).dot(bmotion.increment);
            Z += f(X,t)*Y*bmotion.dt;
            Y += c(X, t)*Y*bmotion.dt + Y*mu(X,t).transpose()*bmotion.increment;
            X += (b(X,t) - sigma(X,t)*mu(X,t))*bmotion.dt + sigma(X,t)*bmotion.increment;
            t += bmotion.dt;
            counterN++;
    }
    
}

float Node::Random_Normal(int random_int){
    if(generate) {

        generate = false;

        return spare;

    } else {

        do {

            u = ((float)random_int/max_mt) * 2.0 - 1.0;

            v = ((float)random_int/max_mt) * 2.0 - 1.0;

            s = u * u + v * v;

        } while (s >= 1.0 || s == 0.0);

        s = sqrt(-2.0 * log(s) / s);

        spare = v * s;

        generate = true;

        return (float)(u * s);

    }

}

void Update_Stat(float sol_0, float xi, float & summ, float & mean,
                 float & sqsumm, float & var, float & std, float & summ_0, 
                 float & sqsumm_0, float & summ_xi, float & sqsumm_xi, 
                 float & var_xi, float & crossumm, int & counter, int & counterN, 
                 int & summN){

    float sol = sol + xi;
    summ += sol;
    sqsumm += sol*sol;
    summ_0 += sol_0;
    sqsumm_0 += sol_0*sol_0;
    summ_xi += xi;
    sqsumm_xi += xi*xi;
    crossumm = xi*sol_0;
    summN += counterN;
    counterN = 0;
    counter ++;
    float aux = 1.0f/counter;
    mean = summ * aux;
    var = sqsumm*aux - mean*mean;
    std = sqrt(aux*var);

}

