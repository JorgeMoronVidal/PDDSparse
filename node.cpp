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
    generate = false;
    X.resize(x0.size());
    N.resize(x0.size());
    tolerance = tol;
    h = discretization;
    seed = random_seed;

}

void Node::init (Eigen::VectorXf X_init, 
                float tol, 
                float discretization, 
                unsigned int random_seed,
                int node_index,
                std::vector<int> interface_index,
                std::vector<int> subdomains,
                Eigen::SparseMatrix<float> H_mtrx,
                Eigen::SparseMatrix<float> psi_m){

    x0 = X_init;
    generate = false;
    X.resize(x0.size());
    N.resize(x0.size());
    tolerance = tol;
    h = discretization;
    seed = random_seed;
    H = H_mtrx;
    i_node = node_index;
    i_interface = interface_index;

}


void Node::Solve_FKAK(BVP bvp, float* params){

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

    float mean = 0.0f, sol_0 = 0.0f, summ = 0.0f, sqsumm = 0.0f, summ_0 = 0.0f, sqsumm_0 = 0.0f,
          xi = 0.0f, summ_xi = 0.0f, sqsumm_xi = 0.0f, var_xi = 0.0f, crossumm = 0.0f,
          sqh = sqrt(h);
        
    /*
    -counter stores the total number of trayectories computed.
    -counterN stores the number of steps per path
    -summN stores the summ of counterN for various trayectories
    */

    int  counter = 0, counterN = 0, summN = 0;

    //random number generator
    thread_local std::mt19937 i_random(seed);
    Eigen::VectorXf increment;
    increment.resize(X.size());

    Eigen::MatrixXf sigma;
    Eigen::VectorXf mu;

    for(int i = 0; i < 1000; i++){

        xi = 0.0f;
        X = x0;
        Y = 1.0f;
        Z = 0.0f;

        do{
            sigma = bvp.sigma.Value(X,N);
            mu = bvp.mu.Value(X,N);
            for(int j = 0; j < increment.size(); j++) increment(j) = Random_Normal(i_random) * sqh;
            xi+= Y*(sigma.transpose()*bvp.F.Value(X,N)).dot(increment);
            Z += bvp.f.Value(X,N)*Y*h;
            Y += bvp.c.Value(X, N)*Y*h + Y*bvp.mu.Value(X,N).transpose()*increment;
            X += (bvp.b.Value(X,N) - sigma*mu)*h + sigma*increment;
            counterN++;
    }while(bvp.boundary.Dist(params, X, E_P,N) < 0.0f);
    X = E_P;
    sol_0 = Y*bvp.g.Value(X,N)+Z;
    Update_Stat(sol_0,xi, summ, mean, sqsumm, summ_0, sqsumm_0, summ_xi,
    sqsumm_xi, var_xi, crossumm, counter, counterN, summN);
    }

    do{

        xi = 0.0f;
        X = x0;
        Y = 1.0f;
        Z = 0.0f;

        do{
            sigma = bvp.sigma.Value(X,N);
            mu = bvp.mu.Value(X,N);

            for(int j = 0; j < increment.size(); j++) increment(j) = Random_Normal(i_random) * sqh;

            xi+= Y*(sigma.transpose()*bvp.F.Value(X,N)).dot(increment);
            Z += bvp.f.Value(X,N)*Y*h;
            Y += bvp.c.Value(X, N)*Y*h + Y*bvp.mu.Value(X,N).transpose()*increment;
            X += (bvp.b.Value(X,N) - sigma*mu)*h + sigma*increment;
            counterN++;
        
    }while(bvp.boundary.Dist(params, X, E_P,N) < 0.0f);
    X = E_P;
    sol_0 = Y*bvp.g.Value(X,N)+Z;
    Update_Stat(sol_0,xi, summ, mean, sqsumm, summ_0, sqsumm_0, summ_xi,
    sqsumm_xi, var_xi, crossumm, counter, counterN, summN);

    //std::cout << std*2.0f/mean << "\n";
    }while(std*2.0f/mean > tolerance*fabs(mean-bvp.u.Value(X,N)));

    float aux = 1.0f/counter;
    covar =  (crossumm - summ_0*summ_xi*aux)*aux;
    var_xi = sqsumm_xi*aux - summ_xi*summ_xi*aux*aux;
    pearson_c = covar/sqrt((sqsumm_0*aux - (summ-summ_xi)*
    (summ-summ_xi)*aux*aux)*var_xi);
    solution = mean;

}
float Node::Random_Normal(std::mt19937 & random_int){
    if(generate) {

        generate = false;

        return spare;

    } else {

        do {

            u = ((float)random_int()/max_mt) * 2.0f - 1.0f;

            v = ((float)random_int()/max_mt) * 2.0f - 1.0f;

            s = u * u + v * v;

        } while (s >= 1.0 || s == 0.0);

        s = sqrt(-2.0 * log(s) / s);

        spare = v * s;

        generate = true;

        return (u * s);

    }

}

void Node::Update_Stat(float sol_0, float xi, float & summ, float & mean,
                 float & sqsumm,  float & summ_0, 
                 float & sqsumm_0, float & summ_xi, float & sqsumm_xi, 
                 float & var_xi, float & crossumm, int & counter, int & counterN, 
                 int & summN){

    float sol = sol_0 + xi;
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

