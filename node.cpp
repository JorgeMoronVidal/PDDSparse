#include "node.hpp"
Node::Node(void){

    Y = 1;
    Z = 0;

}

void Node::init (Eigen::VectorXf X_init, 
                 float tol, 
                 float discretization){

    x0 = X_init;
    X.resize(x0.size());
    N.resize(x0.size());
    tolerance = tol;
    h = discretization;

}

void Node::init (Eigen::VectorXf X_init, 
                float discretization,
                int N_tray, 
                int node_index,
                std::vector<int> interface_index,
                std::vector<int> subdomain){

    x0 = X_init;
    X.resize(x0.size());
    N.resize(x0.size());
    E_P.resize(x0.size());
    N_trayectories = N_tray;
    h = discretization;
    sqh = sqrt(h);
    i_node = node_index;
    i_interface = interface_index;

}


void Node::Solve_FKAK(BVP bvp, 
                      float* params, 
                      gsl_rng *rng){

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
            for(int j = 0; 
                    j < increment.size(); 
                    j++) increment(j) = (float) gsl_ran_gaussian_ziggurat(rng,1)*sqh;
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

            for(int j = 0;
                j < increment.size(); 
                j++) increment(j) = gsl_ran_gaussian_ziggurat(rng,1)*sqh * sqh;

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

void Node::Solve_PDDSparse(BVP bvp,
                       gsl_rng *rng, 
                       float * parameters_stencil,
                       float * parameters_surface,
                       int N_tray,
                       float c2,
                       std::vector< std::vector<float> > & iPsi,
                       std::vector<Eigen::VectorXf> & stencil_position,
                       std::vector<int> & stencil_index,
                       std::vector<float> & G,
                       float & B){

    //Boundary of the stencil
    Boundary sten_boundary;
    sten_boundary._init_(Rectangle2D, Stopping);
    //Random increment for each step
    Eigen::VectorXf increment, mu;
    increment.resize(X.size());
    //H and B are emptied
    G.clear();
    G.resize(stencil_position.size());
    B = 0.0f;
    //Sigma Matrix
    Eigen::MatrixXf sigma;
    //Control variates variable
    float xi;
    //Accumulators for the B component
    float b_1 = 0.0f, b_2 = 0.0f;
    /*Control variable 
    -0 if trayectory still inside
    -1 if trayectory exits by stencil's boundary
    -2 if trayectory exits by problem's boundary*/
    int control;
    //Counters of trayectories depending on where they end
    int count_1 = 0, count_2 = 0;
    //We compute N trayectories 
    for(int i = 0; i < N_tray; i++){

        xi = 0.0f;
        X = x0;
        Y = 1.0f;
        Z = 0.0f;
        control = 0;

        do{
            sigma = bvp.sigma.Value(X,N);
            mu = bvp.mu.Value(X,N);
            for(int j = 0; 
                    j < increment.size(); 
                    j++) increment(j) = (float) gsl_ran_gaussian_ziggurat(rng,1)*sqh;
            
            xi+= Y*(sigma.transpose()*bvp.F.Value(X,N)).dot(increment);
            Z += bvp.f.Value(X,N)*Y*h;
            Y += bvp.c.Value(X, N)*Y*h + Y*bvp.mu.Value(X,N).transpose()*increment;
            X += (bvp.b.Value(X,N) - sigma*mu)*h + sigma*increment;

            if(bvp.boundary.Dist(parameters_surface, X, E_P, N) > 0.0f){
                control = 2;
                break;
            }else

            if(sten_boundary.Dist(parameters_stencil, X, E_P, N) > 0.0f){
                control = 1;
                break;
            }
        }while(control == 0);
        switch (control) {
        
             case 1: 
                for(unsigned int j; j < stencil_index.size(); j++){
                    G[j] += -iPsi[i_node][stencil_index[j]] * 
                    bvp.rbf.Value(E_P,x0,c2) * Y;
                }
                b_1 += Z;
                count_1 ++;
                break;
        
              case 2: 
                b_1 += Z;
                b_2 += bvp.g.Value(E_P,N);
                count_2 ++;
                break;
        
              default : 
              std::cout << "Something went wrong while solving";
        
        }

    }

    for(unsigned int j; j < stencil_index.size(); j++){
        G[j] = G[j]/count_1;
    }

    B = b_1/N_tray + b_2/count_2;

}

/*float Node::Random_Normal(std::mt19937 & random_int){
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

}*/