#include"GMSolver.hpp"
#include"rectangle.hpp"
#define PROJECT_OUTSIDE
GMSolver::GMSolver(void){
    h = 0.001;
}
GMSolver::GMSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed){
        Configure(boundary_value_problem,boundary_parameters,discretization, seed);
}

void GMSolver::Update_Stat(double sol_0, double xi)
{

    sol = sol_0 + xi;
    sums[0] += sol;
    sums[1] += sol*sol;
    sums[2] += sol_0;
    sums[3] += sol_0*sol_0;
    sums[4] += xi;
    sums[5] += xi*xi;
    sums[6] += xi*sol_0;
    sums[7] += pow(sol - sol_a,2);
    N_trayectories ++;
    norm = 1.0/N_trayectories;
    err =fabs(mean-sol_a);
    rerr = err/sol_a;
    mean = sums[0] * norm;
    mse = sums[7] * norm;
    var = sums[1]*norm - mean*mean;
    std = sqrt(norm*var);
}

bool GMSolver::Inside(){
    bool stoppingbc = bvp.boundary.stop(E_P);
    double dist = bvp.boundary.Dist(params, X, E_P, N);

    if(stoppingbc){

      if( dist < -0.5826*(N.transpose()*sigma).norm()*sqrth){

          if (t > 0.0) {

              status = in;

          } else {

              status = time_out;

          }
      } else {

        if (t > 0.0) {

              status = stop;

          } else {

              X = E_P;
              status = time_out;

          }

      }

    }else{
      if( dist <= -0.5826*(N.transpose()*sigma).norm()*sqrth){

        if (t > 0.0) {

              status = in;

          } else {

              status = time_out;

          }
        } else {
            if (t > 0.0) {
          
              status = reflect;

            } else {

              X = E_P;
              status = time_out;

          }
      }
    }

    switch(status){

      case stop:
        ji_t = 0.0;
        X = E_P;
        return false;
        break;

      case reflect:
        ji_t = dist;
        X = E_P;
        return true;
        break;

      case time_out:
        ji_t = 0.0;
        t = 0.0;
        return false;
        break;

      default:
        ji_t = 0.0;
        return true;
        break;
    }
}
void GMSolver::LPG_Step(double rho, Boundary sten_boundary, Stencil stencil){
    if (d_k > -rho){
        do{
            Increment_Update();
            Xp = X + bvp.b.Value(X,t)*h + bvp.sigma.Value(X,t)*increment;
            N_rngcalls += X.size();
            omega =  gsl_ran_exponential(rng,2*h); //Exponential distribution with parameter 1/(2*h)
            uc = (N.transpose()*bvp.sigma.Value(X,t)).dot(increment) +N.transpose().dot(bvp.b.Value(X,t))*h;
            nu = 0.5 *(uc+sqrt(pow((N.transpose()*bvp.sigma.Value(X,t)).norm(),2.0)*omega+pow(uc,2.0)));
            //d_k = bvp.boundary.Dist(params, Xp,E_Pp,Np);
            if (d_k < -0.0) d_k = 0.0;
            ji_t = std::max(0.0,nu+d_k);
            Xp = Xp - ji_t*N;
        }while((Xp - E_P).dot(N)>0.0);
        d_k = sten_boundary.Dist(stencil.stencil_parameters, Xp,E_Pp,Np);
        if(d_k > 0.0){
            //printf("WARNING: The particle didn't enter in the domain  after Lepingle step.\n");
            Xp = E_Pp;
            ji_t += d_k;
            d_k = 0.0;
        }
    } else {
        Increment_Update();
        N_rngcalls += X.size();
        Xp = X + bvp.b.Value(X,t)*h + bvp.sigma.Value(X,t)*increment;
        ji_t = 0.0;
        d_k = bvp.boundary.Dist(params, Xp,E_Pp,Np);
        if(d_k > -0.0){
            Xp = E_Pp;
            ji_t = d_k;
            d_k = 0.0;
            //std::cout << ji_t << std::endl;
        }
    }
}
double GMSolver::Solve(Eigen::VectorXd X0, double T_start, double tolerance){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    N_trayectories = 0;
    N_rngcalls = 0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    sol_a = bvp.u.Value(X0,T_start);
    //Sample test
    for(unsigned int i = 0; i < 100; i++)
    {
        Reset(X0, T_start);
        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside());

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }
        Update_Stat(sol_0,xi);
    }
    //Trayectories are computed till the MSE we want to obtain is achieved
    do{
        Reset(X0, T_start);
        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside());
        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }
        Update_Stat(sol_0,xi);
    }while(2.0*std > tolerance*err);
    norm = 1.0/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*norm)*norm;
    var_xi = sums[5]*norm - sums[4]*sums[4]*norm*norm;
    pearson_c = covar/sqrt((sums[3]*norm - (sums[0]-sums[4])*
    (sums[0]-sums[4])*norm*norm)*var_xi);

    return mean;
}

double GMSolver::Solve(Eigen::VectorXd X0, double tolerance){
    Solve(X0,INFINITY,tolerance);
    return mean;
}
double GMSolver::Solve(Eigen::VectorXd X0, double T_start, unsigned int Ntray){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    sol_a = bvp.u.Value(X0,T_start);
    for(unsigned int i = 0; i < Ntray; i++)
    {
        Reset(X0, T_start);

        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
        }while(Inside());

        if(status == time_out){
          sol_0 = Y*bvp.p.Value(X,0.0f)+Z;
        } else {
          sol_0 = Y*bvp.g.Value(E_P,t)+Z;
        }

        Update_Stat(sol_0,xi);
    }
    norm = 1.0f/N_trayectories;
    covar =  (sums[6] - sums[2]*sums[4]*norm)*norm;
    var_xi = sums[5]*norm - sums[4]*sums[4]*norm*norm;
    pearson_c = covar/sqrt((sums[3]*norm - (sums[0]-sums[4])*
    (sums[0]-sums[4])*norm*norm)*var_xi);

    return mean;
}
double GMSolver::Solve(Eigen::VectorXd X0, unsigned int Ntray){
  return  Solve(X0, INFINITY, Ntray);
}

void GMSolver::Solve(Eigen::VectorXd X0, double c2, Stencil stencil, 
             std::vector<int> & G_j, std::vector<double> & G, double &B, 
             unsigned int Ntray, unsigned int N_job){
Solve(X0, INFINITY, c2, stencil, G_j, G, B, Ntray, N_job);
}
void GMSolver::Solve(Eigen::VectorXd X0, double T_start, double c2,
             Stencil & stencil, std::vector<int> & G_j,
             std::vector<double> & G, double & B,  unsigned int Ntray,
             unsigned int N_job){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    //Boundary of the stencil
    Boundary sten_boundary;
    sten_boundary._init_(Rectangle2D, Stopping);
    //G and B are emptied
    G_j.clear();
    G.clear();
    var_G.clear();
    double b = 0.0, bb = 0.0;
    B = 0.0;
    APL = 0.0;
    stencil.Reset();
    /*Control variable 
    -0 if trayectory still inside
    -1 if trayectory exits by stencil's boundary
    -2 if trayectory exits by problem's boundary*/
    int16_t control;
    Reset(X0,T_start);
    if(bvp.boundary.Dist(params, X0, E_P, N) >= -1E-06 && bvp.boundary.stop(X0)){
                B = bvp.g.Value(X0,t)/(1.0*Ntray/N_job);
                b = B;
                bb = B*B;
                //printf("OUT Node [ %f, %f ] Distance = %f N*sig = %f h = %f",X0[0],X0[1],bvp.boundary.Dist(stencil.stencil_parameters, X0, E_P, N),(N.transpose()*bvp.sigma.Value(X0,t)).norm(),h);
    } else{
      double h_cent;
      h_cent = h;
      if(bvp.boundary.stop(X0)){
        h = std::min(h,
        pow(bvp.boundary.Dist(stencil.stencil_parameters, X0, E_P, N)/(
        (N.transpose()*bvp.sigma.Value(X0,t)).norm()*2*0.5826),2.0));
      }
      sqrth = sqrt(h);
      //printf("Node [ %f, %f ] Distance = %f N*sig = %f h = %f\n",X0[0],X0[1],bvp.boundary.Dist(stencil.stencil_parameters, X0, E_P, N),(N.transpose()*bvp.sigma.Value(X0,t)).norm(),h);
      for(unsigned int i = 0; i < N_job; i++){
        Reset(X0,T_start);
        control = 0;
        sigma = bvp.sigma.Value(X0,t);
        do{
            Increment_Update();
            N_rngcalls += X0.size();
            VR_CV_Step();
            APL += h;
            if(t <= 0){
                control = 3;
            }
            if(bvp.boundary.Dist(params, X, E_P, N) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                control = 2;
                break;
            }
            if(sten_boundary.Dist(stencil.stencil_parameters, X, E_P, N) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                control = 1;
                break;
            }
        } while(control == 0);
        
        switch (control) {        
            case 1:
            #ifdef PROJECT_OUTSIDE
                stencil.Projection(X,E_P);
            #endif
                stencil.G_update(E_P,Y,bvp,c2);
                b += Z + xi;
                bb += pow(Z+xi,2.0);
            break;        
            case 2:
                b += Z + xi + Y*bvp.g.Value(E_P,t);
                bb += pow(Z + xi + Y*bvp.g.Value(E_P,t),2.0);
            break;
            case 3:
                b += Z + xi + Y*bvp.p.Value(X,t);
                bb += pow(Z + xi + Y*bvp.g.Value(E_P,t),2.0);
            break;
            default : 
                std::cout << "Something went wrong while solving";
        }
      }
      stencil.G_return_withrep(G_j, G, var_G,Ntray);
      B = b/N_job;
      var_B = bb/N_job-pow(B,2.0);
      B = b/Ntray;
      APL = APL/(h*N_job);
      h = h_cent;
      sqrth = sqrt(h_cent);
  }
}
void GMSolver::Solve_mix(Eigen::VectorXd X0, double c2,
             double rho, Stencil & stencil, std::vector<int> & G_j,
             std::vector<double> & G, double & B, unsigned int Ntray,
             unsigned int N_job){
                 Solve_mix( X0, INFINITY, c2, rho, stencil,  G_j, G, B, Ntray, N_job);
             }
void GMSolver::Solve_mix(Eigen::VectorXd X0, double T_start, double c2,
             double rho, Stencil & stencil, std::vector<int> & G_j,
             std::vector<double> & G, double & B, unsigned int Ntray,
             unsigned int N_job){
    for(unsigned int i = 0; i < 10; i++) sums[i] = 0.0;
    increment.resize(X0.size());
    X.resize(X0.size());
    E_P.resize(X0.size());
    N.resize(X0.size());
    //Boundary of the stencil
    Boundary sten_boundary;
    sten_boundary._init_(Rectangle2D, Stopping);
    //G and B are emptied
    G_j.clear();
    G.clear();
    var_G.clear();
    double b = 0.0, bb = 0.0;
    B = 0.0;
    APL = 0.0;
    stencil.Reset();
    /*Control variable 
    -0 if trayectory still inside
    -1 if trayectory exits by stencil's boundary
    -2 if trayectory exits by problem's boundary*/
    int16_t control;
    Reset(X0,T_start);
    if((bvp.boundary.Dist(params, X0, E_P, N) >= -1E-06) &&  bvp.boundary.stop(E_P)){
      B = bvp.g.Value(X0,t)/(1.0*Ntray/N_job);
      b = B;
      bb = B*B;
    } else{
      double h_cent = bvp.boundary.Dist(stencil.stencil_parameters, X0, E_P, N);
      h_cent = h;
      if(bvp.boundary.stop(E_P)){
        h = std::min(h,
        pow(bvp.boundary.Dist(stencil.stencil_parameters, X0, E_P, N)/(
        (N.transpose()*bvp.sigma.Value(X0,t)).norm()*2*0.5826),2.0));
        sqrth = sqrt(h);
      }
      Xp.resize(X0.size());
      E_Pp.resize(X0.size());
      Np.resize(X0.size());
      //printf("Node [ %f, %f ] Distance = %f N*sig = %f h = %f\n",X0[0],X0[1],bvp.boundary.Dist(stencil.stencil_parameters, X0, E_P, N),(N.transpose()*bvp.sigma.Value(X0,t)).norm(),h);
      for(unsigned int i = 0; i < N_job; i++){
        Reset(X0,T_start);
        d_k = sten_boundary.Dist(stencil.stencil_parameters, X, E_P, N);
        control = 0;
        do{
            if(bvp.boundary.stop(E_P)){
                sigma = bvp.sigma.Value(X,t);
                Increment_Update();
                N_rngcalls += X0.size();
                VR_CV_Step();
                APL += h;
                if(t <= 0){
                    control = 3;
                }
                if(bvp.boundary.Dist(params, X, E_P, Np) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                    control = 2;
                    break;
                }
                if(sten_boundary.Dist(stencil.stencil_parameters, X, E_P, N) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                    control = 1;
                    break;
                }
            } else {
                LPG_Step(rho,sten_boundary,stencil);
                N_rngcalls += X0.size();
                xi += Y*bvp.F.Value(X,t).dot(increment);
                Z += Y*(bvp.f.Value(X,t)*h + bvp.psi.Value(X,N,t)*ji_t);
                Y += bvp.c.Value(X,t)*Y*h + bvp.varphi.Value(X,N,t) * Y * ji_t;
                X = Xp;
                N = Np;
                E_P = E_Pp;
                t += - h;
                ji_t = 0.0;
                sigma = bvp.sigma.Value(X,t);
                /*if(bvp.boundary.stop(E_P)){
                    if(bvp.boundary.Dist(params, X, E_Pp, Np) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                    control = 2;
                    break;
                    }
                    if(sten_boundary.Dist(stencil.stencil_parameters, X, E_Pp, Np) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                    control = 1;
                    break;
                }
                }else{
                    if(sten_boundary.Dist(stencil.stencil_parameters, X, E_Pp, Np) > -0.5826*(N.transpose()*sigma).norm()*sqrth){
                    control = 1;
                    break;
                    }
                }*/
            }
        } while(control == 0);
        switch (control) {        
            case 1:
                #ifdef PROJECT_OUTSIDE
                stencil.Projection(X,E_P);
                #endif
                stencil.G_update(E_P,Y,bvp,c2);
                b += Z + xi;
                bb += pow(Z+xi,2.0);
            break;        
            case 2:
                b += Z + xi + Y*bvp.g.Value(E_P,t);
                bb += pow(Z + xi + Y*bvp.g.Value(E_P,t),2.0);
            break;
            case 3:
                b += Z + xi + Y*bvp.p.Value(X,t);
                bb += pow(Z + xi + Y*bvp.p.Value(X,0.0),2.0);
            break;
            default : 
                std::cout << "Something went wrong while solving";

        }
      }
      stencil.G_return_withrep(G_j, G, var_G,Ntray);
      B = b/N_job;
      var_B = bb/N_job-pow(B,2.0);
      var_B = var_B/Ntray;
      B = b/Ntray;
      APL = APL/(h*N_job);
      h = h_cent;
      sqrth = sqrt(h_cent);
  }
}
void GMSolver::Test(std::string filename, Eigen::VectorXd X0, double T_start, 
                  double tolerance, double h0,  unsigned int Nsamples){
    FILE *fp;
    fp = fopen(filename.c_str(),"w");
    time_t now = time(NULL);
    char *date = ctime(&now);
    int len = strlen(date);
    date[len-1] = ' ';
    printf("\n");
    printf("**********************************************************\n");
    printf("*** FKAKsolver class file version 0.1    produced by J. Moron ***\n");
    printf("*** C++ elliptic_test on %s         ***\n",date);
    printf("**********************************************************\n");
    printf("\n");
    printf("**********************************************************\n");
    printf("*** Convergence tests ************************************\n");
    printf("****** ****************************************************\n");
    printf("\n    h \t   exa_sol  num_sol   variance   pcoef \t    std       error   ");
    printf(" r.error    mse       cost \n---------------------------------");
    printf("-------------------------------------------------------------\n");
    fprintf(fp,"h,exact_sol,num_sol,variance,pcoef,std,err,rerr,mse,cost\n");
    fclose(fp);
    for(unsigned int i = 0; i < Nsamples; i++){
        h = h0 * pow(2,-1.0*i);
        sqrth = sqrt(h);
        Solve(X0, T_start, tolerance);
        fopen(filename.c_str(),"a");
        printf("%.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e  %.2e \n",
        h,sol_a, mean,var,pearson_c,std,err,rerr,mse,1.0*N_rngcalls);
        fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",
        h,sol_a, mean,var,pearson_c,std,err,rerr,mse,1.0*N_rngcalls);
        fclose(fp);
    }
}

void GMSolver::Test(std::string filename, Eigen::VectorXd X0, double tolerance, 
                  double h0, unsigned int Nsamples){
        Test(filename,X0,INFINITY,tolerance,h0,Nsamples);
    }




