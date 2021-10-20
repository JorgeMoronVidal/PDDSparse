#include "FKACSolver.hpp"
#include "stencil.hpp"
class GMSolver: public EMFKAC {
    private:
        /*-dist: Distance to the boundary*/
        double dist, norm;
        /*-Are bc's stopping?*/
        bool stoppingbc;
        bool Inside();
        /*Updates statistical quantities*/
        void Update_Stat(double sol_0, double xi);
        /*LPG variables*/
        double d_k, omega, uc, nu;
        Eigen::VectorXd E_Pp, Np, Xp;
        void LPG_Step(double rho, Boundary sten_boundary, Stencil stencil);
    public:
       /*Variance of he G matrix and the B vector*/
       std::vector<double> var_G;
       double var_B, APL;
        GMSolver(void);
        /*Class initialization
        -boundary_value_problem is a BVP object which stores all problem's equations
        -surface_parameters stores the parameters for the boundary construction
        -discretization stores the time discretization for the Stochastic process 
        -seed is the RNG seed*/  
        GMSolver(BVP boundary_value_problem,
            std::vector<double> boundary_parameters,
            double discretization,
            unsigned int seed);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator till 2*std < tolerance*err*/
        double Solve(Eigen::VectorXd X0, double T_start, double tolerance);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        double Solve(Eigen::VectorXd X0, double tolerance);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator with N trayectories*/
        double Solve(Eigen::VectorXd X0, double T_start, unsigned int Ntray);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        double Solve(Eigen::VectorXd X0, unsigned int Ntray);
        /*Solves an elliptic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator */
        void Solve(Eigen::VectorXd X0, double c2, Stencil stencil, 
                    std::vector<int> & G_j, std::vector<double> & G,
                    double & B, unsigned int Ntray, unsigned int N_job);
        /*Same but for mixed BC's*/
        void Solve_mix(Eigen::VectorXd X0,double c2,
             double rho, Stencil & stencil, std::vector<int> & G_j,
             std::vector<double> & G,  double & B, unsigned int Ntray, 
             unsigned int N_job);
        /*Solves a parabolic equation on inital point X0 using FKAK formula 
        approximated by Gobet-Menozzi's integrator with N trayectories*/
        void Solve(Eigen::VectorXd X0, double T_start, double c2,
                    Stencil & stencil, std::vector<int> & G_j, 
                    std::vector<double> & G, double &B,  unsigned int Ntray,
                    unsigned int N_job);
        /*Same but for mixed BC's Lepingle + GM integrator*/
        void Solve_mix(Eigen::VectorXd X0, double T_start, double c2,
             double rho, Stencil & stencil, std::vector<int> & G_j,
             std::vector<double> & G, double & B, unsigned int Ntray,
             unsigned int N_job);
        /*Convergence test for a given parabolic BVP*/
        void Test(std::string filename, Eigen::VectorXd X0, double T_start, 
                  double tolerance, double h0, unsigned int Nsamples);
        /*Convergence test for a given elliptic BVP*/
        void Test(std::string filename, Eigen::VectorXd X0, double tolerance, 
                  double h0, unsigned int Nsamples);
};