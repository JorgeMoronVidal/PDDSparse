#ifndef NODE
#define NODE

#include <iostream>
#include <map>
#include <eigen3/Eigen/Core>
#include "lookuptable.hpp"
#include "BVP.hpp"

class Node{
    private:
        Eigen::VectorXf x0;
        void UpdateStat(float &sum , float &sqsum, int &n, int sumsteps, int steps);
    public:

    float h, solution, tolerance, std, var, pearcoef, apl;
    Node(void);
    void init (Eigen::VectorXf X_init, float tol, float discrezitation);
    void Solve(BVP bvp);

};
#endif