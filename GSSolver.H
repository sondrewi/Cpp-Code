//Header file for Gauss-Seidel Solver Class

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <fstream>
#include <Eigen/Dense>
#include <stdexcept>
#include "GSSmooth.H"

class GSSolver {
    private:
        int max_it;
        int sweeps;

    public:
        SparseMat& A; Mesh& mesh; SparseAddress& spa; SLS& sls; GSSmooth& smoother;
        GSSolver(GSSmooth& smoother_object, int sweeps_, int maxit);
        void solve(double relTol);
};
