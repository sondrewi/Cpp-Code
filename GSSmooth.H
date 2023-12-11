//Header file for GSSmooth class that will carry out Gauss-Seidel
//sweeps on a sparse linear system

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <fstream>
#include <Eigen/Dense>
#include <stdexcept>
#include "SLS.H"

class GSSmooth {
    public:
        SparseMat& A; Mesh& mesh; SparseAddress& spa; SLS& sls;
        GSSmooth(SLS& SparseSystem);

        void smooth(int nSweeps);
};
