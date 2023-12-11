//Header file for Sparse Linear System Class

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <fstream>
#include <Eigen/Dense>
#include <stdexcept>
#include "SparseMat.H"

class SLS {
    private:
        Eigen::VectorXd x, b, res, Ax;

    public:
        SparseMat& A; Mesh& mesh; SparseAddress& spa;
        SLS(SparseMat& SM);
        void set_diffusion_mat(double gamma, Eigen::MatrixXd& boundary_values, std::vector<int>& bc_type);
        void set_convection_mat(Eigen::Vector3d& u, Eigen::MatrixXd& boundary_values, std::vector<int>& bc_type, bool up_diff);
        void set_temporal_mat(double dt);
        void set_res();
        double calc_res_norm(int norm);
        double norm_factor();
        double Scaled_norm(int norm);
        double get_x_entry(int i) const;
        void set_x_entry(int i, double value);
        double get_b_entry(int i) const;
        void set_b_entry(int i, double value);
        void add_to_b_entry(int i, double value);
        void set_Ax();
        double calc_b_norm(int norm);
        void reset_x();
        void reset_b();
};
