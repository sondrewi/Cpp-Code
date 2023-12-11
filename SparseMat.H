//Headr file for SparseMat class

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <fstream>
#include <Eigen/Dense>
#include <stdexcept>
#include "SparseAddress.H"

class SparseMat {
    private:
        std::vector<double> entries;

    public:
        SparseAddress& spa;
        Mesh& mesh;
        SparseMat(SparseAddress& sa);
        void add_to_entry(int i, int j, double entry);
        bool symmetric() const;
        bool diagonal_dominance() const;
        bool spd() const;
        double compute_det(int n, int start_row, int start_col, std::vector<int> skip) const;
        double get_entry2D(int i, int j) const;
        double get_entry1D(int i) const;
        Eigen::VectorXd operator*(const Eigen::VectorXd& v) const;
        void reset();
};
