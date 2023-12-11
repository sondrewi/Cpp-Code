//Header file for Sparse Address class

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <fstream>
#include <Eigen/Dense>
#include <stdexcept>
#include <cstdlib>
#include "Practical1/Mesh.H"

class SparseAddress {
    private:
        int nr_equations;
        int matrix_band;

        std::vector<int> col_index;
        std::vector<int> row_index;

    public:
        Mesh& mesh;
        SparseAddress(Mesh& the_mesh);
        int get_nr_eqs() const;
        int nr_off_diag() const;
        int mat_band() const;
        int get_col_index(int) const;
        int get_row_index(int) const;
        std::pair<int,int> oneD_to_twoD(int i) const;
        int twoD_to_oneD(int i, int j, bool trunc = false) const;
};
