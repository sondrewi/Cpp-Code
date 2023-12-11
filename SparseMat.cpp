#include "SparseMat.H"

/*Class to store coefficients of sparse matrix. Constructed using
a SparseAddress class object*/

SparseMat::SparseMat(SparseAddress& sa) : spa(sa), mesh(spa.mesh){

    //Constructor that initialises entries vector
    entries.resize(spa.get_nr_eqs() + spa.nr_off_diag());
    std::fill(entries.begin(), entries.end(), 0);

}

//Function to add to an entry in sparse matrix
void SparseMat::add_to_entry(int i, int j, double entry){
    //Set a non-zero entry (error will be raised if coefficients for zero-entry supplied)
    int index = spa.twoD_to_oneD(i, j);
    entries[index] += entry;
}

//Get entry of original matrix. We allow zero to be returned if
//indices correspond to zero entry of original matrix. Hence function twoD_to_OneD
//from SparseAddress class not utilised
double SparseMat::get_entry2D(int i, int j) const {
    if ((i > spa.get_nr_eqs()) || (j > spa.get_nr_eqs()))  {
        std::cerr << "Index beyond scope of matrix" << std::endl;
        exit(EXIT_FAILURE);
    }

    int index = spa.get_row_index(i);
    bool found_index = false;

    //Loop over indices of col_index known to correspond to row i
    for (int k = index; k < spa.get_row_index(i+1); k++){
        if (spa.get_col_index(k) == j){
            found_index = true;
            index = k;
        }
    }

    if (!found_index) {
        return 0.0;
    }

    else{return entries[index];}
}

double SparseMat::get_entry1D(int i) const {
    return entries[i];
}

//Check if matrix is symmetric
bool SparseMat::symmetric() const{
    bool symmetric = true;

    /*Go along col_index. If entry lies above diagonal, check if it is
    equal to its symmetric counterpart*/
    for (int row_col=1; row_col < spa.get_nr_eqs(); row_col++){

        //Loop over entries below diagonal
        for (int j = 0; j < row_col; j++){
            if(get_entry2D(row_col, j) != get_entry2D(j, row_col)){
                symmetric = false;
                break;
            }
        }
    }
    return symmetric;
}

/*Check if Matrix is diagonally dominant
Diagonal dominance defined as all sums absolute off-diagonals
being less than or equal to their respective diagonal entry
with at least one row having a diagonal entry larger than the sum
of its absolute off-diagonals*/
bool SparseMat::diagonal_dominance() const{
    bool diag_dom = false;

    double diagonal_entry;
    double off_diag_sum = 0;

    //Loop over rows
    for (int i=0; i < spa.get_nr_eqs(); i++){

        //Loop over non-zero entries in row
        for (int j = spa.get_row_index(i); j < spa.get_row_index(i+1); j++){
            if(spa.get_col_index(j) == i){diagonal_entry = entries[j];}
            else{off_diag_sum += abs(entries[j]);}
        }

        //Check for diagonal dominance in given row
        if (off_diag_sum > diagonal_entry){
            diag_dom = false;
            break;
        }
    }

    return diag_dom;
}

//Check for symmetric positive definiteness of matrix
bool SparseMat::spd() const{
    bool sym_pd;
    bool sym = symmetric();
    bool dd = diagonal_dominance();

    //If not symmetric it can't be symmetric positive definite
    if (!sym){sym_pd = false;}

    //If symmetric and diagonally dominant we know it is also
    //symmetric positive definite
    else if (sym && dd){sym_pd = true;}

    /*If symmetric and not diagonally dominant, we use Sylvester's criterion:
    A symmetric matrix is pd if all its leading principals are positive
    i.e. the determinants of all leading sub-matrices are positive*/
    else{
        double det = 0;

        //Loop over all sub-matrices of A where A is nxn
        for (int n=0; n < spa.get_nr_eqs(); n++){
            std::vector<int> skip;
            det = compute_det(n, 0, 0, skip);
            if (det <= 0){
                sym_pd = false;
                break;
            }
        }
    }
    return sym_pd;
}

/*Compute the determinant of nxn submatrix starting at (start_row, start_col) using recursive
method as implied by Laplce's method for determinant calculation
Skip any columns contained in skip vector. Algorithm iterates over the rows of sub-matrix */
double SparseMat::compute_det(int n, int start_row, int start_col, std::vector<int> skip) const{
    double det = 0;
    std::vector<int> new_skip;

    if (start_row + (n - 1) > spa.get_nr_eqs()){
        std::cerr << "Supplied entries not a valid sub-matrix" << std::endl;
        exit(EXIT_FAILURE);
    }

    //Go along entries in start_row contained in sub-matrix
    int start_index = spa.twoD_to_oneD(start_row, start_col, true);
    int end_index = spa.twoD_to_oneD(start_row, start_col + n, true);

    //If n - number of skipped columns is one, determinant is just the entry in
    //row start_row for which the respective column is not contained in skip vector
    if ((n - skip.size()) == 1){

        for (int j = start_index; j < end_index; j++){
            if((!skip.empty()) && (std::find(skip.begin(), skip.end(), spa.get_col_index(j)) != skip.end())){continue;}

            else{
                det = entries[j];
                break;
            }
        }
    }

    //Otherwise
    else{
        int cofactor = 1;

        //Go along first row in submatrix
        for (int i = start_index; i < end_index; i++){
            //If current entry corresponds to column to be skipped, continue
            if ((!skip.empty()) && (std::find(skip.begin(), skip.end(), spa.get_col_index(i)) != skip.end())){continue;}

            //Otherwise, add to determinant current entry times cofactor and determinant of new submatrix.
            //New submatrix must have column of current entry removed (skipped)
            else{
                new_skip = skip;
                new_skip.push_back(spa.get_col_index(i));

                //Compute cofactor from column index in actual matrix
                det += cofactor*entries[i]*compute_det(n, start_row + 1, start_col, new_skip);

                //next non-skipped entry in row will have cofactor opposite to current one
                cofactor *= -1;
            }
        }
    }

    return det;
}

//Multiply a sparse matrix by a vector
Eigen::VectorXd SparseMat::operator*(const Eigen::VectorXd& v) const{
    if (v.size() != spa.get_nr_eqs()){
        std::cerr << "Matrix column number not same as vector entry number" << std::endl;
        exit(EXIT_FAILURE);
    }

    Eigen::VectorXd result = Eigen::VectorXd::Zero(spa.get_nr_eqs());

    for (int i = 0; i < spa.get_nr_eqs(); i++){
        double sum = 0;
        int col;

        //Iterate over rows when multiplying by vector
        for (int j = spa.get_row_index(i); j < spa.get_row_index(i+1); j++){
            col = spa.get_col_index(j);
            sum += entries[j]*v(col);
        }

        result(i) = sum;
    }

    return result;
}

void SparseMat::reset(){
    std::fill(entries.begin(), entries.end(), 0.0);
}

