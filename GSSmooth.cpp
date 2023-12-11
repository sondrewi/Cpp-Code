#include "GSSmooth.H"

//Constructor for Gauss-Seidel smoother
GSSmooth::GSSmooth(SLS& SparseSystem) : sls(SparseSystem), A(SparseSystem.A), mesh(SparseSystem.mesh), spa(SparseSystem.spa) {
}

//Function to carry out predefined number of sweeps of Gauss-Seidel algorithm
void GSSmooth::smooth(int nSweeps){
    int sweep_count = 0;

    while (sweep_count < nSweeps){
        sweep_count += 1;

        //Iterate over the cells (equations)
        for (int i=0; i < mesh.getNumCells(); i++){
            double sum = 0;
            double diag;

            /*In the sparse system Ax = b, let A = (L + D + U), i.e. decomposition
            into lower triangular, diagonal, and upper triangular part. Let x_j
            denote the current guess for the solution vector x (x_j is a vector).
            Let x[i]_(j+1) denote the i-th entry in x_(j+1). Denote by L[i], and U[i]
            the i-th rows of matrices L and U, respectively, and by inv(D)[i,i] the inverse
            of the i-th diagonal entry in D.
            In Gauss-Seidel, the new guess for x is computed by the formula
            x[i]_(j+1) = (inv(D)[i,i])*(b[i] - L[i]*x_(j+1) - U*x_(j).
            For each i, we compute x[i]_(j+1) based only on "updated" elements x[j]_(j+1) for which
            j < i, since L is lower triangular. */

            //Find coefficients belonging to row i
            for (int j = spa.get_row_index(i); j < spa.get_row_index(i+1); j++){

                //If coefficient is on diagonal, store it for later
                if (spa.get_col_index(j) == i){diag = A.get_entry1D(j);}

                //Else, add it to sum L[i]*x_(j+1) + U*x_(j)
                else {sum += A.get_entry1D(j)*sls.get_x_entry(spa.get_col_index(j));}
            }

            //Compute x[i]_(j+1)
            sls.set_x_entry(i, (1/diag)*(sls.get_b_entry(i) - sum));
        }
    }
}
