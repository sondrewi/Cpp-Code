#include "GSSolver.H"

//Constructor for Gauss-Seidel solver. Specify number of sweeps between each time that residual is evaluated and max number of iterations
GSSolver::GSSolver(GSSmooth& smoother_object, int sweeps_, int maxit) : A(smoother_object.A), mesh(smoother_object.mesh), spa(smoother_object.spa), sls(smoother_object.sls), smoother(smoother_object){
    sweeps = sweeps_;
    max_it = maxit;
}

void GSSolver::solve(double relTol){
    int iter_count = 0;
    double prev_resid_norm = sls.calc_res_norm(2);
    double b_norm = sls.calc_b_norm(2);

    std::cout << "Norm res: " << prev_resid_norm << std::endl;

    //While relative tolerance time b-norm is smaller than residual norm or max iterations not exceeded,
    //carry out a new sets of sweeps of Gauss-Seidel method
    while(b_norm * relTol < sls.calc_res_norm(2) && iter_count < max_it){
        iter_count += 1;

        smoother.smooth(sweeps);

        //Set new residual and Ax vector based on current guess x
        sls.set_Ax();
        sls.set_res();
        std::cout << "Reduction in Norm: " << prev_resid_norm - sls.calc_res_norm(2) << std::endl;
        prev_resid_norm = sls.calc_res_norm(2);
    }
}
