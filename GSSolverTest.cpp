/*Test file for GSSolver. Here, we solve the convection-diffusion equation over a 3m*0.5m domain
for some passive scalar (such as temperature). We assume the scalar is uniformaly 1 over the whole
domain at time 0. We assume constant dirichlet boundary of 0 at the left side of the domain, and a
zero-gradient Neumann condition at the four other boundaries. A source term appears at the location
(0.5, 0.0). In the below example it is activated for 5 time-steps every 20th time-step.
Note that while the problem description is 2D, the actual mesh used is 3D, but there is only one
cell in the z-direction across the entire mesh*/

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <Eigen/Dense>
#include "GSSolver.H"

std::string points_file = "/Users/sondrew/Documents/University/Cambridge/Courses/IFD/mesh_files/points_m.dat";
std::string faces_file = "/Users/sondrew/Documents/University/Cambridge/Courses/IFD/mesh_files/faces_m.dat";
std::string cells_file = "/Users/sondrew/Documents/University/Cambridge/Courses/IFD/mesh_files/cells_m.dat";
std::string boundary_file = "/Users/sondrew/Documents/University/Cambridge/Courses/IFD/mesh_files/boundary_m.dat";

int main(void) {
    //Initiate mesh, sparse addressing object, sparse matrix object and sparse linear system
    Mesh mesh(points_file, faces_file, cells_file, boundary_file);
    SparseAddress SA(mesh);
    SparseMat A(SA);
    SLS SparseSystem(A);

    //Define boundary conditions. Note that a boundary condition type 0
    //is defined to correspond to Dirichlet type while 1 will correspond to Neumann
    //(see SLS.cpp set_diffusion_mat function for use of this convention)
    Eigen::MatrixXd boundary_values = Eigen::MatrixXd::Zero(3, 4);

    std::vector<int> bc_type; bc_type.resize(4);
    bc_type[0] = 0; bc_type[1] = 1;
    bc_type[2] = 1; bc_type[3] = 1;

    //Set the uniform velocity field
    Eigen::Vector3d u(1.0,0.0,0.0);

    //Initiate the smoother and solver objects
    //carry out 3 sweeps per solver iteration
    //maximum 20 solver iterations
    GSSmooth smoother(SparseSystem);

    GSSolver solver(smoother, 3, 20);

    //Set up large matrix of solutions. We solve for 1000 timesteps in addition to an initial solution
    Eigen::MatrixXd solutions = Eigen::MatrixXd::Zero(mesh.getNumCells(), 1001);

    for (int i = 0; i < mesh.getNumCells(); i++){
        mesh.set_cell_scalar(i, 1);
        solutions(i, 0) = 1;
    }

    //find cell whose centroid is closest to (0.5,0,0)
    Eigen::Vector3d source_loc(0.5,0.0,0.0);
    int closest = 0;
    double closest_dist = (mesh.get_cell_centroid(0) - source_loc).norm();

    for (int i = 1; i < mesh.getNumCells(); i++){
        if ((mesh.get_cell_centroid(i) - source_loc).norm() < closest_dist){
            closest = i;
            closest_dist = (mesh.get_cell_centroid(i) - source_loc).norm();
        }
    }

    //Define diffusion parameter gamma and time step dt
    double gamma = 0.01, dt = 0.01;

    //Solve system for 1001 timesteps
    for (int i = 1; i < 1001; i++){
        //Reset x and b vector in sparse system to zero
        SparseSystem.reset_x(); SparseSystem.reset_b();

        /*Also reset entries in A matrix (note that while the A matrix will not actually change between time-steps
        in this particular example, this will not be the case for other applications such as with moving meshes,
        for which the size and entries for A may change significantly*/
        A.reset();

        mesh.calc_gradients(bc_type, boundary_values);

        //Discretise temporal, convection and diffusion operators as described in SLS.cpp
        //and add discretisations to A and b as appropriate
        SparseSystem.set_temporal_mat(dt);
        SparseSystem.set_diffusion_mat(gamma, boundary_values, bc_type);
        SparseSystem.set_convection_mat(u, boundary_values, bc_type, true);

        //Account for intermittent source term
        if (i%20 > 15){
            SparseSystem.add_to_b_entry(closest, 400*mesh.get_cell_vol(closest));
        }

        SparseSystem.set_Ax();
        SparseSystem.set_res();

        solver.solve(0.00000001);

        for (int j = 0; j < mesh.getNumCells(); j++){
            if (abs(SparseSystem.get_x_entry(j)) < 1e-9){
                mesh.set_cell_scalar(j, 0);
                solutions(j, i) = 0;
            }

            else{
                mesh.set_cell_scalar(j, SparseSystem.get_x_entry(j));
                solutions(j, i) = SparseSystem.get_x_entry(j);
            }
        }
    }

    //Output file
    std::ofstream output("Practical6/CD_OOP.dat");

    for (int i = 0; i < mesh.getNumCells(); i++){
        output << mesh.get_cell_centroid(i)(0) << "  " << mesh.get_cell_centroid(i)(1) << "  ";
        for (int k=0; k < 1000; k++){
            output << solutions(i, k) << "  ";
        }
        output << solutions(i, 1000) << std::endl;
    }

}
