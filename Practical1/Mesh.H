//Header file for Mesh.cpp file

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <fstream>
#include <Eigen/Dense>
#include <stdexcept>
#include <cstdlib>
#include "BoundaryPatch.H"
#include "Cell.H"
#include "Face.H"
#include "Point.H"

class Mesh {
    private:
        std::vector<Point> points;
        std::vector<Face> faces;
        std::vector<Cell> cells;
        std::vector<BoundaryPatch> boundaryPatches;

        void readPoints(const std::string& fileName);
        void readFaces(const std::string& fileName);
        void readCells(const std::string& fileName);
        void readBoundaryPatches(const std::string& fileName);
        void assign_face_owners();
        void assign_cell_neighbours();
        void assign_boundary_cells();
        Eigen::VectorXd getCellVolumes() const;
        Eigen::Vector3d face_raw_avg(int face_nr) const;
        void Set_Sf(int face_nr);
        std::vector<std::array<int, 2> > triangle_list(int face_nr) const;
        std::vector<Eigen::Vector3d> triangle_midpoints(int face_nr) const;
        Eigen::VectorXd triangle_areas(int face_nr) const;
        void set_face_centroid(int face_nr);
        Eigen::Vector3d cell_raw_avg(int cell_nr) const;
        Eigen::VectorXd pyramid_volumes(int cell_nr) const;
        void set_cell_vol(int cell_nr);
        void set_cell_centroid (int cell_nr);
        void set_face_interpolation_factor(int face_nr);
        void set_delta_coef(int face_nr);

    public:
        // Constructor to read mesh data from files
        Mesh(const std::string& pointsFile, const std::string& facesFile, const std::string& cellsFile, const std::string& boundaryFile);
        int getNumPoints() const;
        int getNumFaces() const;
        int getNumCells() const;
        int getNumBoundaryPatches() const;
        double get_cell_vol(int cell_nr) const;
        double get_cell_scalar(int cell_nr) const;
        bool get_face_empty(int face_nr) const;
        int get_face_owner(int face_nr) const;
        int get_face_nbr(int face_nr) const;
        Eigen::Vector3d get_face_normal(int face_nr) const;
        int face_boundary(int face_nr) const;
        Eigen::Vector3d get_cell_centroid(int cell_nr) const;
        double get_face_interpolation_factor(int face_nr) const;
        double get_delta_coef(int face_nr) const;
        Eigen::Vector3d get_cell_grad(int cell_nr) const;
        std::vector<int> get_cell_nbrs(int cell_nr) const;
        void set_cell_scalar(int cell_nr, double value);
        void calc_gradients(std::vector<int>& bc_type, Eigen::MatrixXd& boundary_values);
};
