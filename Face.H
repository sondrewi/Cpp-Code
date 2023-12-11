#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>

// Define Face class
class Face {
    public:
        Eigen::ArrayXi vertices;
        int ownerCell;
        int neighbourCell;
        Eigen::Vector3d Sf;
        Eigen::Vector3d centroid;
        double f_x;
        double delta;
        int boundary;
        bool empty_face;
};
