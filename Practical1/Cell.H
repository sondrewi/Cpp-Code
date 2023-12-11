#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>

// Define Cell class header
class Cell {
public:
    Eigen::ArrayXi faces;
    std::vector<int> nbrs;
    Eigen::Vector3d centroid;
    double vol;
    double T;
    Eigen::Vector3d grad;
};
