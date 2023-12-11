

// Define BoundaryPatch class header
class BoundaryPatch {
public:
    std::string name;
    Eigen::ArrayXi faces;
    int cell;
    double bc;
};
