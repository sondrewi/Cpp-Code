#include "Mesh.H"

/* This file defines a unstructured mesh class to be used for CFD purposes.
At the base level, the mesh is defined by a set of points (defined in supplied pointsFile)
Faces are made up of specific sets of cells as defined in the facesFile
Volumetric cells are made up of a set of faces as defined in cellsFile
Boundary patches are made up of closed sets of faces (as defined in boundaryFile)

The Mesh class communicates with four other classes point, face, cell, and Boundarypatch

All grids generated will be inherently 3D. We can obtain equivalents
for 1D and 2D models by restricting the number of cells in the z-direction (2D) and y-direction (1D) to 1
Faces pointing in these directions will be designated as "empty" (no boundary condition necessary)

Note that a significant amount of parsing is necessary when reading the files, as the format is
string-based (.dat files)
*/

// Constructor to read mesh data from files
Mesh::Mesh(const std::string& pointsFile, const std::string& facesFile, const std::string& cellsFile, const std::string& boundaryFile) {
    //Read files containing specifications of mesh
    readPoints(pointsFile);
    readFaces(facesFile);
    readCells(cellsFile);
    readBoundaryPatches(boundaryFile);

    // Calculate area-scaled normal vectors to cell faces
    // and face centroids
    //Assign these to the relevant face instances
    for (int i = 0; i < faces.size(); i++){
        Set_Sf(i);
        set_face_centroid(i);
    }

    // Calculate cell centroids and cell volumes
    //Assign these to the relevant cell instances
    for (int i = 0; i < cells.size(); i++){
        set_cell_centroid(i);
        set_cell_vol(i);
    }

    //Assign face owner and neighbour cell (i.e. the cells of whose boundary the face is part)
    //We need to distinguish owner and neighbour cells for consistency purposes when solving
    //associated system of equations. Face owner is the cell for which the face normal is outward-pointing
    assign_face_owners();

    //For each cell assign a set of neighbouring cells
    assign_cell_neighbours();

    /* For each boundary patch assign to it an owner cell
    // This function also ensures that face normals are outward
     pointing for faces that form part of a boundary patch */
    assign_boundary_cells();

    /* Face interpolation factor and delta coefficient
    are attributes of each face. Only assign face-interpolation
    factor for non-boundary patches */
    for (int i = 0; i < faces.size(); i++){
        if (faces[i].boundary < 0){set_face_interpolation_factor(i);}
        set_delta_coef(i);
    }

}

/* Function to read points from a file. The pointsFile will always have the form of:
nrPoints (number of points)
(
(p1_x p1_y p1_z)
(p2_x p2_y p2_z)
...
(pn_x pn_y pn_z)
)
*/
void Mesh::readPoints(const std::string& fileName) {
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the points file." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line; //declare string to hold current line
    std::string str_num; //declare string to representation of current number
    std::getline(inputFile, line); //get first line from file
    int numPoints; //will hold number of vertices

    // Get string representation of nr of points
    for (std::size_t i = 0; i < line.length(); ++i) {
        str_num += line[i];
    }

    // convert to integer
    numPoints = std::stoi(str_num);

    //Resize points vector (vector of points objects declared in header file)
    points.resize(numPoints);

    std::getline(inputFile, line); //skip past line containing only parenthesis

    int j; //declare variable that will be used to iterate over characters in line

    for (int i = 0; i < numPoints; ++i){ //iterate over each line containing vertex
        std::getline(inputFile, line);

        j = 1; //Each line starts with a parenthesis that we can skip past
        str_num = "";
        while (line[j] != ' '){
            str_num += line[j]; //append characters of number to string until space is reached
            j += 1;
        }

        points[i].coord[0] = std::stod(str_num); //Assign to given point its x-coordinate

        j += 1; //skip past space character
        str_num = "";
        while (line[j] != ' '){
            str_num += line[j];
            j += 1;
        }

        points[i].coord[1] = std::stod(str_num); //Assign to given point its y-coordinate

        j += 1;
        str_num = "";
        while ((line[j] != ' ') && (line[j] != ')')){
            str_num += line[j];
            j += 1;
        }

        points[i].coord[2] = std::stod(str_num); //Assign to given point its z-coordinate
    }

    inputFile.close();
}

/* Function to read face definitions from a file.
Let k_i denote the number of points that make up the corners of a face
The facesFile will always have the form of:
nrFaces (number of faces)
(
k_0([space separated list of k_0 points (indexed from zero) for face 0])
k_1([space separated list of k_1 points (indexed from zero) for face 1])
...
k_n ([space separated list of k_n points (indexed from zero) for face n])
)

**Note: it is assumed that points are provided in a counterclockwise order
with respect to the face plane
*/
void Mesh::readFaces(const std::string& fileName) {
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the faces file." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line; //declare string to hold current line
    std::string str_num; //declare string to representation of current number
    std::getline(inputFile, line); //get first line from file
    int numFaces; //will hold number of faces
    int numVertices; //will hold number of verteces on current face

    // Get string representation of nr of faces
    for (std::size_t i = 0; i < line.length(); ++i) {
        str_num += line[i];
    }

    //convert string to integer and resize face vector of face objects (decalred in header)
    numFaces = std::stoi(str_num);
    faces.resize(numFaces);

    std::getline(inputFile, line); //skip past line containing only parenthesis

    int j; //declare variable that will be used to iterate over characters in line

    /*Iteration over each line proceeds much as with the points, but now we populate the "vertices"-vector for
    each face object */
    for (int i = 0; i < numFaces; ++i){ //iterate over each line containing list of vertices
        int j=0; //reset j
        std::getline(inputFile, line);
        str_num = "";

        //Get number of points which make up a face in string form
        while((line[j] != '(') && (line[j] != ' ')){
            str_num += line[j];
            j += 1;
        }

        //Convert to integer and resize vertices vector of the given face object
        numVertices = std::stoi(str_num);
        faces[i].vertices.resize(numVertices);

        j += 1; //skip to first character of first vertex label on this face

        for (int k = 0; k < numVertices; ++k) {
            str_num = "";
            while ((line[j] != ' ') && (line[j] != ')')){
                str_num += line[j];
                j += 1;
            }

            faces[i].vertices[k] = std::stoi(str_num);
            j += 1;
        }

        /* Faces contain integer attribute boundary which indicates the index of the
        patch if positive and is -1 if the face is not on boundary of domain

        They also contain boolean attribute empty_face which indicate if face
        is on a boundary that is being considered with the boundary conditions or not
        (only relevant for 1D and 2D models)

        Initially mark all faces as not being on a boundary and not being empty
        This will be changed when boundary patches are initialised and
        when cell-neighbours are assigned, respectively.
        A face is empty if it is at the boundary of mesh but forms part
        of no boundary patch.
        */
        faces[i].boundary = -1;
        faces[i].empty_face = false;
    }

    inputFile.close();
}


/* Function to read cells from a file
Let k_i denote the number of faces that make up the boundaries of a cell
The cellsFile will always have the form of:
nrcells (number of cells)
(
k_0([space separated list of k_0 faces (indexed from zero) for cell 0])
k_1([space separated list of k_1 faces (indexed from zero) for cell 1])
...
k_n ([space separated list of k_n faces (indexed from zero) for cell n])
)
*/
void Mesh::readCells(const std::string& fileName) {
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the cells file." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line; //declare string to hold current line
    std::string str_num; //declare string to representation of current number
    std::getline(inputFile, line); //get first line from file
    int numCells; //will hold number of cells
    int numFaces; //will hold number of verteces on current Cell

    //get number of cells
    for (std::size_t i = 0; i < line.length(); ++i) {
        str_num += line[i];
    }

    //Resize vector that holds cell objects (declared in header file)
    numCells = std::stoi(str_num);
    cells.resize(numCells);

    std::getline(inputFile, line); //skip past line containing only parenthesis

    int j; //declare variable that will be used to iterate over characters in line

    /* Note that we cannot directly assign to faces their owner and neighbour cells at this stage, as we
    do not know what faces belong to boundary patches. Thus we assign to every face an "owner" and "neighbour"
    value of -1 */

    for (int i = 0; i < numCells; ++i){ //iterate over each line containing list of faces for each cell
        j=0; //reset j
        std::getline(inputFile, line);
        str_num = "";

        // find number of faces belonging to current cell
        while((line[j] != '(') && (line[j] != ' ')){
            str_num += line[j];
            j += 1;
        }

        //Resize faces vector of current cell object
        numFaces = std::stoi(str_num);
        cells[i].faces.resize(numFaces);

        j += 1; //skip to first character of first face index on this cell

        //Iterate over number of faces defining current cell
        for (int k = 0; k < numFaces; ++k) {
            str_num = "";
            while ((line[j] != ' ') && (line[j] != ')')){
                str_num += line[j];
                j += 1;
            }

            //Populate faces vector in current cell object with current face index
            //initially define every neighbour and owner by -1
            //which will signify boundary patch
            cells[i].faces[k] = std::stoi(str_num);
            faces[cells[i].faces[k]].ownerCell = -1;
            faces[cells[i].faces[k]].neighbourCell = -1;
            j += 1;
        }
    }

    inputFile.close();
}

/* Function to read boundary patches from a file
Let k_i denote the number of faces that make up a boundary patch
Let name_i denote the name of patch i
The boundaryFile will always have the form of:
nrpatches (number of boundary patches)
(
name_0
k_0
(
[space separated list of k_0 faces (indexed from zero) for patch 0])
)
...
name_n
k_n
(
([space separated list of k_n faces (indexed from zero) for patch n])
)
)
*/
void Mesh::readBoundaryPatches(const std::string& fileName) {
    std::ifstream inputFile(fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the boundary patches file." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line; //declare string to hold current line
    std::string str_num; //declare string to representation of current number
    int numBoundaries; //will hold number of Boundary patches
    int numFaces; //will hold number of verteces on current Boundary patch

    std::getline(inputFile, line); //get first line from file

    //get number of boundary patches
    for (std::size_t i = 0; i < line.length(); ++i) {
        str_num += line[i];
    }

    numBoundaries = std::stoi(str_num);
    boundaryPatches.resize(numBoundaries);

    std::getline(inputFile, line); //jump to line containing only parenthesis

    int j; //declare variable that will be used to iterate over characters in line

    for (int i = 0; i < numBoundaries; ++i){ //iterate over each block containing name of patch and info on patches
        j=0; //reset j
        std::getline(inputFile, line); //jump to line containing name

        boundaryPatches[i].name = line; //assign name

        std::getline(inputFile, line); //jump to line containing number of faces adjoining

        str_num = "";

        // find number of faces adjoining current boundary patch
        while(line[j]){
            str_num += line[j];
            j += 1;
        }

        numFaces = std::stoi(str_num);
        boundaryPatches[i].faces.resize(numFaces);

        std::getline(inputFile, line); //jump to line containing which faces are adjoining
        std::getline(inputFile, line);
        j = 0; //Reset j

        //for each boundary patch, assign to its faces vector the right face labels
        for (int k = 0; k < numFaces; ++k) {
            str_num = "";
            while ((line[j] != ' ') && (line[j] != ')') && (j < line.length())){
                str_num += line[j];
                j += 1;
            }
            boundaryPatches[i].faces[k] = std::stoi(str_num);

            //Assign to each face adjoining patch the patch number
            //Boundary is a
            faces[std::stoi(str_num)].boundary = i;

            j += 1;
        }
        std::getline(inputFile, line);
    }

    inputFile.close();
}

void Mesh::assign_face_owners(){
    //Need to loop over cells to assign face owners and neighbours
    //as faces do not point to cells in any way
    int numCells = cells.size();
    Eigen::ArrayXi cell_faces;
    Eigen::Vector3d f_normal, f_centroid,  Sf, c_centroid;
    double distanceSq_add_normal, distanceSq_sub_normal;
    /*Declare vectors that will hold face normal, area-normalised centroid,
    area-scaled normal, cell-centroid for adjoining cell. Denote by cf the vector
    from cell centroid to face centroid. We also declare here variables to hold the
    norms |cf + face-normal| and |cf - face-normal|.  */

    for (int i = 0; i < numCells; i++){

        //list of faces for each cell
        Eigen::ArrayXi cell_faces = cells[i].faces;
        int numFaces_cell = cell_faces.size();

        //initiate varibales for current cell
        cell_faces = cells[i].faces;
        c_centroid = cells[i].centroid;

        /*For each face adjoining a cell, if face normal ADDED to cf is larger than
        face normal SUBTRACTED from cf then the given cell is assigned
        as owner. Otherwise, the given cell is assigned as neighbour */
        for (int j = 0; j < numFaces_cell; j++){
            //Get area-scaled face normal
            Eigen::Vector3d Sf = faces[cell_faces[j]].Sf;

            //Get face normal & centroid
            f_normal = Sf/Sf.norm();
            f_centroid = faces[cell_faces[j]].centroid;

            //Compute |cf + face-normal| and |cf - face-normal|
            distanceSq_add_normal = (f_centroid - c_centroid + f_normal).norm();
            distanceSq_sub_normal = (f_centroid - c_centroid - f_normal).norm();

            if (distanceSq_add_normal > distanceSq_sub_normal){
                faces[cell_faces[j]].ownerCell = i;
            }

            else {
                faces[cell_faces[j]].neighbourCell = i;
            }
        }
    }

    /* The above is inconsistent with each face having an owner cell. In order to
    ensure this, we loop over all faces and let owner <- neighbour if owner is -1*/
    //Now make sure each face has an owner. That is, if a face only has a neighbour
    //cell but no owner, make the neighbour the owner and
    for (int j = 0; j < faces.size(); j++){
        if (faces[j].ownerCell < 0){
            faces[j].ownerCell = faces[j].neighbourCell;
            faces[j].neighbourCell = -1;
        }
    }
}

void Mesh::assign_cell_neighbours(){
    //loop over faces and couple face owner/neighbours as cell neighbours
    for (int i = 0; i < faces.size(); i++){
        int cell1, cell2;
        cell1 = faces[i].ownerCell;
        cell2 = faces[i].neighbourCell;

        if ((cell2 < 0) && (faces[i].boundary < 0)){
            //if face is on boundary and is not already assigned to a patch
            //(done during parsing), mark it as empty face.
            faces[i].empty_face = true;
        }

        // Else face marks boundary between two cells. Alter the cells'
        // neighbour vectors accordingly
        else if (faces[i].boundary < 0){
            cells[cell1].nbrs.push_back(cell2);
            cells[cell2].nbrs.push_back(cell1);
        }
    }
}

void Mesh::assign_boundary_cells(){
    //Loop over boundary patches and assign adjoining cell by owner cell of first face
    int a_face, face_nr, cell_nr, numFaces;
    Eigen::Vector3d cell_cent, face_cent, Sf, f_normal;
    double distanceSq_add_normal, distanceSq_sub_normal;

    for (int i = 0; i < boundaryPatches.size(); i++){
        a_face = boundaryPatches[i].faces[0];
        boundaryPatches[i].cell = faces[a_face].ownerCell;

        //also make sure that Sf (area-scaled normal) is outward pointing for every
        //face adjoining a boundary patch:
        numFaces = boundaryPatches[i].faces.size();

        for (int j=0; j < numFaces; j++){
            /* Procedure is the same as when owner/neighbour cells of faces
             were assigned. If Sf is inward pointing, we have
             |face_cent - cell_cent + f_normal| < |face_cent - cell_cent - f_normal|
            */
            face_nr = boundaryPatches[i].faces[j];
            cell_nr = faces[face_nr].ownerCell;
            cell_cent = cells[cell_nr].centroid;
            face_cent = faces[face_nr].centroid;
            Sf = faces[face_nr].Sf;
            f_normal = Sf/Sf.norm();

            distanceSq_add_normal = (face_cent - cell_cent + f_normal).norm();
            distanceSq_sub_normal = (face_cent - cell_cent - f_normal).norm();

            if (distanceSq_add_normal < distanceSq_sub_normal){
                faces[face_nr].Sf *= -1;
            }
        }
    }
}

// Add getter functions for various mesh properties
int Mesh::getNumPoints() const {
    return points.size();
}

int Mesh::getNumFaces() const {
    return faces.size();
}

int Mesh::getNumCells() const {
    return cells.size();
}

int Mesh::getNumBoundaryPatches() const {
    return boundaryPatches.size();
}

std::vector<int> Mesh::get_cell_nbrs(int cell_nr) const{
    return cells[cell_nr].nbrs;
}

double Mesh::get_cell_vol(int cell_nr) const {
    return cells[cell_nr].vol;
}

double Mesh::get_cell_scalar(int cell_nr) const {
    return cells[cell_nr].T;
}

bool Mesh::get_face_empty(int face_nr) const {
    return faces[face_nr].empty_face;
}

int Mesh::get_face_owner(int face_nr) const {
    return faces[face_nr].ownerCell;
}

int Mesh::get_face_nbr(int face_nr) const {
    return faces[face_nr].neighbourCell;
}

Eigen::Vector3d Mesh::get_face_normal(int face_nr) const{
    return faces[face_nr].Sf;
}

int Mesh::face_boundary(int face_nr) const {
    return faces[face_nr].boundary;
}

Eigen::Vector3d Mesh::get_cell_centroid(int cell_nr) const {
    return cells[cell_nr].centroid;
}

double Mesh::get_face_interpolation_factor(int face_nr) const {
    return faces[face_nr].f_x;
}

double Mesh::get_delta_coef(int face_nr) const {
    return faces[face_nr].delta;
}

Eigen::Vector3d Mesh::get_cell_grad(int cell_nr) const {
    return cells[cell_nr].grad;
}

void Mesh::set_cell_scalar(int cell_nr, double value){
    cells[cell_nr].T = value;
}

//Return volume of all cells
Eigen::VectorXd Mesh::getCellVolumes() const {
    int numCells = cells.size();
    Eigen::VectorXd volumes(numCells);

    for (int i = 0; i < numCells; i++){
        volumes[i] = cells[i].vol;
    }

    return volumes;
}

//face_raw_avg is used to find the raw average of the points
//defining face corners
Eigen::Vector3d Mesh::face_raw_avg(int face_nr) const{
    //find number of vertices for current face
    int numVertices = faces[face_nr].vertices.size();

    //declare vector that will hold RAW AVERAGE of all vertices
    //and integer variable that will hold a point index
    Eigen::Vector3d raw_avg(0.0,0.0,0.0);
    int vertex_lable;

    //compute raw average
    for (int i = 0; i < numVertices; i += 1){
        vertex_lable =  faces[face_nr].vertices[i];
        raw_avg += (points[vertex_lable].coord / numVertices);
    }

    return raw_avg;
}

/*Set Sf is used to compute face normal scaled by face area. Face will not be a perfect
plane due to numerical inaccuracy. Hence, calculate face area by summing areas of
triangles formed by adjacent points and face raw average (on a given face). The DIRECTION
of Sf is found by summing the normal vectors of said triangles. This sum is then normalised
and subsequently scaled by the total area*/
void Mesh::Set_Sf(int face_nr){
    std::vector<std::array<int, 2> > tri_list = triangle_list(face_nr);
    Eigen::Vector3d raw_avg = face_raw_avg(face_nr);
    Eigen::Vector3d Sf(0,0,0);

    //Variable to contain total face area
    double areas = 0;

    //Define two "vectors" that will be used to compute normal to
    //each triangle of the face
    //In practice just subtract third vertex from first and second
    Eigen::Vector3d vec1, vec2;

    /*tri_list contains sets of pairs of adjacent points/face corners
    two points will only form a pair if they are adjacent when going
    round the face counter-clockwise*/
    for (int i = 0; i < tri_list.size(); i++){
        //Vectors from each of points in a pair to raw face average
        vec1 = points[tri_list[i][0]].coord - raw_avg;
        vec2 = points[tri_list[i][1]].coord - raw_avg;

        //Vector normal to triangle formed by pair and face raw avg
        Sf += vec1.cross(vec2);

        //Add to total area of face
        areas += vec1.cross(vec2).norm()/2;
    }

    Sf.normalize();
    Sf = areas * Sf;

    //Set Sf attribute of given face
    faces[face_nr].Sf = Sf;
}

/*triangle_list produces a vector of pairs of adjacent points/face corners.
Two points will only form a pair if they are adjacent when going
round the face counter-clockwise*/
std::vector<std::array<int, 2> > Mesh::triangle_list(int face_nr) const{
    Eigen::Vector3d raw_avg = face_raw_avg(face_nr);

    //find number of vertices for current face
    int numVertices = faces[face_nr].vertices.size();

    // declare array that will hold list of vertices for non-overlapping triangles for face
    // where vertex-points are identified by their indices. Only two vertices are
    // listed for each triangle since the third will always be the raw_avg
    std::vector<std::array<int, 2> > tri_list; tri_list.resize(numVertices);

    for(int i = 0; i < numVertices - 1; i++){
        tri_list[i][0] = faces[face_nr].vertices[i];
        tri_list[i][1] = faces[face_nr].vertices[i+1];
    }

    //Last vertex must be coupled with first
    tri_list[numVertices - 1][0] = faces[face_nr].vertices[numVertices - 1];
    tri_list[numVertices - 1][1] = faces[face_nr].vertices[0];

    return tri_list;
}

//Find midpoint of each triangle of which face is composed.
//This function will be used when calculating face centroid
std::vector<Eigen::Vector3d> Mesh::triangle_midpoints(int face_nr) const{
    //Last vertex is always the raw average of all vertices
    std::vector<std::array<int, 2> > tri_list = triangle_list(face_nr);
    int numVertices = faces[face_nr].vertices.size();
    Eigen::Vector3d raw_avg = face_raw_avg(face_nr);
    std::vector<Eigen::Vector3d > mpoints; mpoints.resize(numVertices);

    //midpoint of each triangle is just average of vertices
    for (int i = 0; i < numVertices; i++){
        mpoints[i] = (points[tri_list[i][0]].coord + points[tri_list[i][1]].coord + raw_avg)/3.0;
    }

    return mpoints;
}

void Mesh::set_face_centroid(int face_nr) {
    //Compute face centroid by taking weighted mean of midpoints
    Eigen::Vector3d centroid(0.0,0.0,0.0);
    std::vector<std::array<int, 2> > tri_list = triangle_list(face_nr);
    int numVertices = faces[face_nr].vertices.size();
    std::vector<Eigen::Vector3d > midpoints = triangle_midpoints(face_nr);
    Eigen::Vector3d raw_avg = face_raw_avg(face_nr);

    //Define two "vectors" that will be used to compute normal to
    //each triangle of the face
    //In practice just subtract third vertex from first and second
    Eigen::Vector3d vec1, vec2;

    //Define variable for area of given triangle in face and variable for
    //total area
    double tri_area, total_area = 0;

    //Loop over pairs in tri_list and calculated midpoints to compute
    //weighted average
    for (int i = 0; i < numVertices; i++){
        //Vectors from each of points in a pair to raw face average
        vec1 = points[tri_list[i][0]].coord - raw_avg;
        vec2 = points[tri_list[i][1]].coord - raw_avg;

        //area of triangle
        tri_area = vec1.cross(vec2).norm();
        total_area += tri_area;

        centroid += (tri_area*midpoints[i]);
    }

    centroid /= total_area;

    faces[face_nr].centroid = centroid;
}

//Find raw average of face centroids around cell
//This is just a simple average where we loop over faces
//adjoining cell
Eigen::Vector3d Mesh::cell_raw_avg(int cell_nr) const{
    int numFaces = cells[cell_nr].faces.size();
    Eigen::ArrayXi cell_faces = cells[cell_nr].faces;
    Eigen::Vector3d raw_avg(0.0,0.0,0.0), face_cent;

    for (int i = 0; i < numFaces; i++){
        face_cent = faces[cell_faces[i]].centroid;
        raw_avg += face_cent;
    }

    raw_avg /= numFaces;

    return raw_avg;
}

/* Cell centroid will be calculated as the average of face-centroids,
weighted by the volume of the pyramid formed between face and cell raw average.
The given function returns a vector pyramid volumes for a given cell, indexed in
order of the face indices for that cell*/
Eigen::VectorXd Mesh::pyramid_volumes(int cell_nr) const{
    int numFaces = cells[cell_nr].faces.size();
    Eigen::ArrayXi cell_faces = cells[cell_nr].faces;
    Eigen::Vector3d raw_avg = cell_raw_avg(cell_nr);
    Eigen::VectorXd pyr_vols(numFaces);

    Eigen::Vector3d face_cent;
    double pyramid_height;

    for (int i = 0; i < numFaces; i++){
        face_cent = faces[cell_faces[i]].centroid;
        pyramid_height = (face_cent - raw_avg).norm();

        //Volume of pyramid for which base is polygon is calculated using standard formula
        pyr_vols[i] = (1.0/3.0)*pyramid_height*(faces[cell_faces[i]].Sf.norm());
    }

    return pyr_vols;
}

//Cell volume is stored as attribute of each cell. Computed
//simply as sum of pyramid volumes for given cell
void Mesh::set_cell_vol(int cell_nr) {
    double cell_volume = 0.0;
    int numFaces = cells[cell_nr].faces.size();
    Eigen::VectorXd pyr_volumes(numFaces);
    pyr_volumes = pyramid_volumes(cell_nr);

    for (int i=0; i<numFaces; i++){
        cell_volume += pyr_volumes[i];
    }

    cells[cell_nr].vol = cell_volume;
}

/* Cell centroid will be calculated as the average of face-centroids,
weighted by the volume of the pyramid formed between face and cell raw average.*/
void Mesh::set_cell_centroid (int cell_nr) {
    int numFaces = cells[cell_nr].faces.size();
    Eigen::VectorXd pyr_volumes = pyramid_volumes(cell_nr);
    Eigen::Vector3d the_centroid(0.0,0.0,0.0);
    Eigen::ArrayXi cell_faces;
    cell_faces = cells[cell_nr].faces;
    Eigen::VectorXd f_centroid;
    double cell_volume = 0.0;

    //Compute weighted average
    for (int i=0; i<numFaces; i++){
        f_centroid = faces[cell_faces[i]].centroid;
        the_centroid += pyr_volumes[i]*(f_centroid);
        cell_volume += pyr_volumes[i];
    }

    the_centroid /= cell_volume;
    cells[cell_nr].centroid = the_centroid;
}

/*Set face interpolation factor for given face. When interpolating cell variable ø onto
faces, we use formula (f_x)*ø_ownerCell + (1 - f_x)*ø_nbrCell, where f_x is interpolation
factor. Let Nf denote vector between neighbour cell centroid and face centroid. Let PN
denote vector between owner cell centroid and neighbour cell centroid.
We then set f_x = |Nf|/|PN|. We will only invoke this function for faces not on boundary */
void Mesh::set_face_interpolation_factor(int face_nr){
    int cellP = faces[face_nr].ownerCell;
    int cellN = faces[face_nr].neighbourCell;

    //Get centroids
    Eigen::Vector3d face_cent = faces[face_nr].centroid;
    Eigen::Vector3d cellP_cent = cells[cellP].centroid;
    Eigen::Vector3d cellN_cent = cells[cellN].centroid;

    //Calculate as outlined above
    double f_x = (face_cent - cellN_cent).norm()/(cellP_cent - cellN_cent).norm();
    faces[face_nr].f_x = f_x;
}

/*Set delta coefficient for given face. If face is not on boundary,
delta coefficient is simply 1/(distance between adjoining cell centres).
Otherwise, delta coefficient is 1/(distance from owner cell centroid to face centroid)*/
void Mesh::set_delta_coef(int face_nr){
    int cellP = faces[face_nr].ownerCell;
    int cellN = faces[face_nr].neighbourCell;

    //If face is on boundary (has no neighbour cell)
    if (cellN < 0){
        Eigen::Vector3d cellP_cent = cells[cellP].centroid;
        Eigen::Vector3d face_cent = faces[face_nr].centroid;
        double delta = 1/(cellP_cent - face_cent).norm();
        faces[face_nr].delta = delta;
    }

    //If face is not on boundary of domain. Note, empty faces
    //included here, although their delta coefficient will not be used
    else{
        Eigen::Vector3d cellP_cent = cells[cellP].centroid;
        Eigen::Vector3d cellN_cent = cells[cellN].centroid;
        double delta = 1/(cellP_cent - cellN_cent).norm();
        faces[face_nr].delta = delta;
    }
}

/*We calculate gradients at a cell P by fitting a plane:
Define by d_N the vector going from centroid of P to the
centroid of its neighbour N. Then the error when using the
gradient to interpolate ø_N is
e_N = ø_N - (ø_P + d_N dot (grad(ø))_P)

We seek to minimise the sum of square weighted errors for the
neighbours of P:
e^2_P = sum_N((w_N e_N)^2), where w_N = 1/|d_N|.
The weight is included so that we assign greater importance to
closer neighbours (it is more important that the interpolation onto those
neighbours is better using the gradient since they are closer)

This problem is equivalent to minimising the squared 2-norm of
(M (grad(ø))_P - v) where M is the matrix whose rows are the vectors w_N*d_N
and v is the vector whose N-th entry is w_N*(ø_N - ø_P)

Note that for a boundary cell P, we define the problem in an equivalent way,
with d_N denoting the vector from cell centroid to boundary face centroid
and ø_N denotes the value ø_b at the boundary either taken directly (Dirichlet)
or interpolated (Neumann)

For a cell in 3-dimensional mesh, the matrix M will
have rank at least 3 (due to convexity of cells). We know that the
unique minimum in (grad(ø))_P of the squared 2-norm of (M (grad(ø))_P - v)
can then be found as (grad(ø))_P = inv(M^T M) M^T b, where inv(M^T M) M^T is
the Moore-Penrose pseudo-inverse.

In a "2-dimensional" mesh, we can ignore the z-gradient. For a non-boundary
cell, M will have two columns and rank of at least 2. Gradient can then be
calculated in an equivalent way.

In 1 dimension M will be a vector and the problem is trivial.
*/
void Mesh::calc_gradients(std::vector<int>& bc_type, Eigen::MatrixXd& boundary_values){
    int nr_cells = cells.size();
    Eigen::Matrix3d G, G_inv;
    Eigen::MatrixXd G_sub, G_sub_inv, A;
    Eigen::Vector3d d_N, gradient, face_normal; //may need face normal for boundary face intensive value calculation
    Eigen::ArrayXi cell_faces;
    int N, boundary; //index of neighbour cell and of boundary patch if a face is at boundary
    double ø_N, ø_i;

    //Loop over cells to calculate gradient at each centroid
    for (int i = 0; i < nr_cells; i++){
        //ø_P will stay fixed throughout calculation gradient at given cell i
        //(intensive value at cell centroid)
        ø_i = cells[i].T;

        G = Eigen::Matrix3d::Zero(3,3);
        gradient = Eigen::Vector3d::Zero(3);

        //Extract the list of faces for a cell
        cell_faces = cells[i].faces;


        //Note that the matrix-matrix product M^T M above will take the form
        //G = M^T M = sum_N(w^2_N d_N (d_N)^T)
        //Loop over faces
        for (int j = 0; j < cell_faces.size(); j++){
            if(!faces[cell_faces[j]].empty_face){

                //If a face is not at boundary, d_N is distance between cell centroids
                //d_N computed by subtracting current cell centroid (index i) from neighbouring centroid
                //We account for this by multiplying d_N by factor (-1)^(nbr == i)
                if (faces[cell_faces[j]].boundary < 0){
                    if (i == faces[cell_faces[j]].ownerCell){N = faces[cell_faces[j]].neighbourCell;}
                    else {N = faces[cell_faces[j]].ownerCell;}

                    d_N = cells[N].centroid - cells[i].centroid;
                }

                //If face is at boundary, d_N is distance from cell centroid to face centroid
                else{
                    d_N = faces[cell_faces[j]].centroid - cells[i].centroid;
                }

                //Add to matrix G = M^T M
                G += pow(d_N.norm(), 2)*(d_N * d_N.transpose());
            }
        }


        //Here, we deal with the issue that G may be singular, so that we only invert
        //the submatrix of nonzero entries
        std::vector<int> nonZeroRows;

        //Find non-zero rows of G
        for (int k=0; k < 3; k++){
            if (G(k,k) > 1e-10){
                nonZeroRows.push_back(k);
            }
        }

        //Suppress small entries of G to zero
        for (int k=0; k < 3; k++){
            for (int l=0; l < 3; l++){
                if (G(k,l) < 1e-12){
                    G(k,l) = 0;
                }
            }
        }

        //We will find the invertible sub-matrix of G, G_sub
        //by making the multiplication G_sub = AGA^T
        A = Eigen::MatrixXd::Zero(nonZeroRows.size(), 3);

        //Matrix A is constructed as follows from a zero 3x3 matrix:
        //For k in (0,1,2):
        //  If the kth row of G is non-zero, let A(k,k) = 1
        //Delete all rows of A consisting entirely of zeros
        for (int k=0; k < nonZeroRows.size(); k++){
            A(k, nonZeroRows[k]) = 1.0;
        }

        //Calculate the sub-matrix G_sub
        G_sub = (A*G)*(A.transpose());

        //Cheap inversion of G_sub as it always has rank <= 3
        G_sub_inv = G_sub.inverse();

        //Extend G_sub to size 3x3 by undoing truncation
        G_inv = (A.transpose()*G_sub_inv)*A;

        //Again loop over faces to compute inv(M^T M)(M^T)b = inv(M^T M)sum_N(w^2_N * (ø_N - ø_P) * d_N)
        for (int j = 0; j < cell_faces.size(); j++){
            if(!faces[cell_faces[j]].empty_face){
                //If a face is not at boundary, d_N is distance between cell centroids
                //d_N should subtract current cell centroid (index i) from neighbouring centroid
                //We account for this by multiplying d_N by factor (-1)^(nbr == i)
                if (faces[cell_faces[j]].boundary < 0){
                    if (i == faces[cell_faces[j]].ownerCell){N = faces[cell_faces[j]].neighbourCell;}
                    else {N = faces[cell_faces[j]].ownerCell;}

                    d_N = cells[N].centroid - cells[i].centroid;

                    //If face is not at boundary, we can take ø_N directly from neighbour cell centroid
                    ø_N = cells[N].T;
                }

                //If face is at boundary, d_N is distance from cell centroid to face centroid
                else{
                    d_N = faces[cell_faces[j]].centroid - cells[i].centroid;

                    //we need to check face boundary condition here to find ø_N = ø_b
                    //If face boundary condition is dirichlet we can take value directly
                    boundary = faces[cell_faces[j]].boundary;
                    if (bc_type[boundary] == 0){ø_N = boundary_values(0, boundary);}

                    //If face boundary condition is Neumann g_b, we interpolate face value as
                    //ø_N = ø_b = ø_P + |d_N|*(n dot g_b) (where n denotes face normal)
                    else{
                        face_normal = faces[cell_faces[j]].Sf/faces[cell_faces[j]].Sf.norm();
                        ø_N = ø_i + (d_N.norm()*(face_normal.dot(boundary_values.col(boundary))));
                    }
                }

                gradient += (pow(d_N.norm(),2) * (ø_N - ø_i) * d_N);
            }
        }

        //Now calculate final gradient by left-multiplying with G_inv = inv(M^T M)
        gradient = G_inv*gradient;

        //We suppress excessively small gradient components
        for (int k =0; k < 3; k++){
            if (std::abs(gradient(k)) < 1e-10){
                gradient(k) = 0;
            }
        }

        cells[i].grad = gradient;
    }
}

