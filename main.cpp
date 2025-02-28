#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>

#include "common.hpp"
#include "element.hpp"
#include "solver.hpp"

#define NSDIMS 2 // number of spatial dimensions

/** TODO LIST:
 * 1. Create assembly operator routine, mapping element DOF to global DOF (even for just one element) - DONE
 * 2. Assign global DOF based on displacement BCs (list the unknown displacements first, then the known) - DONE
 * 3. Test Linear FEM assignment with new assembly operator- DONE
 * 4. Implement Newton-Raphson loop
 * 5. Test Linear FEM assignment with Newton-Raphson loop and two displacement steps
 * 6. Implement material model and nonlinear stress-strain relation
 * 7. Make necessary changes to Newton-Raphson loop to accomodate new stress-strain relation and test
 */



// define the material (given in the problem statement)
const double E = 100;
const double mu = 0.25;
const double density = 1; 

// Eigen::Matrix<double,8,8> stiffnessMatrix();
// Eigen::Vector<double,8> internalForce(const Eigen::Vector<double,8>& d);
// Eigen::Vector<double,8> newtonRaphson(const Eigen::Vector<double,8>& F_ext, const Eigen::Vector<double,8>& d0);

void readInputData(const std::string& filename,
                   std::vector<Node>& nodes, std::vector<ElementNodes>& element_nodes,
                   std::vector<DisplacementBC>& displacement_BCs, std::vector<ForceBC>& force_BCs,
                   std::unique_ptr<Material>& material)
{
    std::ifstream infile(filename);
    if (!infile.good())
    {
        std::cerr << "ERROR: Problem reading input file. Check the file path. Aborting..." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        char c;
        if (line[0] == 'n')
        {
            int i;
            double x, y;
            if (!(iss >> c >> i >> x >> y)) { assert(0); }
            nodes.emplace_back(x, y);
        }
        else if (line[0] == 'e')
        {
            int n1, n2, n3, n4;
            if (!(iss >> c >> n1 >> n2 >> n3 >> n4)) { assert(0); }
            element_nodes.emplace_back(n1, n2, n3, n4);
        }
        else if (line[0] == 'd')
        {
            int i;
            char ax;
            double d;
            if (!(iss >> c >> i >> ax >> d)) { assert(0); }
            Axis axis;
            if (ax == 'x') { axis = Axis::X; }
            else if (ax == 'y') { axis = Axis::Y; }
            else { assert(0); }

            displacement_BCs.emplace_back(i, axis, d);
        }
        else if (line[0] == 'f')
        {
            int i;
            char ax;
            double d;
            if (!(iss >> c >> i >> ax >> d)) { assert(0); }
            Axis axis;
            if (ax == 'x') { axis = Axis::X; }
            else if (ax == 'y') { axis = Axis::Y; }
            else { assert(0); }

            force_BCs.emplace_back(i, axis, d);
        }
        else if (line[0] == 'm')
        {
            std::string material_str;
            if (!(iss >> c >> material_str)) { assert(0); }

            if (material_str == "PlaneStrain")
            {
                material = std::make_unique<PlaneStrainMaterial>(100, 0.25); // parameters from linear FEM assignment
            }
            else if (material_str == "HW3")
            {
                material = std::make_unique<HW3Material>(40, -50, -30); // parameters from HW3 assignment
            }
            else
            {
                assert(0);
            }
        }
    }
}

int main(int argc, char** argv) 
{   

    assert(argc == 2);
    std::string input_filename(argv[1]);

    // create nodes of mesh
    std::vector<Node> nodes;
    std::vector<ElementNodes> element_nodes;
    std::vector<DisplacementBC> displacement_BCs;
    std::vector<ForceBC> force_BCs;
    std::unique_ptr<Material> material;

    readInputData(input_filename, nodes, element_nodes, displacement_BCs, force_BCs, material);

    std::cout << "Read in " << nodes.size() << " nodes, " << element_nodes.size() << " elements, " <<
        displacement_BCs.size() << " displacement BCs, and " << force_BCs.size() << " force BCs from input file." << std::endl;

    Solver solver(nodes, element_nodes, displacement_BCs, force_BCs, material.get());
    solver.solve(10);
    solver.evaluateElementAtIntegrationPoints(0);
}