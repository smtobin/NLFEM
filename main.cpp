#include <Eigen/Dense>
#include <iostream>

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

// define the vertices of the quadrilateral 2D element
const Eigen::Vector2d x1(0.5,0.5);
const Eigen::Vector2d x2(-0.5,0.5);
const Eigen::Vector2d x3(-0.5,-0.5);
const Eigen::Vector2d x4(0.5,-0.5);

// Eigen::Matrix<double,8,8> stiffnessMatrix();
// Eigen::Vector<double,8> internalForce(const Eigen::Vector<double,8>& d);
// Eigen::Vector<double,8> newtonRaphson(const Eigen::Vector<double,8>& F_ext, const Eigen::Vector<double,8>& d0);
void applyNodalForce(int ind, double force, Eigen::Vector<double,8>& R);
void applyDisplacementBC(int ind, double displacement, Eigen::Matrix<double,8,8>& K, Eigen::Vector<double,8>& R);



int main(int /* argc */, char** /* argv */) 
{   
    // create nodes of mesh
    std::vector<Node> nodes;
    nodes.emplace_back(0.5, 0.5);
    nodes.emplace_back(-0.5, 0.5);
    nodes.emplace_back(-0.5, -0.5);
    nodes.emplace_back(0.5, -0.5);

    // create list of element's nodes
    std::vector<ElementNodes> element_nodes;
    element_nodes.emplace_back(0, 1, 2, 3);

    // create displacement BCs
    std::vector<DisplacementBC> displacement_BCs;
    displacement_BCs.emplace_back(3, Axis::Y, 0);
    displacement_BCs.emplace_back(0, Axis::X, 0.01);
    displacement_BCs.emplace_back(1, Axis::X, 0);
    displacement_BCs.emplace_back(2, Axis::X, 0);
    displacement_BCs.emplace_back(3, Axis::X, 0.01);

    // create force BCs
    std::vector<ForceBC> force_BCs;

    Solver solver(nodes, element_nodes, displacement_BCs, force_BCs, density, E, mu);
    solver.solve();
    solver.evaluateElementAtIntegrationPoints(0);
}

void applyNodalForce(int ind, double force, Eigen::Vector<double,8>& R)
{
    R(ind) += force;
}

void applyDisplacementBC(int ind, double displacement, Eigen::Matrix<double,8,8>& K, Eigen::Vector<double,8>& R)
{
    // modify the stiffness matrix and load vector
    K.row(ind) = Eigen::Vector<double,8>::Zero();
    K(ind,ind) = 1;
    for (int i = 0; i < 8; i++)
    {
        R(i) -= K(i,ind)*displacement;
    }
    R(ind) = displacement;
    K.col(ind) = K.row(ind);
}

// Eigen::Matrix<double,8,8> stiffnessMatrix()
// {
//     Eigen::Matrix<double,8,8> element_K = element.K();
    
//     // manually apply displacement BCs
//     //  - all rows corresponding to unknown displacements first
//     //  - rows corresponding to known displacements
//     Eigen::Matrix<double,8,8> K;
//     K.row(0) = element_K.row(1);
//     K.row(1) = element_K.row(3);
//     K.row(2) = element_K.row(5);
//     K.row(3) = element_K.row(7);
//     K.row(4) = element_K.row(0);
//     K.row(5) = element_K.row(2);
//     K.row(6) = element_K.row(4);
//     K.row(7) = element_K.row(6);

//     return K;
// }

// Eigen::Vector<double,8> internalForce(const Eigen::Vector<double,8>& d)
// {
//     Eigen::Matrix<double, 8, 8> K = element.K();
//     return K*d;
// }

// Eigen::Vector<double,8> newtonRaphson(const Eigen::Vector<double,8>& F_ext, const Eigen::Vector<double,8>& d0)
// {
//     Eigen::Vector<double,8> R = F_ext - internalForce(d0);
//     Eigen::Vector<double,8> d = d0;
    
//     int max_iter = 1000;
//     double tol = 1e-12;
//     for (int i = 0; i < max_iter; i++)
//     {
//         if (R.norm() < R.norm() * tol)
//             break;

//         Eigen::Matrix<double,8,8> K = stiffnessMatrix();
//         Eigen::Vector<double,8> delta_d = K.inverse() * R;
//         d += delta_d;

//         R = F_ext - internalForce(d0);
//     }

//     return d;
// }