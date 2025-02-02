#include <Eigen/Dense>
#include <iostream>

#include "element.hpp"

void applyNodalForce(int ind, double force, Eigen::Vector<double,8>& R);
void applyDisplacementBC(int ind, double displacement, Eigen::Matrix<double,8,8>& K, Eigen::Vector<double,8>& R);

int main(int /* argc */, char** /* argv */) 
{
    // define the material (given in the problem statement)
    const double E = 100;
    const double mu = 0.25;
    const double density = 1; 

    // define the vertices of the quadrilateral 2D element
    const Eigen::Vector2d x1(0.5,0.5);
    const Eigen::Vector2d x2(-0.5,0.5);
    const Eigen::Vector2d x3(-0.5,-0.5);
    const Eigen::Vector2d x4(0.5,-0.5);
    
    // create the element
    QuadElement element(x1, x2, x3, x4, density, E, mu);

    Eigen::Matrix<double, 8, 8> K = element.K();
    Eigen::Vector<double, 8> R({0, 0, 0, 0, 0, 0, 0, 0});

    // prescribe displacements
    //  u1 = 0.01
    //  u2 = 0
    //  u3 = 0
    //  u4 = 0.01
    applyDisplacementBC(0, 0.01, K, R); // u1
    applyDisplacementBC(2, 0, K, R);    // u2
    applyDisplacementBC(4, 0, K, R);    // u3
    applyDisplacementBC(6, 0.01, K, R); // u4

    // std::cout << "Modified K:\n" << K << std::endl;
    // std::cout << "Modified R:\n" << R << std::endl;

    // solve KU = R using Eigen's LLT solver
    Eigen::Vector<double, 8> U = K.llt().solve(R);

    // evaluate stress and strain at integration points
    const std::vector<double>& integration_points = element.integrationPoints();
    for (const auto& ri : integration_points)
    {
        for (const auto& sj : integration_points)
        {
            std::cout << "\n === At integration point (r=" << ri << ", s=" << sj << ") ===" << std::endl;
            Eigen::Matrix<double, 3, 8> B = element.B(ri, sj);
            Eigen::Vector3d strain_vec = B*U;
            Eigen::Vector3d stress_vec = element.C()*B*U;
            std::cout << "  (e_xx, e_yy, e_xy):       " << strain_vec[0] << ", " << strain_vec[1] << ", " << strain_vec[2] << std::endl;
            std::cout << "  (sig_xx, sig_yy, sig_xy): " << stress_vec[0] << ", " << stress_vec[1] << ", " << stress_vec[2] << std::endl;
        }
    }
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