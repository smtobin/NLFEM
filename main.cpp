#include <Eigen/Dense>
#include <iostream>

#include "element.hpp"

int main(int /* argc */, char** /* argv */) 
{
    // define the material (given in the problem statement)
    const double E = 100;
    const double mu = 0.25; 

    const Eigen::Vector2d x1(0.5,0.5);
    const Eigen::Vector2d x2(-0.5,0.5);
    const Eigen::Vector2d x3(-0.5,-0.5);
    const Eigen::Vector2d x4(0.5,-0.5);
    const double density = 1;
    QuadElement element(x1, x2, x3, x4, density, E, mu);

    Eigen::Matrix<double, 8, 8> K = element.K();
    Eigen::Vector<double, 8> R({0, 0, 0, 0, 0, 0, 0, 0});

    // prescribe displacements
    //  u1 = 0.01
    //  u2 = 0
    //  u3 = 0
    //  u4 = 0.01

    // u1
    K.row(0) = Eigen::Vector<double,8>({1,0,0,0,0,0,0,0});
    for (int i = 0; i < 8; i++)
    {
        R(i) -= K(i,0)*0.01;
    }
    R(0) = 0.01;
    K.col(0) = K.row(0);

    // u2
    K.row(2) = Eigen::Vector<double,8>({0,0,1,0,0,0,0,0});
    R(2) = 0;
    K.col(2) = K.row(2);

    // u3
    K.row(4) = Eigen::Vector<double,8>({0,0,0,0,1,0,0,0});
    R(4) = 0;
    K.col(4) = K.row(4);

    // u4
    K.row(6) = Eigen::Vector<double,8>({0,0,0,0,0,0,1,0});
    R(1) -= K(1,6)*0.01;
    R(3) -= K(3,6)*0.01;
    R(5) -= K(5,6)*0.01;
    R(7) -= K(7,6)*0.01;
    R(6) = 0.01;
    K.col(6) = K.row(6);

    std::cout << "Modified K:\n" << K << std::endl;
    std::cout << "Modified R:\n" << R << std::endl;

    Eigen::Vector<double, 8> U = K.llt().solve(R);
    std::cout << "\nK:\n" << element.K() << std::endl;
    std::cout << "\nM:\n" << element.M() << std::endl;
    std::cout << "\nU:\n" << U << std::endl;

    //////

    Eigen::Matrix<double, 3, 8> B = element.B(-0.577, -0.577);
    // Eigen::Vector<double, 8> u_hat({0.5, 0, 0, 0, 0, 0, 0, 0});
    std::cout << "\n" << B*U << std::endl;
    std::cout << "\n" << element.C()*B*U << std::endl;
}