#include <Eigen/Dense>
#include <iostream>

#include "element.hpp"

int main(int argc, char **argv) 
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

    std::cout << element.M() << std::endl;
    std::cout << element.K() << std::endl;

}