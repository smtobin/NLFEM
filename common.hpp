#pragma once

#define NSDIMS 2    // number of spatial dimensions
#define NR_TOL 1e-20 // termination tolerance for Newton-Raphson
#define NR_MAX_ITER 20 // max number of iterations for Newton-Raphson loop

#include <Eigen/Dense>

typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Vector<double, 6> Vector6d;
typedef Eigen::Vector<double, NSDIMS> Node;
typedef Eigen::Vector<int, 4> ElementNodes;


enum Axis
{
    X = 0,
    Y = 1,
    Z = 2
};

struct PlasticState
{
    Vector6d elastic_strain;    // elastic strain
    // Vector6d plastic_strain;    // plastic strain
    double alpha; // isotropic hardening
    Vector6d beta; // back stress
};