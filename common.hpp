#pragma once

#define NSDIMS 2    // number of spatial dimensions

#include <Eigen/Dense>

typedef Eigen::Vector<double, NSDIMS> Node;
typedef Eigen::Vector<int, 4> ElementNodes;


enum Axis
{
    X = 0,
    Y = 1,
    Z = 2
};