#pragma once

#include "common.hpp"
#include "material.hpp"

#include <iostream>
#include <cassert>

// 2D Quadrilateral element with plane strain
class QuadElement
{
    public:
    QuadElement()
    {}

    QuadElement(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2, const Eigen::Vector2d& x3, const Eigen::Vector2d& x4, const Material* material)
        : _x1(x1), _x2(x2), _x3(x3), _x4(x4),
          _material(material),
          _integration_points({-1.0/std::sqrt(3), 1.0/std::sqrt(3)}),
          _integration_weights({1.0, 1.0}),
          _global_DOF()
    {
    }

    int numNodes() const { return 4; }

    void setGlobalDOF(const std::vector<int>& global_dof)
    {
        assert(global_dof.size() == 8);
        _global_DOF = global_dof;
    }

    const std::vector<int>& globalDOF() const { return _global_DOF; }

    // evaluates the interpolation matrix H at (r,s)
    Eigen::Matrix<double, 2, 8> H(double r, double s) const;

    // evaluates the strain-displacement matrix B at (r,s)
    Eigen::Matrix<double, 3, 8> B(double r, double s) const;

    // evaluates the Jacobian operator matrix at (r,s)
    Eigen::Matrix2d J(double r, double s) const;

    // evalutes the mass matrix M for the element
    // Eigen::Matrix<double, 8, 8> M() const;

    // evalutes the stiffness matrix K for the element for given nodal displacements
    Eigen::Matrix<double, 8, 8> K(const Eigen::VectorXd& d_e) const;

    // evaulates the internal force vector for the element for given nodal displacements
    Eigen::Vector<double, 8> internalForce(const Eigen::VectorXd& d_e) const;

    // evalutes the deformation gradient at (r,s) for given nodal displacements
    Eigen::Matrix2d deformationGradient(double r, double s, const Eigen::VectorXd& d_e) const;

    const Material* material() const { return _material; }

    const std::vector<double>& integrationPoints() const
    {
        return _integration_points;
    }

    const std::vector<double>& integrationWeights() const
    {
        return _integration_weights;
    }

    private:
    Eigen::Vector2d _x1;
    Eigen::Vector2d _x2;
    Eigen::Vector2d _x3;
    Eigen::Vector2d _x4;

    const Material* _material;

    std::vector<double> _integration_points;
    std::vector<double> _integration_weights;

    std::vector<int> _global_DOF;
};