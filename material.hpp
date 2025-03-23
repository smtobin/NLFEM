#pragma once

#include "common.hpp"
#include <iostream>

class Material
{
    public:
    virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Matrix2d& F) const = 0;

};

// class PlaneStrainMaterial : public Material
// {
//     public:
//     explicit PlaneStrainMaterial(double E, double mu)
//         : _E(E), _mu(mu)
//     {}

//     virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Vector3d&  strain ) const override
//     {
//         // elasticity matrix for plane strain
//         Eigen::Matrix3d D_mat;
//         D_mat << (1-_mu), _mu, 0,
//                     _mu, (1-_mu), 0,
//                     0, 0, 0.5*(1-2*_mu);
//         D_mat *= _E / ( (1+_mu) * (1-2*_mu));

//         Eigen::Vector3d stress_vec;
//         stress_vec = D_mat * strain;
//         return std::make_pair(stress_vec, D_mat);
//     }

//     private:
//     double _E; // Young's modulus
//     double _mu; // Poisson's ratio
// };

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

// class HW3Material : public Material
// {
//     public:
//     explicit HW3Material(double a, double b, double c)
//         : _a(a), _b(b), _c(c)
//     {}

//     virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Vector3d& strain) const override
//     {
//         const double ex = strain[0];
//         const double ey = strain[1];
//         const double exy = 0.5*strain[2];

//         // elasticity matrix
//         Eigen::Matrix3d D_mat;

//         D_mat(0,0) = 2*_a + 6*_b*(ex + ey) + 6*_c*ex;
//         D_mat(0,1) = 6*_b*(ex + ey);
//         D_mat(0,2) = 3*_c*(exy);

//         D_mat(1,0) = 6*_b*(ex + ey);
//         D_mat(1,1) = 2*_a + 6*_b*(ex + ey) + 6*_c*ey;
//         D_mat(1,2) = 3*_c*exy;

//         D_mat(2,0) = 3*_c*exy;
//         D_mat(2,1) = 3*_c*exy;
//         D_mat(2,2) = _a + 1.5*_c*(ex+ey);

//         // stress (from stress-strain relationship)
//         Eigen::Vector3d stress_vec;
//         stress_vec[0] = 2*_a*ex + 3*_b*(ex + ey)*(ex + ey) + 3*_c*(ex*ex + exy*exy);
//         stress_vec[1] = 2*_a*ey + 3*_b*(ex + ey)*(ex + ey) + 3*_c*(ey*ey + exy*exy);
//         stress_vec[2] = 2*_a*exy + 3*_c*exy*(ex + ey);

//         return std::make_pair(stress_vec, D_mat);
//     }

//     private:
//     double _a;
//     double _b;
//     double _c;
// };

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class MidtermMaterial : public Material
{
    public:
    explicit MidtermMaterial(double E, double nu)
    {
        _lambda = (E*nu) / ( (1+nu) * (1-2*nu) );
        _mu = E / ( 2*(1+nu) );
    }

    virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Matrix2d& F) const override
    {
        // calculate Green strain
        // const Eigen::Matrix2d E = 0.5*(F.transpose()*F - Eigen::Matrix2d::Identity());

        // calculate Jacobian J = det(F)
        const double J = F.determinant();

        // calculate 2nd PK stress
        const Eigen::Matrix2d inv_C = (F.transpose()*F).inverse();
        const Eigen::Matrix2d S = _mu * (Eigen::Matrix2d::Identity() - inv_C) + _lambda * (J - 1) * J * inv_C;

        // calculate Cauchy stress from 2nd PK stress
        const Eigen::Matrix2d sigma = (1/J) * F * S * F.transpose();

        // elasticity matrix
        Eigen::Matrix3d D_mat;

        // elasticity matrix - from c_ijkl
        D_mat(0,0) = 2*_mu*(1/J) + _lambda;
        D_mat(0,1) = _lambda*(2*J - 1);
        D_mat(0,2) = 0;

        D_mat(1,0) = _lambda*(2*J-1);
        D_mat(1,1) = 2*_mu*(1/J) + _lambda;
        D_mat(1,2) = 0;

        D_mat(2,0) = 0;
        D_mat(2,1) = 0;
        D_mat(2,2) = _mu*(1/J) - _lambda*(J - 1);

        // stress as a vector
        Eigen::Vector3d stress_vec;
        stress_vec[0] = sigma(0,0);
        stress_vec[1] = sigma(1,1);
        stress_vec[2] = sigma(0,1);

        return std::make_pair(stress_vec, D_mat);
    }

    private:
    double _lambda;
    double _mu;
};