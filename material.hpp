#pragma once

#include "common.hpp"

class Material
{
    public:
    virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Vector3d& strain) const = 0;

};

class PlaneStrainMaterial : public Material
{
    public:
    explicit PlaneStrainMaterial(double E, double mu)
        : _E(E), _mu(mu)
    {}

    virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Vector3d&  strain ) const override
    {
        // elasticity matrix for plane strain
        Eigen::Matrix3d D_mat;
        D_mat << (1-_mu), _mu, 0,
                    _mu, (1-_mu), 0,
                    0, 0, 0.5*(1-2*_mu);
        D_mat *= _E / ( (1+_mu) * (1-2*_mu));

        Eigen::Vector3d stress_vec;
        stress_vec = D_mat * strain;
        return std::make_pair(stress_vec, D_mat);
    }

    private:
    double _E; // Young's modulus
    double _mu; // Poisson's ratio
};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

class HW3Material : public Material
{
    public:
    explicit HW3Material(double a, double b, double c)
        : _a(a), _b(b), _c(c)
    {}

    virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Vector3d& strain) const override
    {
        const double ex = strain[0];
        const double ey = strain[1];
        const double exy = 0.5*strain[2];

        // elasticity matrix
        Eigen::Matrix3d D_mat;

        D_mat(0,0) = 2*_a + 6*_b*(ex + ey) + 6*_c*ex;
        D_mat(0,1) = 6*_b*(ex + ey);
        D_mat(0,2) = 3*_c*(exy);

        D_mat(1,0) = 6*_b*(ex + ey);
        D_mat(1,1) = 2*_a + 6*_b*(ex + ey) + 6*_c*ey;
        D_mat(1,2) = 3*_c*exy;

        D_mat(2,0) = 3*_c*exy;
        D_mat(2,1) = 3*_c*exy;
        D_mat(2,2) = _a + 1.5*_c*(ex+ey);

        // stress (from stress-strain relationship)
        Eigen::Vector3d stress_vec;
        stress_vec[0] = 2*_a*ex + 3*_b*(ex + ey)*(ex + ey) + 3*_c*(ex*ex + exy*exy);
        stress_vec[1] = 2*_a*ey + 3*_b*(ex + ey)*(ex + ey) + 3*_c*(ey*ey + exy*exy);
        stress_vec[2] = 2*_a*exy + 3*_c*exy*(ex + ey);

        return std::make_pair(stress_vec, D_mat);
    }

    private:
    double _a;
    double _b;
    double _c;
};