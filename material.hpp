#pragma once

#include "common.hpp"
#include <iostream>

class Material
{
    public:
    virtual std::tuple<Vector6d, Matrix6d, PlasticState> materialSubroutine(const Vector6d& strain, const PlasticState& last_state) const = 0;

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

//     virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Vector6d& strain, const PlasticState& last_state) const override
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

// class MidtermMaterial : public Material
// {
//     public:
//     explicit MidtermMaterial(double E, double nu)
//     {
//         _lambda = (E*nu) / ( (1+nu) * (1-2*nu) );
//         _mu = E / ( 2*(1+nu) );
//     }

//     virtual std::pair<Eigen::Vector3d, Eigen::Matrix3d> materialSubroutine(const Eigen::Matrix2d& F) const override
//     {
//         // calculate Green strain
//         // const Eigen::Matrix2d E = 0.5*(F.transpose()*F - Eigen::Matrix2d::Identity());

//         // calculate Jacobian J = det(F)
//         const double J = F.determinant();

//         // calculate 2nd PK stress
//         const Eigen::Matrix2d inv_C = (F.transpose()*F).inverse();
//         const Eigen::Matrix2d S = _mu * (Eigen::Matrix2d::Identity() - inv_C) + _lambda * (J - 1) * J * inv_C;

//         // calculate Cauchy stress from 2nd PK stress
//         const Eigen::Matrix2d sigma = (1/J) * F * S * F.transpose();

//         // elasticity matrix
//         Eigen::Matrix3d D_mat;

//         // elasticity matrix - from c_ijkl
//         D_mat(0,0) = 2*_mu*(1/J) + _lambda;
//         D_mat(0,1) = _lambda*(2*J - 1);
//         D_mat(0,2) = 0;

//         D_mat(1,0) = _lambda*(2*J-1);
//         D_mat(1,1) = 2*_mu*(1/J) + _lambda;
//         D_mat(1,2) = 0;

//         D_mat(2,0) = 0;
//         D_mat(2,1) = 0;
//         D_mat(2,2) = _mu*(1/J) - _lambda*(J - 1);

//         // stress as a vector
//         Eigen::Vector3d stress_vec;
//         stress_vec[0] = sigma(0,0);
//         stress_vec[1] = sigma(1,1);
//         stress_vec[2] = sigma(0,1);

//         return std::make_pair(stress_vec, D_mat);
//     }

//     private:
//     double _lambda;
//     double _mu;
// };


class HW6Material : public Material
{
    public:
    explicit HW6Material(double E, double nu, double sigma_y, double beta, double H_bar_prime, double dt)
        : _E(E), _nu(nu), _sigma_y(sigma_y), _beta(beta), _H_bar_prime(H_bar_prime), _dt(dt)
    {}

    virtual std::tuple<Vector6d, Matrix6d, PlasticState> materialSubroutine(const Vector6d& strain, const PlasticState& last_state) const override
    {
        return radialReturn(strain, last_state);
    }

    /** Given new strain and the last plastic state, returns the
     * 1. 6x1 new stress
     * 2. 6x6 consistent tangent moduli
     * 3. the new plastic state
     */
    std::tuple<Vector6d, Matrix6d, PlasticState> radialReturn(const Vector6d& new_strain, const PlasticState& last_state) const
    {
        const double G = _E / (2 * (1+_nu));  // shear modulus
        const double K = _E / (3 * (1-2*_nu));    // bulk modulus

        // 1. Compute trial elastic stress
        Eigen::Matrix3d strain_mat_new = strainVoigt2Mat(new_strain);
        Eigen::Matrix3d e_p_old = strainVoigt2Mat(last_state.dev_plastic_strain);
        Eigen::Matrix3d beta_old = stressVoigt2Mat(last_state.beta);
        double alpha_old = last_state.alpha;

        Eigen::Matrix3d e_new = strain_mat_new - 1.0/3.0 * strain_mat_new.trace() * Eigen::Matrix3d::Identity();   // deviatoric part of new strain
        Eigen::Matrix3d s_new_trial = 2*G*(e_new - e_p_old);
        Eigen::Matrix3d xi_new_trial = s_new_trial - beta_old;

        // 2. Check yield condition
        double f_new_trial = xi_new_trial.norm() - std::sqrt(2.0/3.0) * _isotropicHardeningModulus(alpha_old);
        // std::cout << "s_trial_11: " << s_new_trial << "   ";
        // std::cout << "f_trial: " << f_new_trial << "   ";

        PlasticState new_state;
        // we are not plastically deforming
        if (f_new_trial <= 0)
        {
            new_state = last_state;

            Vector6d stress = stressMat2Voigt(K * strain_mat_new.trace() * Eigen::Matrix3d::Identity() + s_new_trial);

            Matrix6d D_ep;
            D_ep << 1-_nu, _nu, _nu, 0, 0, 0,
                    _nu, 1-_nu, _nu, 0, 0, 0,
                    _nu, _nu, 1-_nu, 0, 0, 0,
                    0, 0, 0, 0.5*(1-2*_nu), 0, 0,
                    0, 0, 0, 0, 0.5*(1-2*_nu), 0,
                    0, 0, 0, 0, 0, 0.5*(1-2*_nu);
            D_ep *= _E / ( (1+_nu) * (1-2*_nu) );
            return std::make_tuple(stress, D_ep, new_state);
        }

        // std::cout << "PLASTIC";

        // plastic deformation, project elastic state onto plastic yield surface
        // 3. Compute n_n+1 and delta_gamma
        Eigen::Matrix3d n_new = xi_new_trial / xi_new_trial.norm();
        double dgamma = f_new_trial / ( 2*G * ( 1 + _H_bar_prime / (3*G)) );

        // 4. Update quantities
        double alpha_new = alpha_old + std::sqrt(2.0/3.0) * dgamma;
        Eigen::Matrix3d beta_new = beta_old + std::sqrt(2.0/3.0) * (_kinematicHardeningModulus(alpha_new) - _kinematicHardeningModulus(alpha_old)) * n_new;
        Eigen::Matrix3d e_p_new = e_p_old + dgamma * n_new;
        Eigen::Matrix3d sigma_new = K * strain_mat_new.trace() * Eigen::Matrix3d::Identity() + s_new_trial - 2*G*dgamma*n_new;

        // 5. Compute consistent elastoplastic tangent moduli
        double theta = 1 - (2*G*dgamma) / xi_new_trial.norm();
        double theta_bar = 1 / (1 + (_isotropicHardeningModulusPrime(dgamma) + _kinematicHardeningModulusPrime(dgamma)) / (3*G) ) - (1 - theta);
        
        Matrix6d D_ep;
        D_ep(0,0) = K + 4.0/3.0*G*theta - 2*G*theta_bar*n_new(0,0)*n_new(0,0);
        D_ep(0,1) = K - 2.0/3.0*G*theta - 2*G*theta_bar*n_new(0,0)*n_new(1,1);
        D_ep(0,2) = K - 2.0/3.0*G*theta - 2*G*theta_bar*n_new(0,0)*n_new(2,2);
        D_ep(0,3) = -2*G*theta_bar*n_new(0,0)*n_new(1,2);
        D_ep(0,4) = -2*G*theta_bar*n_new(0,0)*n_new(0,2);
        D_ep(0,5) = -2*G*theta_bar*n_new(0,0)*n_new(0,1);

        D_ep(1,0) = D_ep(0,1);
        D_ep(1,1) = K + 4.0/3.0*G*theta - 2*G*theta_bar*n_new(1,1)*n_new(1,1);
        D_ep(1,2) = K - 2.0/3.0*G*theta - 2*G*theta_bar*n_new(1,1)*n_new(2,2);
        D_ep(1,3) = -2*G*theta_bar*n_new(1,1)*n_new(1,2);
        D_ep(1,4) = -2*G*theta_bar*n_new(1,1)*n_new(0,2);
        D_ep(1,5) = -2*G*theta_bar*n_new(1,1)*n_new(0,1);

        D_ep(2,0) = D_ep(0,2);
        D_ep(2,1) = D_ep(1,2);
        D_ep(2,2) = K + 4.0/3.0*G*theta - 2*G*theta_bar*n_new(2,2)*n_new(2,2);
        D_ep(2,3) = -2*G*theta_bar*n_new(2,2)*n_new(1,2);
        D_ep(2,4) = -2*G*theta_bar*n_new(2,2)*n_new(0,2);
        D_ep(2,5) = -2*G*theta_bar*n_new(2,2)*n_new(0,1);

        D_ep(3,0) = D_ep(0,3);
        D_ep(3,1) = D_ep(1,3);
        D_ep(3,2) = D_ep(2,3);
        D_ep(3,3) = G*theta - 2*G*theta_bar*n_new(1,2)*n_new(1,2);
        D_ep(3,4) = -2*G*theta_bar*n_new(1,2)*n_new(0,2);
        D_ep(3,5) = -2*G*theta_bar*n_new(1,2)*n_new(0,1);

        D_ep(4,0) = D_ep(0,4);
        D_ep(4,1) = D_ep(1,4);
        D_ep(4,2) = D_ep(2,4);
        D_ep(4,3) = D_ep(3,4);
        D_ep(4,4) = G*theta - 2*G*theta_bar*n_new(0,2)*n_new(0,2);
        D_ep(4,5) = -2*G*theta_bar*n_new(0,2)*n_new(0,1);

        D_ep(5,0) = D_ep(0,5);
        D_ep(5,1) = D_ep(1,5);
        D_ep(5,2) = D_ep(2,5);
        D_ep(5,3) = D_ep(3,5);
        D_ep(5,4) = D_ep(4,5);
        D_ep(5,5) = G*theta - 2*G*theta_bar*n_new(0,1)*n_new(0,1);

        // assemble new state struct
        new_state.alpha = alpha_new;
        new_state.beta = stressMat2Voigt(beta_new);
        Vector6d stress = stressMat2Voigt(sigma_new);
        new_state.dev_plastic_strain = strainMat2Voigt(e_p_new);
        return std::make_tuple(stress, D_ep, new_state);
    }

    Eigen::Matrix3d strainVoigt2Mat(const Vector6d& strain_vec) const
    {
        Eigen::Matrix3d strain_mat;
        strain_mat <<       strain_vec[0], 0.5*strain_vec[5], 0.5*strain_vec[4],
                            0.5*strain_vec[5], strain_vec[1], 0.5*strain_vec[3],
                            0.5*strain_vec[4], 0.5*strain_vec[3], strain_vec[2];

        return strain_mat;
    }

    Vector6d strainMat2Voigt(const Eigen::Matrix3d& strain_mat) const
    {
        Vector6d strain_vec;
        strain_vec[0] = strain_mat(0,0);
        strain_vec[1] = strain_mat(1,1);
        strain_vec[2] = strain_mat(2,2);
        strain_vec[3] = 2*strain_mat(1,2);
        strain_vec[4] = 2*strain_mat(0,2);
        strain_vec[5] = 2*strain_mat(0,1);

        return strain_vec;
    }

    Eigen::Matrix3d stressVoigt2Mat(const Vector6d& stress_vec) const
    {
        Eigen::Matrix3d stress_mat;
        stress_mat <<       stress_vec[0], stress_vec[5], stress_vec[4],
                            stress_vec[5], stress_vec[1], stress_vec[3],
                            stress_vec[4], stress_vec[3], stress_vec[2];
        return stress_mat;
    }

    Vector6d stressMat2Voigt(const Eigen::Matrix3d& stress_mat) const
    {
        Vector6d stress_vec;
        stress_vec[0] = stress_mat(0,0);
        stress_vec[1] = stress_mat(1,1);
        stress_vec[2] = stress_mat(2,2);
        stress_vec[3] = stress_mat(1,2);
        stress_vec[4] = stress_mat(0,2);
        stress_vec[5] = stress_mat(0,1);

        return stress_vec;
    }

    private:
    double _isotropicHardeningModulus(double alpha) const
    {
        return _sigma_y + _beta * _H_bar_prime * alpha;
    }

    double _kinematicHardeningModulus(double alpha) const
    {
        return (1-_beta) * _H_bar_prime * alpha;
    }

    double _isotropicHardeningModulusPrime(double dgamma) const
    {
        return _beta * _H_bar_prime * dgamma / _dt * std::sqrt(2.0/3.0);
    }

    double _kinematicHardeningModulusPrime(double dgamma) const
    {
        return (1-_beta) * _H_bar_prime * dgamma / _dt * std::sqrt(2.0/3.0);
    }

    double _E;
    double _nu;
    double _sigma_y;
    double _beta;
    double _H_bar_prime;
    double _dt;
};