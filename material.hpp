#pragma once

#include "common.hpp"
#include <iostream>

class Material
{
    public:
    virtual std::tuple<Vector6d, Matrix6d, PlasticState> materialSubroutine(const Vector6d& elastic_trial_strain, const Eigen::Matrix3d& betr_n1, const PlasticState& last_state) const= 0;

    static Eigen::Matrix3d strainVoigt2Mat(const Vector6d& strain_vec)
    {
        Eigen::Matrix3d strain_mat;
        strain_mat <<       strain_vec[0], 0.5*strain_vec[5], 0.5*strain_vec[4],
                            0.5*strain_vec[5], strain_vec[1], 0.5*strain_vec[3],
                            0.5*strain_vec[4], 0.5*strain_vec[3], strain_vec[2];

        return strain_mat;
    }

    static Vector6d strainMat2Voigt(const Eigen::Matrix3d& strain_mat)
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

    static Eigen::Matrix3d stressVoigt2Mat(const Vector6d& stress_vec)
    {
        Eigen::Matrix3d stress_mat;
        stress_mat <<       stress_vec[0], stress_vec[5], stress_vec[4],
                            stress_vec[5], stress_vec[1], stress_vec[3],
                            stress_vec[4], stress_vec[3], stress_vec[2];
        return stress_mat;
    }

    static Vector6d stressMat2Voigt(const Eigen::Matrix3d& stress_mat)
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

class FinalMaterial : public Material
{
    public:
    explicit FinalMaterial(double E, double nu, double sigma_y, double beta, double H_bar_prime)
    : _E(E), _nu(nu), _sigma_y(sigma_y), _beta(beta), _H_bar_prime(H_bar_prime)
    {}

    virtual std::tuple<Vector6d, Matrix6d, PlasticState> materialSubroutine(const Vector6d& elastic_trial_strain, const Eigen::Matrix3d& betr_n1, const PlasticState& last_state) const override
    {
        return radialReturn(elastic_trial_strain, betr_n1, last_state);
    }

    /** Given new strain and the last plastic state, returns the
     * 1. 6x1 new stress
     * 2. 6x6 consistent tangent moduli
     * 3. the new plastic state
     */
    std::tuple<Vector6d, Matrix6d, PlasticState> radialReturn(const Vector6d& elastic_trial_strain, const Eigen::Matrix3d& betr_n1, const PlasticState& last_state) const
    {
        const double G = _E / (2 * (1+_nu));  // shear modulus
        const double K = _E / (3 * (1-2*_nu));    // bulk modulus

        // 1. Compute trial elastic stress
        Eigen::Matrix3d elastic_trial_strain_mat = strainVoigt2Mat(elastic_trial_strain);
        const double elastic_strain_v = 1.0/3.0 * elastic_trial_strain_mat.trace();
        Eigen::Matrix3d e_e_trial = elastic_trial_strain_mat - 1.0/3.0 * elastic_trial_strain_mat.trace() * Eigen::Matrix3d::Identity();
        Eigen::Matrix3d beta_old = stressVoigt2Mat(last_state.beta);
        double alpha_old = last_state.alpha;

        // Eigen::Matrix3d e_new = strain_mat_new - 1.0/3.0 * strain_mat_new.trace() * Eigen::Matrix3d::Identity();   // deviatoric part of new strain
        Eigen::Matrix3d s_new_trial = 2*G*e_e_trial;
        const double p_new = K * elastic_trial_strain_mat.trace();
        Eigen::Matrix3d xi_new_trial = s_new_trial - beta_old;

        // 2. Check yield condition
        double f_new_trial = xi_new_trial.norm() - std::sqrt(2.0/3.0) * _isotropicHardeningModulus(alpha_old);
        // std::cout << "s_trial_11: " << s_new_trial << "   ";
        // std::cout << "f_trial: " << f_new_trial << "   ";

        PlasticState new_state;
        Matrix6d D_ep;
        Vector6d stress;
        // we are not plastically deforming
        if (f_new_trial <= 0)
        {
            new_state = last_state;
            new_state.elastic_strain = elastic_trial_strain;

            stress = stressMat2Voigt(K * elastic_trial_strain_mat.trace() * Eigen::Matrix3d::Identity() + s_new_trial);

            
            D_ep << 1-_nu, _nu, _nu, 0, 0, 0,
                    _nu, 1-_nu, _nu, 0, 0, 0,
                    _nu, _nu, 1-_nu, 0, 0, 0,
                    0, 0, 0, 0.5*(1-2*_nu), 0, 0,
                    0, 0, 0, 0, 0.5*(1-2*_nu), 0,
                    0, 0, 0, 0, 0, 0.5*(1-2*_nu);
            D_ep *= _E / ( (1+_nu) * (1-2*_nu) );
            // return std::make_tuple(stress, D_ep, new_state);
        }
        else
        {

            // plastic deformation, project elastic state onto plastic yield surface
            // 3. Compute n_n+1 and delta_gamma
            Eigen::Matrix3d n_new = xi_new_trial / xi_new_trial.norm();
            double dgamma = f_new_trial / ( 2*G * ( 1 + _H_bar_prime / (3*G)) );

            // 4. Update quantities
            double alpha_new = alpha_old + std::sqrt(2.0/3.0) * dgamma;
            Eigen::Matrix3d beta_new = beta_old + std::sqrt(2.0/3.0) * (_kinematicHardeningModulus(alpha_new) - _kinematicHardeningModulus(alpha_old)) * n_new;
            // Eigen::Matrix3d e_p_new = e_p_old + dgamma * n_new;
            Eigen::Matrix3d e_e_new = e_e_trial - dgamma * n_new;
            Eigen::Matrix3d sigma_new = p_new * Eigen::Matrix3d::Identity() + s_new_trial - 2*G*dgamma*n_new;

            // 5. Compute consistent elastoplastic tangent moduli
            double theta = 1 - (2*G*dgamma) / xi_new_trial.norm();
            double theta_bar = 1 / (1 + (_isotropicHardeningModulusPrime() + _kinematicHardeningModulusPrime()) / (3*G) ) - (1 - theta);
            
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
            new_state.elastic_strain = strainMat2Voigt(e_e_new + elastic_strain_v*Eigen::Matrix3d::Identity());
            new_state.alpha = alpha_new;
            new_state.beta = stressMat2Voigt(beta_new);
            stress = stressMat2Voigt(sigma_new);
        }

        Matrix6d D = _getCijkl691(betr_n1, stress, D_ep);
        // new_state.dev_plastic_strain = strainMat2Voigt(e_p_new);

        // std::cout << "New elastic strain:\n" << new_state.elastic_strain << std::endl;

        return std::make_tuple(stress, D, new_state);
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

    double _isotropicHardeningModulusPrime() const
    {
        return _beta * _H_bar_prime;
    }

    double _kinematicHardeningModulusPrime() const
    {
        return (1-_beta) * _H_bar_prime;
    }

    Matrix6d _getL(const Eigen::Matrix3d& V, const Eigen::Vector3d& D, const Eigen::Vector3d& Dy, const Eigen::Vector3d& Ddy) const
    {
        // computes d_ln(be)/d_be using formulas from deSouza book, Appendix 1
        // Usage: L = getL(Vn1,betr_v,eetr_v2,deetr_v);
        //    Vn1 = eigenvector basis of be_trial
        //    betr_v = 3 eigenvalues of be_trial
        //      These should be obtained from [Vn1,betr_v] = eig(betr_n1) using
        //      Matlab's eigenvalue solver
        //    eetr_v2 = diag(log(betr_v));
        //    deetr_v = diag(1./(betr_v));
    
        Matrix6d L = Matrix6d::Zero();
        Matrix6d Is = Matrix6d::Identity();
        Vector6d One(6);
        One << 1, 1, 1, 0, 0, 0;
    
        Eigen::Vector3d d = D;
        Eigen::Matrix3d x = V * D.asDiagonal() * V.transpose();
        Vector6d xv;
        xv << x(0, 0), x(1, 1), x(2, 2), x(0, 1), x(1, 2), x(2, 0);
        Eigen::Vector3d dy = Dy;
        Eigen::Vector3d ddy = Ddy;
    
        Eigen::Matrix<double, 6, 3> E;
        E << V(0, 0) * V(0, 0), V(0, 1) * V(0, 1), V(0, 2) * V(0, 2),
             V(1, 0) * V(1, 0), V(1, 1) * V(1, 1), V(1, 2) * V(1, 2),
             V(2, 0) * V(2, 0), V(2, 1) * V(2, 1), V(2, 2) * V(2, 2),
             V(0, 0) * V(1, 0), V(0, 1) * V(1, 1), V(0, 2) * V(1, 2),
             V(1, 0) * V(2, 0), V(1, 1) * V(2, 1), V(1, 2) * V(2, 2),
             V(2, 0) * V(0, 0), V(2, 1) * V(0, 1), V(2, 2) * V(0, 2);
    
        Matrix6d dx2dx;
        dx2dx << 2 * x(0, 0), 0, 0, x(0, 1), 0, x(0, 2),
                  0, 2 * x(1, 1), 0, x(0, 1), x(1, 2), 0,
                  0, 0, 2 * x(2, 2), 0, x(1, 2), x(2, 0),
                  x(0, 1), x(0, 1), 0, (x(0, 0) + x(1, 1)) / 2, x(0, 2) / 2, x(1, 2) / 2,
                  0, x(1, 2), x(1, 2), x(0, 2) / 2, (x(1, 1) + x(2, 2)) / 2, x(0, 1) / 2,
                  x(0, 2), 0, x(0, 2), x(1, 2) / 2, x(0, 1) / 2, (x(0, 0) + x(2, 2)) / 2;
    
        double tol = 1e-10;
        Eigen::MatrixXd abc(3, 2);
        abc << 2, 3,
               3, 1,
               1, 2;
    
        if (abs(d(0) - d(1)) > tol && abs(d(0) - d(2)) > tol && abs(d(1) - d(2)) > tol) {
            for (int a = 0; a < 3; a++) {
                int b = abc(a, 0) - 1;
                int c = abc(a, 1) - 1;
    
                double xa = d(a);
                double xb = d(b);
                double xc = d(c);
                double yfac = dy(a) / ((xa - xb) * (xa - xc));
    
                L += yfac * (dx2dx - (xb + xc) * Is - (2 * xa - xb - xc) * E.col(a) * E.col(a).transpose() -
                              (xb - xc) * (E.col(b) * E.col(b).transpose() - E.col(c) * E.col(c).transpose())) +
                              ddy(a) * E.col(a) * E.col(a).transpose();
            }
        } else if (abs(d(0) - d(1)) <= tol && abs(d(0) - d(2)) <= tol) {
            L = ddy(0) * Is;
        } else if (abs(d(1) - d(2)) <= tol) {
            int a = 0;
            int c = abc(a, 1) - 1;
            double xa = d(a);
            double xc = d(c);
    
            double xaxc = (xa - xc);
            double xaxc2 = xaxc * xaxc;
            double yayc = (dy(a) - dy(c));
            double s1 = yayc / xaxc2 - ddy(c) / xaxc;
            double s2 = 2 * xc * yayc / xaxc2 - (xa + xc) / xaxc * ddy(c);
            double s3 = 2 * yayc / (xaxc2 * xaxc) - (ddy(a) + ddy(c)) / xaxc2;
            double s4 = xc * s3;
            double s6 = xc * s4;
    
            L = s1 * dx2dx - s2 * Is - s3 * (xv * xv.transpose()) + s4 * (xv * One.transpose() + One * xv.transpose()) - s6 * (One * One.transpose());
        } else if (abs(d(0) - d(2)) <= tol) {
            int a = 1;
            int c = abc(a, 1) - 1;
            double xa = d(a);
            double xc = d(c);
    
            double xaxc = (xa - xc);
            double xaxc2 = xaxc * xaxc;
            double yayc = (dy(a) - dy(c));
            double s1 = yayc / xaxc2 - ddy(c) / xaxc;
            double s2 = 2 * xc * yayc / xaxc2 - (xa + xc) / xaxc * ddy(c);
            double s3 = 2 * yayc / (xaxc2 * xaxc) - (ddy(a) + ddy(c)) / xaxc2;
            double s4 = xc * s3;
            double s6 = xc * s4;
    
            L = s1 * dx2dx - s2 * Is - s3 * (xv * xv.transpose()) + s4 * (xv * One.transpose() + One * xv.transpose()) - s6 * (One * One.transpose());
        } else if (abs(d(0) - d(1)) <= tol) {
            int a = 2;
            int c = abc(a, 1) - 1;
            double xa = d(a);
            double xc = d(c);
    
            double xaxc = (xa - xc);
            double xaxc2 = xaxc * xaxc;
            double yayc = (dy(a) - dy(c));
            double s1 = yayc / xaxc2 - ddy(c) / xaxc;
            double s2 = 2 * xc * yayc / xaxc2 - (xa + xc) / xaxc * ddy(c);
            double s3 = 2 * yayc / (xaxc2 * xaxc) - (ddy(a) + ddy(c)) / xaxc2;
            double s4 = xc * s3;
            double s6 = xc * s4;
    
            L = s1 * dx2dx - s2 * Is - s3 * (xv * xv.transpose()) + s4 * (xv * One.transpose() + One * xv.transpose()) - s6 * (One * One.transpose());
        }
    
        return L;
    }

    Matrix6d _getCijkl691(const Eigen::Matrix3d& betr_n1, const Vector6d& s_n1, const Matrix6d& c_n1) const
    {
        // Compute spectral decomposition
        Eigen::EigenSolver<Eigen::Matrix3d> solver(betr_n1);
        Eigen::Matrix3d Vn1 = solver.eigenvectors().real();
        Eigen::Vector3d betr_v = solver.eigenvalues().real();
    
        // Compute Be_n = exp(2*ee_n)
        Eigen::Vector3d eetr_v2 = betr_v.array().log();
        // Eigen::Vector3d eetr_v = 0.5 * eetr_v2;
        Eigen::Vector3d deetr_v = betr_v.array().inverse();
    
        // Compute extra 4-order tensors L and B
        Matrix6d L = _getL(Vn1, betr_v, eetr_v2, deetr_v);
    
        Eigen::Vector<double, 9> vec;
        vec << 0, 1, 2, 3, 4, 5, 3, 4, 5;
        Eigen::Matrix<double, 9, 9> C_n12 = c_n1(vec, vec);
        Eigen::Matrix<double, 9, 9> L2 = L(vec, vec);
    
        Eigen::Matrix3d x = betr_n1;
        Eigen::Matrix<double, 9, 9> B;
        B << 2 * x(0, 0), 0, 0, x(0, 1), 0, x(0, 2), x(0, 1), 0, x(0, 2),
             0, 2 * x(1, 1), 0, x(0, 1), x(1, 2), 0, x(0, 1), x(1, 2), 0,
             0, 0, 2 * x(2, 2), 0, x(1, 2), x(0, 2), 0, x(1, 2), x(0, 2),
             x(0, 1), x(0, 1), 0, x(0, 0) / 2 + x(1, 1) / 2, x(0, 2) / 2, x(1, 2) / 2, x(0, 0) / 2 + x(1, 1) / 2, x(0, 2) / 2, x(1, 2) / 2,
             0, x(1, 2), x(1, 2), x(0, 2) / 2, x(1, 1) / 2 + x(2, 2) / 2, x(0, 1) / 2, x(0, 2) / 2, x(1, 1) / 2 + x(2, 2) / 2, x(0, 1) / 2,
             x(0, 2), 0, x(0, 2), x(1, 2) / 2, x(0, 1) / 2, x(0, 0) / 2 + x(2, 2) / 2, x(1, 2) / 2, x(0, 1) / 2, x(0, 0) / 2 + x(2, 2) / 2,
             x(0, 1), x(0, 1), 0, x(0, 0) / 2 + x(1, 1) / 2, x(0, 2) / 2, x(1, 2) / 2, x(0, 0) / 2 + x(1, 1) / 2, x(0, 2) / 2, x(1, 2) / 2,
             0, x(1, 2), x(1, 2), x(0, 2) / 2, x(1, 1) / 2 + x(2, 2) / 2, x(0, 1) / 2, x(0, 2) / 2, x(1, 1) / 2 + x(2, 2) / 2, x(0, 1) / 2,
             x(0, 2), 0, x(0, 2), x(1, 2) / 2, x(0, 1) / 2, x(0, 0) / 2 + x(2, 2) / 2, x(1, 2) / 2, x(0, 1) / 2, x(0, 0) / 2 + x(2, 2) / 2;
    
        Eigen::Matrix<double, 9, 9> c_n1_9x9 = 0.5 * C_n12 * L2 * B;
        Vector6d sigma2 = s_n1;
        Eigen::Matrix<double, 9, 9> Smat;
        Smat << sigma2(0), 0, 0, sigma2(3) / 2, 0, sigma2(5) / 2, sigma2(3) / 2, 0, -sigma2(5) / 2,
                0, sigma2(1), 0, sigma2(3) / 2, sigma2(4) / 2, 0, -sigma2(3) / 2, sigma2(4) / 2, 0,
                0, 0, sigma2(2), 0, sigma2(4) / 2, sigma2(5) / 2, 0, -sigma2(4) / 2, sigma2(5) / 2,
                sigma2(3) / 2, sigma2(3) / 2, 0, sigma2(0) / 4 + sigma2(1) / 4, sigma2(5) / 4, sigma2(4) / 4, sigma2(1) / 4 - sigma2(0) / 4, sigma2(5) / 4, -sigma2(4) / 4,
                0, sigma2(4) / 2, sigma2(4) / 2, sigma2(5) / 4, sigma2(1) / 4 + sigma2(2) / 4, sigma2(3) / 4, -sigma2(5) / 4, sigma2(2) / 4 - sigma2(1) / 4, sigma2(3) / 4,
                sigma2(5) / 2, 0, sigma2(5) / 2, sigma2(4) / 4, sigma2(3) / 4, sigma2(0) / 4 + sigma2(2) / 4, sigma2(4) / 4, -sigma2(3) / 4, sigma2(0) / 4 - sigma2(2) / 4,
                sigma2(3) / 2, -sigma2(3) / 2, 0, sigma2(1) / 4 - sigma2(0) / 4, -sigma2(5) / 4, sigma2(4) / 4, sigma2(0) / 4 + sigma2(1) / 4, -sigma2(5) / 4, -sigma2(4) / 4,
                0, sigma2(4) / 2, -sigma2(4) / 2, sigma2(5) / 4, sigma2(2) / 4 - sigma2(1) / 4, -sigma2(3) / 4, -sigma2(5) / 4, sigma2(1) / 4 + sigma2(2) / 4, -sigma2(3) / 4,
                -sigma2(5) / 2, 0, sigma2(5) / 2, -sigma2(4) / 4, sigma2(3) / 4, sigma2(0) / 4 - sigma2(2) / 4, -sigma2(5) / 4, -sigma2(4) / 4, sigma2(0) / 4 + sigma2(2) / 4;
    
        c_n1_9x9 = c_n1_9x9 - 2 * Smat;
        Matrix6d D_n1 = c_n1.block<6,6>(0, 0);
    
        return D_n1;
    }
    
    
    double _E;
    double _nu;
    double _sigma_y;
    double _beta;
    double _H_bar_prime;
};