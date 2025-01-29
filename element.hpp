#ifndef __ELEMENT_HPP
#define __ELEMENT_HPP

#include <Eigen/Dense>
#include <iostream>

// 2D Quadrilateral element with plane strain
class QuadElement
{
    public:
    QuadElement(const Eigen::Vector2d& x1, const Eigen::Vector2d& x2, const Eigen::Vector2d& x3, const Eigen::Vector2d& x4, double density, double E, double mu)
        : _x1(x1), _x2(x2), _x3(x3), _x4(x4),
          _u1(0,0), _u2(0,0), _u3(0,0), _u4(0,0),
          _density(density), _E(E), _mu(mu)
    {

    }

    // evaluates the shape function matrix H at (x,y)
    // Eigen::Matrix<double, 2, 8> H(double x, double y) const
    // {
    //     double h1 = (_x2[0] - x) * (_x4[1] - y) / (_x2[0] - _x1[0]) / (_x4[1] - _x1[1]);
    //     double h2 = (x - _x1[0]) * (_x3[1] - y) / (_x2[0] - _x1[0]) / (_x3[1] - _x2[1]);
    //     double h3 = (_x4[0] - x) * (y - _x2[1]) / (_x4[0] - _x3[0]) / (_x3[1] - _x2[1]);
    //     double h4 = (x - _x3[0]) * (y - _x1[1]) / (_x4[0] - _x3[0]) / (_x4[1] - _x1[1]);

    //     Eigen::Matrix<double, 2, 8> H_mat;
    //     H_mat << h1, h2, h3, h4, 0,  0,  0,  0,
    //              0,  0,  0,  0,  h1, h2, h3, h4;

    //     return H_mat;
    // }
    Eigen::Matrix<double, 2, 8> H(double r, double s) const
    {
        double h1 = 0.25 * (1 + r) * (1 + s);
        double h2 = 0.25 * (1 - r) * (1 + s);
        double h3 = 0.25 * (1 - r) * (1 - s);
        double h4 = 0.25 * (1 + r) * (1 - s);

        Eigen::Matrix<double, 2, 8> H_mat;
        H_mat << h1, h2, h3, h4, 0,  0,  0,  0,
                    0,  0,  0,  0,  h1, h2, h3, h4;

        return H_mat;
    }

    // evalutes the strain-displacement matrix B at (x,y)
    // Eigen::Matrix<double, 3, 8> B(double x, double y) const
    // {
    //     double dh1_dx = -(_x4[1] - y) / (_x2[0] - _x1[0]) / (_x4[1] - _x1[1]);
    //     double dh2_dx = (_x3[1] - y) / (_x2[0] - _x1[0]) / (_x3[1] - _x2[1]);
    //     double dh3_dx = -(y - _x2[1]) / (_x4[0] - _x3[0]) / (_x3[1] - _x2[1]);
    //     double dh4_dx = (y - _x1[1]) / (_x4[0] - _x3[0]) / (_x4[1] - _x1[1]);
    
    //     double dh1_dy = -(_x2[0] - x) / (_x2[0] - _x1[0]) / (_x4[1] - _x1[1]);
    //     double dh2_dy = -(x - _x1[0]) / (_x2[0] - _x1[0]) / (_x3[1] - _x2[1]);
    //     double dh3_dy = (_x4[0] - x) / (_x4[0] - _x3[0]) / (_x3[1] - _x2[1]);
    //     double dh4_dy = (x - _x3[0]) / (_x4[0] - _x3[0]) / (_x4[1] - _x1[1]);
    
    //     Eigen::Matrix<double, 3, 8> B_mat;
    //     B_mat << dh1_dx, dh2_dx, dh3_dx, dh4_dx, 0, 0, 0, 0,
    //              0, 0, 0, 0, dh1_dy, dh2_dy, dh3_dy, dh4_dy,
    //              dh1_dy, dh2_dy, dh3_dy, dh4_dy, dh1_dx, dh2_dx, dh3_dx, dh4_dx;

    //     return B_mat;
    // }
    Eigen::Matrix<double, 3, 8> B(double r, double s) const
    {
        Eigen::Vector4d dh_dr;
        dh_dr(0) = 0.25 * (1 + s);
        dh_dr(1) = -0.25 * (1 + s);
        dh_dr(2) = -0.25 * (1 - s);
        dh_dr(3) = 0.25 * (1 - s);

        Eigen::Vector4d dh_ds;
        dh_ds(0) = 0.25 * (1 + r);
        dh_ds(1) = 0.25 * (1 - r);
        dh_ds(2) = -0.25 * (1 - r);
        dh_ds(3) = -0.25 * (1 + r);

        const Eigen::Matrix2d J_mat = J(r, s);
        const Eigen::Matrix2d J_inv = J_mat.inverse();

        Eigen::Matrix<double, 3, 8> B_mat = Eigen::Matrix<double, 3, 8>::Zero();
        for (int i = 0; i < 4; i++)
            B_mat(0,i) = J_inv(0,0)*dh_dr(i) + J_inv(0,1)*dh_ds(i);
        
        for (int i = 0; i < 4; i++)
            B_mat(1,i+4) = J_inv(1,0)*dh_dr(i) + J_inv(1,1)*dh_ds(i);
        
        for (int i = 0; i < 4; i++)
        {
            B_mat(2,i) = B_mat(1,i+4);
            B_mat(2,i+4) = B_mat(0,i);
        }

        return B_mat;
    }

    Eigen::Matrix2d J(double r, double s) const
    {
        double dh1_dr = 0.25 * (1 + s);
        double dh2_dr = -0.25 * (1 + s);
        double dh3_dr = -0.25 * (1 - s);
        double dh4_dr = 0.25 * (1 - s);

        double dh1_ds = 0.25 * (1 + r);
        double dh2_ds = 0.25 * (1 - r);
        double dh3_ds = -0.25 * (1 - r);
        double dh4_ds = -0.25 * (1 + r);
        const double d_dr = dh1_dr + dh2_dr + dh3_dr + dh4_dr;
        const double d_ds = dh1_ds + dh2_ds + dh3_ds + dh4_ds;

        Eigen::Matrix2d J_mat;
        J_mat << dh1_dr*_x1[0] + dh2_dr*_x2[0] + dh3_dr*_x3[0] + dh4_dr*_x4[0],
                 dh1_dr*_x1[1] + dh2_dr*_x2[1] + dh3_dr*_x3[1] + dh4_dr*_x4[1],
                 dh1_ds*_x1[0] + dh2_ds*_x2[0] + dh3_ds*_x3[0] + dh4_ds*_x4[0],
                 dh1_ds*_x1[1] + dh2_ds*_x2[1] + dh3_ds*_x3[1] + dh4_ds*_x4[1];
        
        return J_mat;
    }

    // evalutes the mass matrix M for the element
    Eigen::Matrix<double, 8, 8> M() const
    {
        std::array<double, 2> gauss_pts({-1.0/std::sqrt(3), 1.0/std::sqrt(3)});
        Eigen::Matrix<double, 8, 8> M_mat = Eigen::Matrix<double, 8, 8>::Zero();
        for (const auto& ri : gauss_pts)
        {
            for (const auto& sj : gauss_pts)
            {
                const Eigen::Matrix2d J_mat = J(ri, sj);
                const Eigen::Matrix<double, 2, 8> H_mat = H(ri, sj);
                M_mat += _density * J_mat.determinant() * H_mat.transpose() * H_mat;
            }
        }

        return M_mat;
    }

    Eigen::Matrix3d C() const
    {
        // elasticity matrix for plane strain
        Eigen::Matrix3d C_mat;
        C_mat <<    _E / (1-_mu*_mu),      _E*_mu / (1-_mu*_mu),   0,
                _E*_mu / (1-_mu*_mu),   _E / (1-_mu*_mu),      0,
                0,                  0,                  0.5*_E*(1-_mu)/(1-_mu*_mu); 
        return C_mat;
    }

    Eigen::Matrix<double, 8, 8> K() const
    {
        std::array<double, 2> gauss_pts({-1.0/std::sqrt(3), 1.0/std::sqrt(3)});
        Eigen::Matrix<double, 8, 8> K_mat = Eigen::Matrix<double, 8, 8>::Zero();
        const Eigen::Matrix3d C_mat = C();
        std::cout << C_mat << std::endl;
        for (const auto& ri : gauss_pts)
        {
            for (const auto& sj : gauss_pts)
            {
                const Eigen::Matrix2d J_mat = J(ri, sj);
                const Eigen::Matrix<double, 3, 8> B_mat = B(ri, sj);
                K_mat += B_mat.transpose() * C_mat * B_mat * J_mat.determinant();
            }
        }

        return K_mat;
    }

    private:
    Eigen::Vector2d _x1;
    Eigen::Vector2d _x2;
    Eigen::Vector2d _x3;
    Eigen::Vector2d _x4;

    Eigen::Vector2d _u1;
    Eigen::Vector2d _u2;
    Eigen::Vector2d _u3;
    Eigen::Vector2d _u4;

    double _density;
    double _E;
    double _mu;
};

#endif // __ELEMENT_HPP