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
          _density(density), _E(E), _mu(mu),
          _integration_points({-1.0/std::sqrt(3), 1.0/std::sqrt(3)}),
          _integration_weights({1.0, 1.0})
    {
    }

    // evaluates the interpolation matrix H at (r,s)
    Eigen::Matrix<double, 2, 8> H(double r, double s) const
    {
        // interpolation functions
        double h1 = 0.25 * (1 + r) * (1 + s);
        double h2 = 0.25 * (1 - r) * (1 + s);
        double h3 = 0.25 * (1 - r) * (1 - s);
        double h4 = 0.25 * (1 + r) * (1 - s);

        Eigen::Matrix<double, 2, 8> H_mat;
        H_mat << h1, 0,  h2, 0,  h3, 0,  h4, 0,
                 0,  h1, 0,  h2, 0,  h3, 0,  h4;

        return H_mat;
    }

    // evaluates the strain-displacement matrix B at (r,s)
    Eigen::Matrix<double, 3, 8> B(double r, double s) const
    {
        // derivative of interpolation funcs wrt r
        Eigen::Vector4d dh_dr;
        dh_dr(0) = 0.25 * (1 + s);
        dh_dr(1) = -0.25 * (1 + s);
        dh_dr(2) = -0.25 * (1 - s);
        dh_dr(3) = 0.25 * (1 - s);

        // derivative of interpolation funcs wrt s
        Eigen::Vector4d dh_ds;
        dh_ds(0) = 0.25 * (1 + r);
        dh_ds(1) = 0.25 * (1 - r);
        dh_ds(2) = -0.25 * (1 - r);
        dh_ds(3) = -0.25 * (1 + r);

        // get the Jacobian operator and invert it
        const Eigen::Matrix2d J_mat = J(r, s);
        const Eigen::Matrix2d J_inv = J_mat.inverse();

        // assemble B
        Eigen::Matrix<double, 3, 8> B_mat = Eigen::Matrix<double, 3, 8>::Zero();
        for (int i = 0; i < 4; i++)
        {
            B_mat(0,2*i) = J_inv(0,0)*dh_dr(i) + J_inv(0,1)*dh_ds(i);
            B_mat(1,2*i+1) = J_inv(1,0)*dh_dr(i) + J_inv(1,1)*dh_ds(i);
        }
            
        
        for (int i = 0; i < 4; i++)
        {
            B_mat(2,2*i) = B_mat(1,2*i+1);
            B_mat(2,2*i+1) = B_mat(0,2*i);
        }

        return B_mat;
    }

    // evaluates the Jacobian operator matrix at (r,s)
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
        Eigen::Matrix<double, 8, 8> M_mat = Eigen::Matrix<double, 8, 8>::Zero();

        for (unsigned i = 0; i < _integration_points.size(); i++)
        {
            const double ri = _integration_points[i];
            const double wi = _integration_weights[i];
            for (unsigned j = 0; j < _integration_points.size(); j++)
            {
                const double sj = _integration_points[j];
                const double wj = _integration_weights[j];
                const Eigen::Matrix2d J_mat = J(ri, sj);
                const Eigen::Matrix<double, 2, 8> H_mat = H(ri, sj);
                M_mat += wi * wj * _density * J_mat.determinant() * H_mat.transpose() * H_mat;
            }
        }

        return M_mat;
    }

    // evaluates the elasticity matrix C for the element
    Eigen::Matrix3d C() const
    {
        // elasticity matrix for plane strain
        Eigen::Matrix3d C_mat;
        C_mat << (1-_mu), _mu, 0,
                 _mu, (1-_mu), 0,
                 0, 0, 0.5*(1-2*_mu);
        C_mat *= _E / ( (1+_mu) * (1-2*_mu)); 
        return C_mat;
    }

    // evalutes the stiffness matrix K for the element
    Eigen::Matrix<double, 8, 8> K() const
    {
        Eigen::Matrix<double, 8, 8> K_mat = Eigen::Matrix<double, 8, 8>::Zero();
        const Eigen::Matrix3d C_mat = C();
        for (unsigned i = 0; i < _integration_points.size(); i++)
        {
            const double ri = _integration_points[i];
            const double wi = _integration_weights[i];
            for (unsigned j = 0; j < _integration_points.size(); j++)
            {
                const double sj = _integration_points[j];
                const double wj = _integration_weights[j];
                const Eigen::Matrix2d J_mat = J(ri, sj);
                const Eigen::Matrix<double, 3, 8> B_mat = B(ri, sj);
                K_mat += wi * wj * B_mat.transpose() * C_mat * B_mat * J_mat.determinant();
            }
        }

        return K_mat;
    }

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

    double _density;
    double _E;
    double _mu;

    std::vector<double> _integration_points;
    std::vector<double> _integration_weights;
};

#endif // __ELEMENT_HPP