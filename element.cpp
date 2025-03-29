#include "element.hpp"

Eigen::Matrix<double, 2, 8> QuadElement::H(double r, double s) const
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

Eigen::Matrix<double, 3, 8> QuadElement::B(double r, double s, const Eigen::VectorXd& d_e) const
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
    const Eigen::Matrix2d J_mat = DeformedJacobian(r, s, d_e);
    // const Eigen::Matrix2d J_mat = UndeformedJacobian(r, s);
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

Eigen::Matrix2d QuadElement::DeformedJacobian(double r, double s, const Eigen::VectorXd& d_e) const
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

    const Eigen::Vector2d x1_cur = _x1 + Eigen::Vector2d(d_e[0], d_e[1]);
    const Eigen::Vector2d x2_cur = _x2 + Eigen::Vector2d(d_e[2], d_e[3]);
    const Eigen::Vector2d x3_cur = _x3 + Eigen::Vector2d(d_e[4], d_e[5]);
    const Eigen::Vector2d x4_cur = _x4 + Eigen::Vector2d(d_e[6], d_e[7]);
    J_mat <<    dh1_dr*x1_cur[0] + dh2_dr*x2_cur[0] + dh3_dr*x3_cur[0] + dh4_dr*x4_cur[0],
                dh1_dr*x1_cur[1] + dh2_dr*x2_cur[1] + dh3_dr*x3_cur[1] + dh4_dr*x4_cur[1],
                dh1_ds*x1_cur[0] + dh2_ds*x2_cur[0] + dh3_ds*x3_cur[0] + dh4_ds*x4_cur[0],
                dh1_ds*x1_cur[1] + dh2_ds*x2_cur[1] + dh3_ds*x3_cur[1] + dh4_ds*x4_cur[1];
    
    return J_mat;
}

Eigen::Matrix2d QuadElement::UndeformedJacobian(double r, double s) const
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
    J_mat <<    dh1_dr*_x1[0] + dh2_dr*_x2[0] + dh3_dr*_x3[0] + dh4_dr*_x4[0],
                dh1_dr*_x1[1] + dh2_dr*_x2[1] + dh3_dr*_x3[1] + dh4_dr*_x4[1],
                dh1_ds*_x1[0] + dh2_ds*_x2[0] + dh3_ds*_x3[0] + dh4_ds*_x4[0],
                dh1_ds*_x1[1] + dh2_ds*_x2[1] + dh3_ds*_x3[1] + dh4_ds*_x4[1];
    
    return J_mat;
}

// Eigen::Matrix<double, 8, 8> QuadElement::M() const
// {
//     Eigen::Matrix<double, 8, 8> M_mat = Eigen::Matrix<double, 8, 8>::Zero();

//     for (unsigned i = 0; i < _integration_points.size(); i++)
//     {
//         const double ri = _integration_points[i];
//         const double wi = _integration_weights[i];
//         for (unsigned j = 0; j < _integration_points.size(); j++)
//         {
//             const double sj = _integration_points[j];
//             const double wj = _integration_weights[j];
//             const Eigen::Matrix2d J_mat = J(ri, sj);
//             const Eigen::Matrix<double, 2, 8> H_mat = H(ri, sj);
//             M_mat += wi * wj * _density * J_mat.determinant() * H_mat.transpose() * H_mat;
//         }
//     }

//     return M_mat;
// }

Eigen::Matrix<double, 8, 8> QuadElement::K(const Eigen::VectorXd& d_e) const
{
    assert(d_e.size() == NSDIMS*numNodes());

    Eigen::Matrix<double, 8, 8> K_c = Eigen::Matrix<double, 8, 8>::Zero();
    Eigen::Matrix<double, 8, 8> K_sigma = Eigen::Matrix<double, 8, 8>::Zero();
    
    for (unsigned i = 0; i < _integration_points.size(); i++)
    {
        const double ri = _integration_points[i];
        const double wi = _integration_weights[i];
        for (unsigned j = 0; j < _integration_points.size(); j++)
        {
            const double sj = _integration_points[j];
            const double wj = _integration_weights[j];
            // find deformation gradient at (r,s) given the current deformation
            const Eigen::Matrix2d F_mat = deformationGradient(ri, sj, d_e);
            const Eigen::Matrix2d J_mat = DeformedJacobian(ri, sj, d_e);
            // const Eigen::Matrix2d J_mat = UndeformedJacobian(ri, sj);
            const Eigen::Matrix<double, 3, 8> B_mat = B(ri, sj, d_e);

            const auto [stress_vec, D_mat] = _material->materialSubroutine(F_mat);
            K_c += wi * wj * B_mat.transpose() * D_mat * B_mat * J_mat.determinant();

            // NEW STUFF
            // derivative of interpolation funcs wrt r
            Eigen::Vector4d dh_dr;
            dh_dr(0) = 0.25 * (1 + sj);
            dh_dr(1) = -0.25 * (1 + sj);
            dh_dr(2) = -0.25 * (1 - sj);
            dh_dr(3) = 0.25 * (1 - sj);

            // derivative of interpolation funcs wrt s
            Eigen::Vector4d dh_ds;
            dh_ds(0) = 0.25 * (1 + ri);
            dh_ds(1) = 0.25 * (1 - ri);
            dh_ds(2) = -0.25 * (1 - ri);
            dh_ds(3) = -0.25 * (1 + ri);

            // get the Jacobian operator and invert it
            const Eigen::Matrix2d J_inv = J_mat.inverse();

            // derivative of the shape functions wrt current coordiantes, x
            Eigen::Matrix<double, 2, 4> dH_dx = Eigen::Matrix<double, 2, 4>::Zero();
            for (int i = 0; i < 4; i++)
            {
                dH_dx(0,i) = J_inv(0,0)*dh_dr(i) + J_inv(0,1)*dh_ds(i);
                dH_dx(1,i) = J_inv(1,0)*dh_dr(i) + J_inv(1,1)*dh_ds(i);
            }

            Eigen::Matrix2d sigma_mat;
            sigma_mat(0,0) = stress_vec[0];
            sigma_mat(1,1) = stress_vec[1];
            sigma_mat(0,1) = stress_vec[2];
            sigma_mat(1,0) = sigma_mat(0,1);
            for (int a = 0; a < 4; a++)
            {
                Eigen::Vector2d grad_Na = dH_dx.col(a);
                for (int b = 0; b < 4; b++)
                {
                    Eigen::Vector2d grad_Nb = dH_dx.col(b);
                    double val = wi * wj * (grad_Na.transpose() * sigma_mat * grad_Nb)(0,0) * J_mat.determinant();
                    K_sigma.block<2,2>(2*a, 2*b) += val * Eigen::Matrix2d::Identity();
                }
            }
        }
    }

    // std::cout << "K_c:\n" << K_c << "\nK_sigma:\n" << K_sigma << std::endl;

    return K_c + K_sigma;
}

Eigen::Vector<double, 8> QuadElement::internalForce(const Eigen::VectorXd& d_e) const
{
    assert(d_e.size() == NSDIMS*numNodes());

    Eigen::Vector<double, 8> R_vec = Eigen::Vector<double, 8>::Zero();

    for (unsigned i = 0; i < _integration_points.size(); i++)
    {
        const double ri = _integration_points[i];
        const double wi = _integration_weights[i];
        for (unsigned j = 0; j < _integration_points.size(); j++)
        {
            const double sj = _integration_points[j];
            const double wj = _integration_weights[j];
            // find deformation gradient at (r,s) given the current deformation
            const Eigen::Matrix2d F_mat = deformationGradient(ri, sj, d_e);
            const Eigen::Matrix2d J_mat = DeformedJacobian(ri, sj, d_e);
            const Eigen::Matrix<double, 3, 8> B_mat = B(ri, sj, d_e);

            const auto [stress_vec, D_mat] = _material->materialSubroutine(F_mat);
            R_vec += wi * wj * B_mat.transpose() * stress_vec * J_mat.determinant();
        }
    }

    return R_vec;
}

Eigen::Matrix2d QuadElement::deformationGradient(double r, double s, const Eigen::VectorXd& d_e) const
{
    assert(d_e.size() == NSDIMS*numNodes());
    Eigen::Matrix2d F;

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
    // const Eigen::Matrix2d J_mat = J(r, s, d_e);
    const Eigen::Matrix2d J_undeformed = UndeformedJacobian(r, s);
    const Eigen::Matrix2d J_inv = J_undeformed.inverse();

    // derivative of the shape functions wrt current coordiantes, x
    Eigen::Matrix<double, 2, 4> dH_dx = Eigen::Matrix<double, 2, 4>::Zero();
    for (int i = 0; i < 4; i++)
    {
        dH_dx(0,i) = J_inv(0,0)*dh_dr(i) + J_inv(0,1)*dh_ds(i);
        dH_dx(1,i) = J_inv(1,0)*dh_dr(i) + J_inv(1,1)*dh_ds(i);
    }

    Eigen::Vector4d d_x(d_e[0], d_e[2], d_e[4], d_e[6]);
    Eigen::Vector4d d_y(d_e[1], d_e[3], d_e[5], d_e[7]);

    Eigen::Matrix2d du_dX;
    du_dX.row(0) = dH_dx * d_x;
    du_dX.row(1) = dH_dx * d_y;
    

    F = Eigen::Matrix2d::Identity() + du_dX;

    return F;
}