/** Implements the Radial Return algorithm for HW 6 */

#include "common.hpp"
#include "thirdparty/gplot++.h"

#include <iostream>

constexpr double sigma_y = 22;  // yield stress
constexpr double beta = 2.5/6;
constexpr double H_bar_prime = 6;
constexpr double E = 1200;
constexpr double nu = 0.25;
constexpr double G = E / (2 * (1+nu));  // shear modulus
constexpr double K = E / (3 * (1-2*nu));    // bulk modulus

constexpr double dt = 0.05;
constexpr double t_end = 3;

Eigen::Matrix3d strainVoigt2Mat(const Vector6d& strain_vec)
{
    Eigen::Matrix3d strain_mat;
    strain_mat <<       strain_vec[0], 0.5*strain_vec[5], 0.5*strain_vec[4],
                        0.5*strain_vec[5], strain_vec[1], 0.5*strain_vec[3],
                        0.5*strain_vec[4], 0.5*strain_vec[3], strain_vec[2];

    return strain_mat;
}

Vector6d strainMat2Voigt(const Eigen::Matrix3d& strain_mat)
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

Eigen::Matrix3d stressVoigt2Mat(const Vector6d& stress_vec)
{
    Eigen::Matrix3d stress_mat;
    stress_mat <<       stress_vec[0], stress_vec[5], stress_vec[4],
                        stress_vec[5], stress_vec[1], stress_vec[3],
                        stress_vec[4], stress_vec[3], stress_vec[2];
    return stress_mat;
}

Vector6d stressMat2Voigt(const Eigen::Matrix3d& stress_mat)
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

Vector6d strain(double t)
{
    const double eps_1 = 0.06;
    const double eps_2 = 0.09;

    Eigen::Matrix3d mat_1, mat_2;
    mat_1 << 1, 0, 0,
             0, 1, 0,
             0, 0, 0.8;
    mat_2 << 0, 2, 0,
             2, 1, 0,
             0, 0, 5.0/3.0;
    Eigen::Matrix3d strain_mat = eps_1 * t * mat_1 + eps_2 * std::sin(t) * mat_2;
    

    return strainMat2Voigt(strain_mat);
}

double isotropicHardeningModulus(double alpha)
{
    return sigma_y + beta * H_bar_prime * alpha;
}

double isotropicHardeningModulusPrime(double dgamma)
{
    return beta * H_bar_prime * dgamma / dt * std::sqrt(2.0/3.0);
}

double kinematicHardeningModulus(double alpha)
{
    return (1-beta) * H_bar_prime * alpha;
}

double kinematicHardeningModulusPrime(double dgamma)
{
    return (1-beta) * H_bar_prime * dgamma / dt * std::sqrt(2.0/3.0);
}

/** Given new strain and the last plastic state, returns the
 * 1. 6x1 new stress
 * 2. 6x6 consistent tangent moduli
 * 3. the new plastic state
 */
std::tuple<Vector6d, Matrix6d, PlasticState> radialReturn(const Vector6d& new_strain, const PlasticState& last_state)
{
    // 1. Compute trial elastic stress
    Eigen::Matrix3d strain_mat_new = strainVoigt2Mat(new_strain);
    Eigen::Matrix3d e_p_old = strainVoigt2Mat(last_state.dev_plastic_strain);
    Eigen::Matrix3d beta_old = stressVoigt2Mat(last_state.beta);
    double alpha_old = last_state.alpha;

    Eigen::Matrix3d e_new = strain_mat_new - 1.0/3.0 * strain_mat_new.trace() * Eigen::Matrix3d::Identity();   // deviatoric part of new strain
    Eigen::Matrix3d s_new_trial = 2*G*(e_new - e_p_old);
    Eigen::Matrix3d xi_new_trial = s_new_trial - beta_old;

    // 2. Check yield condition
    double f_new_trial = xi_new_trial.norm() - std::sqrt(2.0/3.0) * isotropicHardeningModulus(alpha_old);

    std::cout << "f_trial: " << f_new_trial << std::endl;

    PlasticState new_state;
    // we are not plastically deforming
    if (f_new_trial <= 0)
    {
        new_state.alpha = alpha_old;
        new_state.beta = last_state.beta;

        Vector6d stress = stressMat2Voigt(K * strain_mat_new.trace() * Eigen::Matrix3d::Identity() + s_new_trial);
        new_state.dev_plastic_strain = strainMat2Voigt(e_p_old);

        Matrix6d D_ep;
        D_ep << 1-nu, nu, nu, 0, 0, 0,
                nu, 1-nu, nu, 0, 0, 0,
                nu, nu, 1-nu, 0, 0, 0,
                0, 0, 0, 0.5*(1-2*nu), 0, 0,
                0, 0, 0, 0, 0.5*(1-2*nu), 0,
                0, 0, 0, 0, 0, 0.5*(1-2*nu);
        D_ep *= E / ( (1+nu) * (1-2*nu) );
        return std::make_tuple(stress, D_ep, new_state);
    }

    // plastic deformation, project elastic state onto plastic yield surface
    // 3. Compute n_n+1 and delta_gamma
    Eigen::Matrix3d n_new = xi_new_trial / xi_new_trial.norm();
    double dgamma = f_new_trial / ( 2*G * ( 1 + H_bar_prime / (3*G)) );

    std::cout << "dgamma: " << dgamma << std::endl;

    // 4. Update quantities
    double alpha_new = alpha_old + std::sqrt(2.0/3.0) * dgamma;
    std::cout << "alpha new: " << alpha_new << std::endl;
    Eigen::Matrix3d beta_new = beta_old + std::sqrt(2.0/3.0) * (kinematicHardeningModulus(alpha_new) - kinematicHardeningModulus(alpha_old)) * n_new;
    Eigen::Matrix3d e_p_new = e_p_old + dgamma * n_new;
    Eigen::Matrix3d sigma_new = K * strain_mat_new.trace() * Eigen::Matrix3d::Identity() + s_new_trial - 2*G*dgamma*n_new;

    // 5. Compute consistent elastoplastic tangent moduli
    double theta = 1 - (2*G*dgamma) / xi_new_trial.norm();
    double theta_bar = 1 / (1 + (isotropicHardeningModulusPrime(dgamma) + kinematicHardeningModulusPrime(dgamma)) / (3*G) ) - (1 - theta);
    
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

    // test that f=0
    Eigen::Matrix3d s_new_trial_test = 2*G*(e_new - e_p_new);
    Eigen::Matrix3d xi_new_trial_test = s_new_trial_test - beta_new;

    // 2. Check yield condition
    double f_new_trial_test = xi_new_trial_test.norm() - std::sqrt(2.0/3.0) * isotropicHardeningModulus(alpha_new);
    std::cout << "After updates: f=" << f_new_trial_test << std::endl;

    return std::make_tuple(stress, D_ep, new_state);
}


int main() 
{   
    std::vector<PlasticState> plastic_states;
    std::vector<Vector6d> stresses;
    std::vector<Matrix6d> D_eps;
    std::vector<double> times;

    // set initial plastic state
    PlasticState plastic_state;
    plastic_state.dev_plastic_strain = Vector6d::Zero();
    plastic_state.beta = Vector6d::Zero();
    plastic_state.alpha = 0;
    plastic_states.push_back(plastic_state);
    
    for (double t = 0.0; t <= t_end; t+=dt)
    {
        Vector6d cur_strain = strain(t);
        const auto [stress, D_ep, new_plastic_state] = radialReturn(cur_strain, plastic_state);
        plastic_state = new_plastic_state;
        stresses.push_back(stress);
        plastic_states.push_back(new_plastic_state);
        D_eps.push_back(D_ep); 
        times.push_back(t);
    }

    // plot using GNU plot
    Gnuplot plt{};

    // plt.sendcommand("set terminal wxt size 1500,500"); 
    plt.multiplot(1, 4, "HW6 (a)");
    // plot 1 - sigma vs time
    plt.set_title("sigma vs. time");
    plt.set_xlabel("Time");
    plt.set_ylabel("sigma");
    std::vector<double> sigma_11, sigma_22, sigma_12, sigma_33;
    for (const auto& stress : stresses) { sigma_11.push_back(stress[0]); sigma_22.push_back(stress[1]); sigma_12.push_back(stress[5]); sigma_33.push_back(stress[2]); }
    plt.plot(times, sigma_11, "sigma_{11}");
    plt.plot(times, sigma_22, "sigma_{22}");
    plt.plot(times, sigma_12, "sigma_{12}");
    plt.plot(times, sigma_33, "sigma_{33}");
    plt.show();
    // plot 2 - alpha vs time
    plt.set_title("alpha vs time");
    plt.set_xlabel("time");
    plt.set_ylabel("alpha");
    std::vector<double> alphas;
    for (const auto& state : plastic_states) { alphas.push_back(state.alpha); }
    plt.plot(times, alphas, "alpha");
    plt.show();
    // plot 3 - displacement vs time
    plt.set_title("beta vs. time");
    plt.set_xlabel("time");
    plt.set_ylabel("beta");
    std::vector<double> beta_11, beta_22, beta_12, beta_33;
    for (const auto& state : plastic_states) { beta_11.push_back(state.beta[0]); beta_22.push_back(state.beta[1]); beta_12.push_back(state.beta[5]); beta_33.push_back(state.beta[2]); }
    plt.plot(times, beta_11, "beta_{11}");
    plt.plot(times, beta_22, "beta_{22}");
    plt.plot(times, beta_12, "beta_{12}");
    plt.plot(times, beta_33, "beta_{33}");
    plt.show();
    // plot 3 - moduli vs time
    plt.set_title("C_{ep} vs. time");
    plt.set_xlabel("time");
    plt.set_ylabel("C_{ep}");
    std::vector<double> C_1111, C_2222, C_1212;
    for (const auto& D_ep : D_eps) { C_1111.push_back(D_ep(0,0)); C_2222.push_back(D_ep(1,1)); C_1212.push_back(D_ep(5,5)); }
    plt.plot(times, C_1111, "C_{1111}");
    plt.plot(times, C_2222, "C_{2222}");
    plt.plot(times, C_1212, "C_{1212}");
    plt.show();
}