#include "solver.hpp"
#include "thirdparty/gplot++.h"

#include <iomanip>

Solver::Solver( const std::vector<Node> mesh_nodes, const std::vector<ElementNodes> mesh_element_nodes,
                const std::vector<DisplacementBC>& displacement_BCs, const std::vector<ForceBC>& force_BCs,
                const Material* material)
    : _mesh_nodes(mesh_nodes), _mesh_element_nodes(mesh_element_nodes),
     _displacement_BCs(displacement_BCs), _force_BCs(force_BCs),
     _material(material)
     
{
    // verify that displacement BCs and force BCs do not overlap
    for (const auto& disp_bc : displacement_BCs)
    {
        for (const auto& force_bc : force_BCs)
        {
            assert(disp_bc.inputDOF() != force_bc.inputDOF());
        }
    }

    // create elements
    for (unsigned i = 0; i < _mesh_element_nodes.size(); i++)
    {
        const Eigen::Vector4i& e = _mesh_element_nodes[i];
        _elements.push_back(QuadElement(_mesh_nodes[e[0]],
                                        _mesh_nodes[e[1]],
                                        _mesh_nodes[e[2]],
                                        _mesh_nodes[e[3]],
                                        _material));
    }

    // number the global DOF
    _setupDOF();

}

void Solver::_setupDOF()
{
    // initialize the DOF mapping to all -1, so we know which input DOF have been assigned
    _input_to_global_DOF.resize(numDOF(), -1);

    // first start with the DOF at displacement boundary, and put these at the end of DOF numbering
    for (unsigned i = 0; i < _displacement_BCs.size(); i++)
    {
        const DisplacementBC& disp_BC = _displacement_BCs[i];
        _input_to_global_DOF[disp_BC.inputDOF()] = numDOF() - numKnownDisplacements() + i; 
    }
    // then fill in the rest of the DOF
    int cur_DOF = 0;
    for (int i = 0; i < numDOF(); i++)
    {
        // if -1, then this input DOF has not been assigned a global DOF numbering yet
        if (_input_to_global_DOF[i] == -1)
        {
            _input_to_global_DOF[i] = cur_DOF;
            cur_DOF++;
        }
    }

    // set element nodal to global mappings
    for (unsigned i = 0; i < _elements.size(); i++)
    {
        std::vector<int> element_nodal_to_global_DOF(NSDIMS*_elements[i].numNodes());
        const Eigen::Vector4i& e = _mesh_element_nodes[i];
        for (int k = 0; k < 4; k++)
        {
            element_nodal_to_global_DOF[NSDIMS*k] = _input_to_global_DOF[e[k]*NSDIMS];
            element_nodal_to_global_DOF[NSDIMS*k+1] = _input_to_global_DOF[e[k]*NSDIMS+1];
        }
        _elements[i].setGlobalDOF(element_nodal_to_global_DOF);
    }
}

void Solver::solve(int num_load_steps)
{
    // initialize U and R
    _d_global = Eigen::VectorXd::Zero(numDOF());
    _F_ext_global = Eigen::VectorXd::Zero(numDOF());

    // HW6 plotting
    std::vector<double> sigma_11, sigma_22, sigma_12, alpha, beta_11, beta_22, beta_12, C_1111, C_2222, C_1212, times;

    for (int k = 0; k < num_load_steps; k++)
    {
        // put in known displacements and forces
        for (unsigned i = 0; i < _displacement_BCs.size(); i++)
        {
            int global_DOF = _input_to_global_DOF[_displacement_BCs[i].inputDOF()];
            _d_global[global_DOF] = (k+1) * _displacement_BCs[i].displacement / num_load_steps;
        }

        for (unsigned i = 0; i < _force_BCs.size(); i++)
        {
            int global_DOF = _input_to_global_DOF[_force_BCs[i].inputDOF()];
            _F_ext_global[global_DOF] = (k+1) * _force_BCs[i].force / num_load_steps;
        }

        Eigen::VectorXd U_K_globalnown = _d_global(Eigen::seq(numUnknownDisplacements(), numDOF()-1));
        Eigen::VectorXd R_K_globalnown = _F_ext_global(Eigen::seq(0, numUnknownDisplacements()-1));

        // _assembleStiffnessMatrix();

        std::cout << "\n=== Load Step " << k+1 << " ===" << std::endl;
        _newtonRaphson(_F_ext_global, _d_global);

        // after we've converged, update plastic state for each element (save the latest plastic state for the next time step)
        for (unsigned i = 0; i < _elements.size(); i++)
        {
            _elements[i].updatePlasticState();
        }

        printElementNodalDisplacements(0);
        // break;

        // midterm-specific plotting
        const QuadElement& element = _elements[0];
        const std::vector<double>& integration_points = element.integrationPoints();
        const std::vector<int>& element_global_DOF = element.globalDOF();

        // get element displacement vector
        Eigen::VectorXd U_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
        for (int i = 0; i < NSDIMS*element.numNodes(); i++)
        {
            U_element[i] = _d_global[element_global_DOF[i]];
        }

        double ri = integration_points[1];
        double sj = integration_points[1];
        const Eigen::Matrix2d J_mat = element.UndeformedJacobian(ri, sj);
        const Eigen::Matrix<double, 3, 8> B_mat = element.B(ri, sj, U_element);
        Eigen::Vector3d strain = B_mat * U_element;
        Vector6d strain_6d;
        strain_6d << strain, Eigen::Vector3d::Zero();

        const auto [stress_vec_6d, D_mat_6d, new_plastic_state] = _material->materialSubroutine(strain_6d, element.lastPlasticState(1,1));

        sigma_11.push_back(stress_vec_6d[0]); sigma_22.push_back(stress_vec_6d[1]); sigma_12.push_back(stress_vec_6d[5]);
        alpha.push_back(new_plastic_state.alpha);
        beta_11.push_back(new_plastic_state.beta[0]); beta_22.push_back(new_plastic_state.beta[1]); beta_12.push_back(new_plastic_state.beta[5]);
        C_1111.push_back(D_mat_6d(0,0)); C_2222.push_back(D_mat_6d(1,1)); C_1212.push_back(D_mat_6d(5,5));
        times.push_back(k*0.01);
    }

    // plot using GNU plot
    Gnuplot plt{};

    plt.sendcommand("set terminal wxt size 1500,500"); 
    plt.multiplot(1, 4, "HW6 (a)");
    // plot 1 - sigma vs time
    plt.set_title("sigma vs. time");
    plt.set_xlabel("Time");
    plt.set_ylabel("sigma");
    plt.plot(times, sigma_11, "sigma_{11}");
    plt.plot(times, sigma_22, "sigma_{22}");
    plt.plot(times, sigma_12, "sigma_{12}");
    plt.show();
    // plot 2 - alpha vs time
    plt.set_title("alpha vs time");
    plt.set_xlabel("time");
    plt.set_ylabel("alpha");
    plt.plot(times, alpha, "alpha");
    plt.show();
    // plot 3 - displacement vs time
    plt.set_title("beta vs. time");
    plt.set_xlabel("time");
    plt.set_ylabel("beta");
    plt.plot(times, beta_11, "beta_{11}");
    plt.plot(times, beta_22, "beta_{22}");
    plt.plot(times, beta_12, "beta_{12}");
    plt.show();
    // plot 3 - moduli vs time
    plt.set_title("C_{ep} vs. time");
    plt.set_xlabel("time");
    plt.set_ylabel("C_{ep}");
    plt.plot(times, C_1111, "C_{1111}");
    plt.plot(times, C_2222, "C_{2222}");
    plt.plot(times, C_1212, "C_{1212}");
    plt.show();
}

void Solver::printElementNodalDisplacements(int element_index) const
{
    const QuadElement& element = _elements[element_index];
    const std::vector<int>& element_global_DOF = element.globalDOF();

    // get element displacement vector
    Eigen::VectorXd U_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
    for (int i = 0; i < NSDIMS*element.numNodes(); i++)
    {
        U_element[i] = _d_global[element_global_DOF[i]];
    }

    std::cout << "\n\n=== Element " << element_index << " displacements ===" << std::endl;
    std::cout << "\tx\t\t\ty\t\n----------------\t----------------" << std::endl;
    std::cout << std::setprecision(15);
    for (int i = 0; i < U_element.size(); i+=2)
    {
        std::cout << std::setw(16) << U_element[i] << "\t" << U_element[i+1] << std::endl;
    }

}

void Solver::_newtonRaphson(const Eigen::VectorXd& F_ext, const Eigen::VectorXd& d0)
{
    // compute initial residual
    // _assembleInternalForceVector(d0);
    // Eigen::VectorXd res = F_ext(Eigen::seq(0,numUnknownDisplacements()-1)) - _F_int_global(Eigen::seq(0, numUnknownDisplacements()-1));
    double res0_norm = 0;

    Eigen::VectorXd d = d0;
    Eigen::VectorXd res(numUnknownDisplacements());

    std::cout << "Iter\tRes\n----\t---" << std::endl; 

    for (int i = 0; i < NR_MAX_ITER; i++)
    {
        // recompute internal force at new d
        _assembly(d);

        res = F_ext(Eigen::seq(0,numUnknownDisplacements()-1)) - _F_int_global(Eigen::seq(0, numUnknownDisplacements()-1));
        if (i == 0)
            res0_norm = res.squaredNorm();

        std::cout << i << "\t" << res.squaredNorm() << std::endl;

        // check if we've converged
        if (res.squaredNorm() < res0_norm * NR_TOL)
        {
            break;
        }

        // delta d = K^-1*R
        // we only do this for the unknown displacements, which are the first numUnknownDisplacements() rows in K and the residual
        Eigen::MatrixXd K_unknown = _K_global.block(0,0,numUnknownDisplacements(), numUnknownDisplacements());
        Eigen::VectorXd delta_d = K_unknown.llt().solve(res);
        // std::cout << "\ndelta_d\n" << delta_d << std::endl;
        d(Eigen::seq(0,numUnknownDisplacements()-1)) += delta_d;
    }

    // update displacements once we've converged
    _d_global = d;
}

void Solver::_assembly(const Eigen::VectorXd& d)
{
    _F_int_global = Eigen::VectorXd::Zero(numDOF());
    _K_global = Eigen::MatrixXd::Zero(numDOF(), numDOF());

    // assemble global internal force vector
    for (const auto& element : _elements)
    {
        const std::vector<int>& element_global_DOF = element.globalDOF();
        // get element displacement vector
        Eigen::VectorXd d_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
        for (int i = 0; i < NSDIMS*element.numNodes(); i++)
        {
            d_element[i] = d[element_global_DOF[i]];
        }

        const auto [element_F_int, element_K] = element.elementSubroutine(d_element);
        for (int i = 0; i < 8; i++)
        {
            const int i_global = element_global_DOF[i];
            _F_int_global(i_global) += element_F_int(i);

            for (int j = 0; j < 8; j++)
            {  
                const int j_global = element_global_DOF[j];
                _K_global(i_global, j_global) += element_K(i, j);
            }
        }
    }
}