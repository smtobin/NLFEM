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
    Eigen::VectorXd last_d_global = _d_global;
    _F_ext_global = Eigen::VectorXd::Zero(numDOF());

    // plotting
    std::vector<double> sigma_11, sigma_22, sigma_12, alpha, times;

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
        _newtonRaphson(_F_ext_global, _d_global, last_d_global);

        // after we've converged, update plastic state for each element (save the latest plastic state for the next time step)
        for (unsigned i = 0; i < _elements.size(); i++)
        {
            _elements[i].updatePlasticState();
        }

        printElementNodalDisplacements(0);

        last_d_global = _d_global;

        // plotting
        const QuadElement& element = _elements[0];
        
        Eigen::Vector3d stress_vec3 = element.stressAtIP(1,1);
        PlasticState new_plastic_state = element.lastPlasticState(1,1);
        sigma_11.push_back(stress_vec3[0]); sigma_22.push_back(stress_vec3[1]); sigma_12.push_back(stress_vec3[2]);
        alpha.push_back(new_plastic_state.alpha);
        times.push_back(k*0.1);
    }

    // plot using GNU plot
    Gnuplot plt{};

    plt.sendcommand("set terminal x11 size 1000,500"); 
    plt.multiplot(1, 2, "Plots for Final Project - Uniaxial Tension");
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

void Solver::_newtonRaphson(const Eigen::VectorXd& F_ext, const Eigen::VectorXd& d0, const Eigen::VectorXd& d_old)
{
    
    double res0_norm = 0;

    Eigen::VectorXd d = d0;
    Eigen::VectorXd res(numUnknownDisplacements());

    std::cout << "Iter\tRes\n----\t---" << std::endl; 

    for (int i = 0; i < NR_MAX_ITER; i++)
    {
        // recompoute internal force at new d
        _assembly(d, d_old);

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
        d(Eigen::seq(0,numUnknownDisplacements()-1)) += delta_d;
    }

    // update displacements once we've converged
    _d_global = d;
}

void Solver::_assembly(const Eigen::VectorXd& d_new, const Eigen::VectorXd& d_old)
{
    _F_int_global = Eigen::VectorXd::Zero(numDOF());
    _K_global = Eigen::MatrixXd::Zero(numDOF(), numDOF());

    // assemble global internal force vector
    for (const auto& element : _elements)
    {
        const std::vector<int>& element_global_DOF = element.globalDOF();
        // get element displacement vector
        Eigen::VectorXd d_e_new = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
        Eigen::VectorXd d_e_old = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
        for (int i = 0; i < NSDIMS*element.numNodes(); i++)
        {
            d_e_new[i] = d_new[element_global_DOF[i]];
            d_e_old[i] = d_old[element_global_DOF[i]];
        }

        const auto [element_F_int, element_K] = element.elementSubroutine(d_e_new, d_e_old);
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

// void Solver::_assembleStiffnessMatrix(const Eigen::VectorXd& d)
// {
//     // initialize K to all zeros
//     _K_global = Eigen::MatrixXd::Zero(numDOF(), numDOF());

//     // assemble global stiffness matrix
//     for (const auto& element : _elements)
//     {
//         const std::vector<int>& element_global_DOF = element.globalDOF();
//         // get element displacement vector
//         Eigen::VectorXd d_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
//         for (int i = 0; i < NSDIMS*element.numNodes(); i++)
//         {
//             d_element[i] = d[element_global_DOF[i]];
//         }
        
//         const Eigen::Matrix<double, 8, 8> element_K = element.K(d_element);
//         for (int i = 0; i < 8; i++)
//         {
//             const int i_global = element_global_DOF[i];
//             for (int j = 0; j < 8; j++)
//             {  
//                 const int j_global = element_global_DOF[j];
//                 _K_global(i_global, j_global) += element_K(i, j);
//             }
//         }
//     }
// }

// void Solver::_assembleInternalForceVector(const Eigen::VectorXd& d)
// {
//     // initial _F_int_global to all zeros
//     _F_int_global = Eigen::VectorXd::Zero(numDOF());

//     // assemble global internal force vector
//     for (const auto& element : _elements)
//     {
//         const std::vector<int>& element_global_DOF = element.globalDOF();
//         // get element displacement vector
//         Eigen::VectorXd d_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
//         for (int i = 0; i < NSDIMS*element.numNodes(); i++)
//         {
//             d_element[i] = d[element_global_DOF[i]];
//         }

//         const Eigen::Vector<double,8> element_F_int = element.internalForce(d_element);
//         for (int i = 0; i < 8; i++)
//         {
//             const int i_global = element_global_DOF[i];
//             _F_int_global(i_global) += element_F_int(i);
//         }
//     }
// }