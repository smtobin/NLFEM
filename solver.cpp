#include "solver.hpp"

Solver::Solver( const std::vector<Node> mesh_nodes, const std::vector<ElementNodes> mesh_element_nodes,
                const std::vector<DisplacementBC>& displacement_BCs, const std::vector<ForceBC>& force_BCs)
    : _mesh_nodes(mesh_nodes), _mesh_element_nodes(mesh_element_nodes),
     _displacement_BCs(displacement_BCs), _force_BCs(force_BCs)
     
{
    // create material

    // TODO: better way to switch material depending on hw assignment
    // for Linear FEM code
    _material = std::make_unique<PlaneStrainMaterial>(100, 0.25);

    // for HW 3
    // _material = std::make_unique<HW3Material>(40, -50, -30);


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
                                        _material.get()));
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

void Solver::solve()
{
    // initialize U and R
    _U = Eigen::VectorXd::Zero(numDOF());
    _R = Eigen::VectorXd::Zero(numDOF());

    // put in known displacements and forces
    for (unsigned i = 0; i < _displacement_BCs.size(); i++)
    {
        int global_DOF = _input_to_global_DOF[_displacement_BCs[i].inputDOF()];
        _U[global_DOF] = _displacement_BCs[i].displacement;
    }

    for (unsigned i = 0; i < _force_BCs.size(); i++)
    {
        int global_DOF = _input_to_global_DOF[_force_BCs[i].inputDOF()];
        _R[global_DOF] = _force_BCs[i].force;
    }

    Eigen::VectorXd U_known = _U(Eigen::seq(numUnknownDisplacements(), numDOF()-1));
    Eigen::VectorXd R_known = _R(Eigen::seq(0, numUnknownDisplacements()-1));

    _assembleStiffnessMatrix();

    // _newtonRaphson(_R, _U);

    // partition K into parts that are known and unknown
    Eigen::MatrixXd K_ff = _K.block(0,0,numUnknownDisplacements(), numUnknownDisplacements());
    Eigen::MatrixXd K_fu = _K.block(0,numUnknownDisplacements(), numUnknownDisplacements(), numKnownDisplacements());
    Eigen::MatrixXd K_uu = _K.block(numUnknownDisplacements(), numUnknownDisplacements(), numKnownDisplacements(), numKnownDisplacements());

    // solve for unknown displacements
    _U(Eigen::seq(0, numUnknownDisplacements()-1)) = K_ff.llt().solve(R_known - K_fu*U_known);
    // solve for unknown reaction forces
    _R(Eigen::seq(numUnknownDisplacements(), numDOF()-1)) = K_fu.transpose() * _U(Eigen::seq(0, numUnknownDisplacements()-1)) + K_uu*U_known;
}

void Solver::evaluateElementAtIntegrationPoints(int element_index)
{
    const QuadElement& element = _elements[element_index];
    const std::vector<double>& integration_points = element.integrationPoints();
    const std::vector<int>& element_global_DOF = element.globalDOF();

    // get element displacement vector
    Eigen::VectorXd U_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
    for (int i = 0; i < NSDIMS*element.numNodes(); i++)
    {
        U_element[i] = _U[element_global_DOF[i]];
    }

    // evaluate stress and strain at element integration points
    for (const auto& ri : integration_points)
    {
        for (const auto& sj : integration_points)
        {
            std::cout << "\n === At integration point (r=" << ri << ", s=" << sj << ") ===" << std::endl;
            Eigen::Matrix<double, 3, 8> B = element.B(ri, sj);
            Eigen::Vector3d strain_vec = B*U_element;
            Eigen::Vector3d stress_vec = element.material()->D(strain_vec)*B*U_element;
            std::cout << "  (e_xx, e_yy, e_xy):       " << strain_vec[0] << ", " << strain_vec[1] << ", " << strain_vec[2] << std::endl;
            std::cout << "  (sig_xx, sig_yy, sig_xy): " << stress_vec[0] << ", " << stress_vec[1] << ", " << stress_vec[2] << std::endl;
        }
    }
}

void Solver::_newtonRaphson(const Eigen::VectorXd& F_ext, const Eigen::VectorXd& d0)
{
    // Eigen::VectorXd F_int = Eigen::VectorXd::Zero(F_ext.size());
    // for linear problems, F_int = K*d
    _assembleStiffnessMatrix();
    Eigen::VectorXd F_int = _K * d0;
    Eigen::VectorXd res = F_int - F_ext;
    std::cout << "res0:\n" << res << std::endl;
    double res0_norm = res.norm();

    Eigen::VectorXd d = d0;

    int max_iter = 1000;
    double tol = 1e-6;
    for (int i = 0; i < max_iter; i++)
    {
        if (res.norm() < res0_norm * tol)
            break;

        _assembleStiffnessMatrix();
        Eigen::VectorXd delta_d = _K.inverse() * res;
        std::cout << "delta_d:\n" << delta_d << std::endl;
        d += delta_d;

        res = F_ext - _K*d;
        break;
    }

    _U = d;
    std::cout << "U:\n" << d << std::endl;
}

void Solver::_assembleStiffnessMatrix()
{
    // initialize K to all zeros
    _K = Eigen::MatrixXd::Zero(numDOF(), numDOF());

    // assemble global stiffness matrix
    for (const auto& element : _elements)
    {
        const std::vector<int>& element_global_DOF = element.globalDOF();
        // get element displacement vector
        Eigen::VectorXd d_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
        for (int i = 0; i < NSDIMS*element.numNodes(); i++)
        {
            d_element[i] = _U[element_global_DOF[i]];
        }
        
        const Eigen::Matrix<double, 8, 8> element_K = element.K(d_element);
        for (int i = 0; i < 8; i++)
        {
            const int i_global = element_global_DOF[i];
            for (int j = 0; j < 8; j++)
            {  
                const int j_global = element_global_DOF[j];
                _K(i_global, j_global) += element_K(i, j);
            }
        }
    }
}