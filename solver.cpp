#include "solver.hpp"

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

        _newtonRaphson(_F_ext_global, _d_global);
    }


    // Test out deformation gradient

    const QuadElement& element = _elements[0];
    const std::vector<double>& integration_points = element.integrationPoints();
    const std::vector<int>& element_global_DOF = element.globalDOF();

    // get element displacement vector
    Eigen::VectorXd U_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
    for (int i = 0; i < NSDIMS*element.numNodes(); i++)
    {
        U_element[i] = _d_global[element_global_DOF[i]];
    }
    
    Eigen::Matrix2d F = element.deformationGradient(0,0, U_element);
    std::cout << "Deformation gradient:\n" << F << std::endl;
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
        U_element[i] = _d_global[element_global_DOF[i]];
    }

    // evaluate stress and strain at element integration points
    for (const auto& ri : integration_points)
    {
        for (const auto& sj : integration_points)
        {
            std::cout << std::setprecision(5);
            std::cout << "\n === At integration point (r=" << ri << ", s=" << sj << ") ===" << std::endl;
            // find deformation gradient at (r,s) given the current deformation
            const Eigen::Matrix2d F_mat = element.deformationGradient(ri, sj, U_element);
            const auto [stress_vec, D_mat] = element.material()->materialSubroutine(F_mat);

            // calculate Green strain
            const Eigen::Matrix2d E_mat = 0.5*(F_mat.transpose()*F_mat - Eigen::Matrix2d::Identity());
            const Eigen::Vector3d E_vec(E_mat(0,0), E_mat(1,1), 0.5*E_mat(0,1));
            std::cout << std::setprecision(15);
            std::cout << "  (e_xx, e_yy, e_xy):       " << E_vec[0] << ", " << E_vec[1] << ", " << E_vec[2] << std::endl;
            std::cout << "  (sig_xx, sig_yy, sig_xy): " << stress_vec[0] << ", " << stress_vec[1] << ", " << stress_vec[2] << std::endl;
        }
    }
}

void Solver::_newtonRaphson(const Eigen::VectorXd& F_ext, const Eigen::VectorXd& d0)
{
    // compute initial residual
    _assembleInternalForceVector(d0);
    Eigen::VectorXd res = F_ext(Eigen::seq(0,numUnknownDisplacements()-1)) - _F_int_global(Eigen::seq(0, numUnknownDisplacements()-1));
    double res0_norm = res.squaredNorm();

    Eigen::VectorXd d = d0;

    for (int i = 0; i < NR_MAX_ITER; i++)
    {
        // check if we've converged
        if (res.squaredNorm() < res0_norm * NR_TOL)
        {
            break;
        }

        _assembleStiffnessMatrix(d);

        // delta d = K^-1*R
        // we only do this for the unknown displacements, which are the first numUnknownDisplacements() rows in K and the residual
        Eigen::MatrixXd K_unknown = _K_global.block(0,0,numUnknownDisplacements(), numUnknownDisplacements());
        Eigen::VectorXd delta_d = K_unknown.llt().solve(res);
        d(Eigen::seq(0,numUnknownDisplacements()-1)) += delta_d;

        // recompute internal force at new d
        _assembleInternalForceVector(d);

        res = F_ext(Eigen::seq(0,numUnknownDisplacements()-1)) - _F_int_global(Eigen::seq(0, numUnknownDisplacements()-1));
    }

    // update displacements once we've converged
    _d_global = d;
}

void Solver::_assembleStiffnessMatrix(const Eigen::VectorXd& d)
{
    // initialize K to all zeros
    _K_global = Eigen::MatrixXd::Zero(numDOF(), numDOF());

    // assemble global stiffness matrix
    for (const auto& element : _elements)
    {
        const std::vector<int>& element_global_DOF = element.globalDOF();
        // get element displacement vector
        Eigen::VectorXd d_element = Eigen::VectorXd::Zero(NSDIMS*element.numNodes());
        for (int i = 0; i < NSDIMS*element.numNodes(); i++)
        {
            d_element[i] = d[element_global_DOF[i]];
        }
        
        const Eigen::Matrix<double, 8, 8> element_K = element.K(d_element);
        for (int i = 0; i < 8; i++)
        {
            const int i_global = element_global_DOF[i];
            for (int j = 0; j < 8; j++)
            {  
                const int j_global = element_global_DOF[j];
                _K_global(i_global, j_global) += element_K(i, j);
            }
        }
    }
}

void Solver::_assembleInternalForceVector(const Eigen::VectorXd& d)
{
    // initial _F_int_global to all zeros
    _F_int_global = Eigen::VectorXd::Zero(numDOF());

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

        const Eigen::Vector<double,8> element_F_int = element.internalForce(d_element);
        for (int i = 0; i < 8; i++)
        {
            const int i_global = element_global_DOF[i];
            _F_int_global(i_global) += element_F_int(i);
        }
    }
}