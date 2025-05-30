#pragma once

#include "common.hpp"
#include "element.hpp"
#include "material.hpp"

#include <vector>
#include <memory>

/** Represents a prescribed displacement BC */
struct DisplacementBC
{
    int node_index;         // index of the node (in the input mesh)
    Axis axis;              // which axis the displacement BC applies to (X, Y, or Z)
    double displacement;    // the value of the prescribed displacement 

    // constructor for convenience
    DisplacementBC(int node_index_, Axis axis_, double displacement_)
        : node_index(node_index_), axis(axis_), displacement(displacement_)
    {}

    // returns the input DOF numbering
    int inputDOF() const
    {
        return node_index*NSDIMS + axis;
    }
};

/** Represents a prescribed force BC */
struct ForceBC
{
    int node_index;     // index of the node (in the input mesh)
    Axis axis;          // which axis the force BC applies to (X, Y, or Z)
    double force;       // the value of the prescribed force

    // constructor for convenience
    ForceBC(int node_index_, Axis axis_, double force_)
        : node_index(node_index_), axis(axis_), force(force_)
    {}

    // returns the input DOF numbering
    int inputDOF() const
    {
        return node_index*NSDIMS + axis;
    }
};


/** Performs assembly and solves the FEM equations. */
class Solver
{
    public:
    /** Constructor takes in input mesh nodes and elements, and displacement and force boundary conditions. */
    explicit Solver(const std::vector<Node> mesh_nodes, const std::vector<ElementNodes> mesh_element_nodes,
                    const std::vector<DisplacementBC>& displacement_BCs, const std::vector<ForceBC>& force_BCs,
                    const Material* material);
    
    int numKnownDisplacements() const { return _displacement_BCs.size(); }
    int numUnknownDisplacements() const { return numDOF() - numKnownDisplacements(); }
    int numDOF() const { return NSDIMS*_mesh_nodes.size(); }

    const std::vector<QuadElement>& elements() const { return _elements; }

    /** Solve the global finite element equations */
    void solve(int num_load_steps=1);

    /** Print out element nodal displacements. */
    void printElementNodalDisplacements(int element_index) const;

    private:
    /** Helper method to set up the DOF numbering based on the input mesh and BCs */
    void _setupDOF();

    /** Newton-Raphson */
    void _newtonRaphson(const Eigen::VectorXd& F_ext, const Eigen::VectorXd& d0, const Eigen::VectorXd& d_old);

    /** Assembles the global stiffness matrix, updating the class variable _K */
    // void _assembleStiffnessMatrix(const Eigen::VectorXd& d);

    /** Assembles the global internal force vector, updating the class variable _N */
    // void _assembleInternalForceVector(const Eigen::VectorXd& d);

    /** Assembles both the global stiffness matrix and global internal force vector */
    void _assembly(const Eigen::VectorXd& d_new, const Eigen::VectorXd& d_old);

    private:
    /** Stores the nodes in the input mesh. */
    std::vector<Node> _mesh_nodes;

    /** Stores the elements' node indices in the input mesh.
     * Right now, elements are just integer 4-vectors, corresponding to the node
     * indices for a quad element.
     */
    std::vector<ElementNodes> _mesh_element_nodes;

    /** Stores the elements in the mesh. */
    std::vector<QuadElement> _elements;

    /** Stores the imposed displacement boundary conditions. */
    std::vector<DisplacementBC> _displacement_BCs;

    /** Stores the imposed force boundary conditions. */
    std::vector<ForceBC> _force_BCs;

    /** The material for the mesh (right now the same for all elements) */
    const Material* _material;

    /** Maps "input" DOF (i.e. DOF numbered according to the input mesh nodes)
     * to some new global DOF numbering.
     * 
     * For example,
     * The first element in the vector corresponds to the global DOF numbering for the
     * x-axis DOF of the 1st input mesh node.
     * 
     */
    std::vector<int> _input_to_global_DOF;

    /** Global stiffness matrix */
    Eigen::MatrixXd _K_global;

    /** Global displacement vector */
    Eigen::VectorXd _d_global;

    /** Global external force vector */
    Eigen::VectorXd _F_ext_global;

    /** Global internal force vector */
    Eigen::VectorXd _F_int_global;

    /** Material parameters */
    double _density;
    double _E;
    double _mu;

};