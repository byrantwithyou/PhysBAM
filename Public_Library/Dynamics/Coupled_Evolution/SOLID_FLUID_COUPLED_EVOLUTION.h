//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLED_EVOLUTION
//#####################################################################
#ifndef __SOLID_FLUID_COUPLED_EVOLUTION__
#define __SOLID_FLUID_COUPLED_EVOLUTION__

#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Topology/TOPOLOGY_POLICY.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class TV> class FLUIDS_PARAMETERS_UNIFORM;
template<class TV> class FLUID_COLLECTION;
template<class TV> class POISSON_COLLIDABLE_UNIFORM;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV_input>
class SOLID_FLUID_COUPLED_EVOLUTION:public NEWMARK_EVOLUTION<TV_input>
{
    typedef TV_input TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef ARRAY<PAIR<int,T> > FACE_WEIGHT_ELEMENTS;typedef ARRAY<FACE_WEIGHT_ELEMENTS*,FACE_INDEX<TV::m> > T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> > >::TYPE T_ARRAYS_STRUCTURE_SIMPLEX_LIST;
    typedef typename MESH_POLICY<TV::dimension-1>::MESH T_BOUNDARY_MESH;
    typedef typename MESH_POLICY<TV::dimension>::MESH T_MESH;typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension>::OBJECT T_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_BOUNDARY_SIMPLEX;typedef VECTOR<int,TV::dimension> T_BOUNDARY_ELEMENT;typedef VECTOR<int,TV::dimension+1> T_ELEMENT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX T_SIMPLEX;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename TV::SPIN T_SPIN;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension-1>::OBJECT T_THIN_SHELL;typedef typename T_THIN_SHELL::MESH T_THIN_SHELL_MESH;
    typedef VECTOR<int,TV::dimension> T_THIN_SHELL_ELEMENT;typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_THIN_SHELL_SIMPLEX;
    typedef FIELD_PROJECTOR<SPARSE_MATRIX_ENTRY<T>,int,&SPARSE_MATRIX_ENTRY<T>::j> T_PROJECTED_INDEX;
    typedef FIELD_PROJECTOR<SPARSE_MATRIX_ENTRY<T>,T,&SPARSE_MATRIX_ENTRY<T>::a> T_PROJECTED_VALUE;
protected:
    typedef NEWMARK_EVOLUTION<TV> BASE;
    using BASE::solid_body_collection;using BASE::solids_parameters;using BASE::example_forces_and_velocities;
    using BASE::B_full;using BASE::rigid_B_full;using BASE::repulsions;using BASE::rigid_deformable_collisions;using BASE::Initialize_World_Space_Masses;
    using BASE::world_space_rigid_mass_inverse;using BASE::world_space_rigid_mass;using BASE::solids_evolution_callbacks;using BASE::krylov_vectors;
    using BASE::X_save;using BASE::rigid_frame_save;using BASE::V_save;using BASE::rigid_velocity_save;using BASE::rigid_angular_momentum_save;

    static const int rows_per_rigid_body=TV::dimension+T_SPIN::dimension;

public:
    int rigid_body_count;
    bool print_matrix_rhs_and_solution;
protected:
    ARRAY<int> kinematic_rigid_bodies;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies;
    T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS dual_cell_weights;
    T_FACE_ARRAYS_FACE_WEIGHT_ELEMENTS rigid_body_dual_cell_weights;
    T_FACE_ARRAYS_SCALAR dual_cell_fluid_volume;
    T_FACE_ARRAYS_BOOL dual_cell_contains_solid;
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters;
    FLUID_COLLECTION<TV>& fluid_collection;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > nodal_fluid_mass;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > rigid_body_fluid_mass;
    ARRAY<TV> rigid_body_updated_center_of_mass;
    ARRAY<SYMMETRIC_MATRIX<T,TV::SPIN::m> > rigid_body_fluid_inertia;
    ARRAY<TV> ar_full,z_full,zaq_full;ARRAY<TWIST<TV> > rigid_ar_full,rigid_z_full,rigid_zaq_full; // extra vectors for conjugate residual
    T_FACE_ARRAYS_SCALAR solid_projected_face_velocities_star;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> coupled_vectors;

    POISSON_COLLIDABLE_UNIFORM<TV>* Get_Poisson()
    {return (fluids_parameters.compressible ? fluids_parameters.euler->euler_projection.elliptic_solver : fluids_parameters.incompressible->projection.poisson_collidable);}

    T_FACE_ARRAYS_SCALAR& Get_Face_Velocities()
    {return fluids_parameters.compressible ? fluids_parameters.euler->euler_projection.face_velocities:fluid_collection.incompressible_fluid_collection.face_velocities;}

    const GRID<TV>& Get_Grid()
    {return (fluids_parameters.compressible ? fluids_parameters.euler->grid : fluids_parameters.incompressible->projection.p_grid);}

public:
    
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > J_deformable;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > J_rigid;
    ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > J_rigid_kinematic;
    ARRAY<ARRAY<TV_INT> > matrix_index_to_cell_index_array;

    T_ARRAYS_SCALAR& Get_Pressure()
    {return (fluids_parameters.compressible ? fluids_parameters.euler->euler_projection.p : fluids_parameters.incompressible->projection.p);}
    
    bool Negative(const GRID<TV>& grid,const int axis,const TV_INT& face_index,const T_ARRAYS_SCALAR& phi)
    {// this check is potentially expensive! leave inefficient for the moment
    TV_INT cells[GRID<TV>::number_of_cells_per_node];
    for(int node=0;node<GRID<TV>::number_of_nodes_per_face;node++){
        grid.Cells_Neighboring_Node(grid.Face_Node_Index(axis,face_index,node),cells);
        for(int c=0;c<GRID<TV>::number_of_cells_per_node;c++) if(phi(cells[c])<=0) return true;}
    return false;}

    bool Simulate_Fluids() const
    {return (solids_fluids_parameters.mpi_solid_fluid || fluids_parameters.simulate) && (fluids_parameters.smoke || fluids_parameters.fire || fluids_parameters.water || fluids_parameters.sph || fluids_parameters.compressible);}

    bool Simulate_Solids() const
    {return ((solids_fluids_parameters.mpi_solid_fluid || solid_body_collection.deformable_body_collection.simulate) && solid_body_collection.deformable_body_collection.particles.Size()) || ((solids_fluids_parameters.mpi_solid_fluid || solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies) && solid_body_collection.rigid_body_collection.rigid_body_particles.Size());}

    SOLID_FLUID_COUPLED_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,
        EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities_input,FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_input,
        FLUID_COLLECTION<TV>& fluid_collection_input,SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters);
    virtual ~SOLID_FLUID_COUPLED_EVOLUTION();

//#####################################################################
    T Get_Density_At_Face(const int axis,const TV_INT& face_index);
    void Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update) PHYSBAM_OVERRIDE;
    void Transfer_Momentum_And_Set_Boundary_Conditions(const T time,GENERALIZED_VELOCITY<TV>* B=0);
    void Set_Dirichlet_Boundary_Conditions(const T time);
    void Compute_W(const T current_position_time);
    void Compute_Coupling_Terms_Deformable(const T_ARRAYS_INT& cell_index_to_matrix_index,const ARRAY<INTERVAL<int> >& interior_regions,const int colors);
    void Compute_Coupling_Terms_Rigid(const T_ARRAYS_INT& cell_index_to_matrix_index,const ARRAY<INTERVAL<int> >& interior_regions,const int colors);
    void Add_Nondynamic_Solids_To_Right_Hand_Side(ARRAY<ARRAY<T> >& right_hand_side,const ARRAY<INTERVAL<int> >& interior_regions,const int colors);
    void Apply_Pressure(const T dt,const T time);
    void Average_Solid_Projected_Face_Velocities_For_Energy_Update(const T_FACE_ARRAYS_SCALAR& solid_projected_face_velocities_star,const T_FACE_ARRAYS_SCALAR& solid_projected_face_velocities_np1,T_FACE_ARRAYS_SCALAR& face_velocities);
    void Apply_Solid_Boundary_Conditions(const T time,const bool use_pseudo_velocities,T_FACE_ARRAYS_SCALAR& face_velocities);
    void Set_Neumann_and_Dirichlet_Boundary_Conditions(const T time);
//#####################################################################
};
}
#endif
