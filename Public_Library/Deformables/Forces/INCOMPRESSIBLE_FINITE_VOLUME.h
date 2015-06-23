//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FINITE_VOLUME
//#####################################################################
#ifndef __INCOMPRESSIBLE_FINITE_VOLUME__
#define __INCOMPRESSIBLE_FINITE_VOLUME__

#include <Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Math_Tools/FACTORIAL.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Topology/TOPOLOGY_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLES_COLLISIONS_FORWARD.h>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;
template<class TV> class COLLISION_PARTICLE_STATE;
template<class TV> struct PRECOMPUTE_PROJECT;
namespace{template<class TV,int d> class POISSON_SYSTEM;}

template<class TV>
class INCOMPRESSIBLE_FINITE_VOLUME_BASE
{
    typedef typename TV::SCALAR T;
public:
    virtual ~INCOMPRESSIBLE_FINITE_VOLUME_BASE(){}
    virtual void Make_Incompressible(const T dt,const bool correct_volume)=0;
    virtual void Test_System()=0;
    virtual void Set_Neumann_Boundary_Conditions(const ARRAY<COLLISION_PARTICLE_STATE<TV> >* particle_states,TRIANGLE_REPULSIONS<TV>* repulsions)=0;
};

template<class TV_input,int d>
class INCOMPRESSIBLE_FINITE_VOLUME:public DEFORMABLES_FORCES<TV_input>,public INCOMPRESSIBLE_FINITE_VOLUME_BASE<TV_input>
{
    typedef TV_input TV;
    typedef typename TV::SCALAR T;
    typedef MATRIX<T,TV::m,d> T_MATRIX;
    typedef typename MESH_POLICY<d-1>::MESH T_BOUNDARY_MESH;
    typedef typename MESH_POLICY<d>::MESH T_MESH;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_OBJECT;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;

    typedef DEFORMABLES_FORCES<TV> BASE;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;
    typedef KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<T> > > KRYLOV_VECTOR_T;
public:
    using BASE::particles;

    STRAIN_MEASURE<TV,d>& strain_measure;
    bool disable_projection;
    ARRAY<T> boundary_pressures; // outside pressures per node (non-boundary nodes ignored); empty array for zero boundary pressure
    T minimum_volume_recovery_time_scale;
    int max_cg_iterations;
    MPI_SOLIDS<TV>* mpi_solids;
    bool merge_at_boundary;
    bool use_neumann;
    bool use_self_moving_projection;
    bool use_rigid_clamp_projection;
    bool use_diagonal_preconditioner;
    // TODO: make the following private (unless they disappear)
    ARRAY<int> neumann_boundary_count;
    TRIANGLE_REPULSIONS<TV>* repulsions;

    struct PROJECTION_DATA
    {
        ARRAY<TV> neumann_boundary_normals;
        ARRAY<int> neumann_boundary_nodes;
        ARRAY<int> neumann_boundary_nodes_isolated;
        ARRAY<int> fixed_nodes;
        ARRAY<REPULSION_PAIR<TV> > point_face_pairs;
        ARRAY<REPULSION_PAIR<TV> > edge_edge_pairs;
        ARRAY<PRECOMPUTE_PROJECT<TV> > point_face_precomputed;
        ARRAY<PRECOMPUTE_PROJECT<TV> > edge_edge_precomputed;
    };
    PROJECTION_DATA projection_data;
    ARRAY<ARRAY<int> > node_regions;
    mutable T project_amount;
private:
    ARRAY<VECTOR<int,2> > boundary_to_element; // for each boundary element, which interior element it's in, and which Bs column to use
    ARRAY<T> rest_volumes_full,volumes_full;
    mutable ARRAY<T> saved_volumes_full;
    T total_rest_volume,total_volume;
    ARRAY<KRYLOV_VECTOR_BASE<T>*> cg_vectors;
    ARRAY<T> pressure_full,divergence_full;
    ARRAY<T> diagonal_preconditioner_full;
    mutable ARRAY<TV> gradient_full;
    ARRAY<T_MATRIX> Bs_per_node;
    ARRAY<TV> boundary_normals;
    FORCE_ELEMENTS force_elements;
    FORCE_ELEMENTS force_boundary_elements;
    FORCE_ELEMENTS force_dynamic_particles;
    ARRAY<int> force_dynamic_particles_list;
    friend class POISSON_SYSTEM<TV,d>;
public:

    INCOMPRESSIBLE_FINITE_VOLUME(STRAIN_MEASURE<TV,d>& strain_measure);
    virtual ~INCOMPRESSIBLE_FINITE_VOLUME();

    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override {}
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override {}
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override {}
    T CFL_Strain_Rate() const override {return FLT_MAX;}
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override {}

//#####################################################################
    static INCOMPRESSIBLE_FINITE_VOLUME* Create(T_OBJECT& object,const bool verbose);
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Make_Incompressible(const T dt,const bool correct_volume) override;
    void Test_System() override;
    T Max_Relative_Velocity_Error();
    void Save_Volumes();
    void Check_Improvement();
private:
    void Gradient(ARRAY_VIEW<const T> pressure,ARRAY<TV>& gradient) const;
    void Negative_Divergence(ARRAY_VIEW<const TV> V,ARRAY<T>& divergence) const;
    void Project_All_Clamping_Constraints(ARRAY_VIEW<TV> field,const PROJECTION_DATA& data) const;
    void Project_All_Isolated_Clamping_Constraints(ARRAY_VIEW<TV> field,const PROJECTION_DATA& data) const;
    void Project_Vector_Field(ARRAY_VIEW<TV> field) const;
    void Set_Neumann_Boundary_Conditions(const ARRAY<COLLISION_PARTICLE_STATE<TV> >* particle_states,TRIANGLE_REPULSIONS<TV>* repulsions) override;
    void Diagonal_Elements(ARRAY<T>& D) const;
    void Update_Preconditioner();
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
//#####################################################################
};

template<class T_OBJECT> INCOMPRESSIBLE_FINITE_VOLUME<typename T_OBJECT::VECTOR_T,T_OBJECT::MESH::dimension>*
Create_Incompressible_Finite_Volume(T_OBJECT& object,const bool verbose=true)
{
    typedef typename T_OBJECT::VECTOR_T TV;typedef typename TV::SCALAR T;static const int d=T_OBJECT::MESH::dimension;
    INCOMPRESSIBLE_FINITE_VOLUME<TV,d>* fvm=INCOMPRESSIBLE_FINITE_VOLUME<TV,d>::Create(object,verbose);
    fvm->Limit_Time_Step_By_Strain_Rate(true,(T).1);fvm->Use_Rest_State_For_Strain_Rate(true);
    return fvm;
}
}
#endif
