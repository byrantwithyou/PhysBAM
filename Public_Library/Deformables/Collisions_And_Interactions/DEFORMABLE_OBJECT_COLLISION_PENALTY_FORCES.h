//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES
//#####################################################################
#ifndef __DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES__
#define __DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES__

#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Deformables/Forces/COLLISION_FORCE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGULATED_SURFACE;
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES:public COLLISION_FORCE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    using DEFORMABLES_FORCES<TV>::particles;using COLLISION_FORCE<TV>::coefficient_of_friction;
    DEFORMABLE_PARTICLES<TV> &undeformed_particles;
    T_OBJECT& collision_body;
    T_SURFACE& undeformed_triangulated_surface;
    T_SURFACE& triangulated_surface;
    IMPLICIT_OBJECT<TV>& implicit_surface;
    ARRAY<int> colliding_particles;
    ARRAY<int> closest_surface_triangle;
    ARRAY<ARRAY<VECTOR<int,TV::m> > > extra_surface_triangles;

    T stiffness;
    T separation_parameter;

    ARRAY<VECTOR<int,TV::m+1> > penetrating_particles;
    T pe;
    ARRAY<VECTOR<TV,TV::m+1> > grad_pe;
    ARRAY<VECTOR<VECTOR<MATRIX<T,TV::m>,TV::m+1>,TV::m+1> > H_pe;
    ARRAY<VECTOR<T,TV::m+1> > stored_weights;

    DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles,
        DEFORMABLE_PARTICLES<TV>& undeformed_particles,T_OBJECT& collision_body,
        T_SURFACE& undeformed_triangulated_surface,IMPLICIT_OBJECT<TV>& implicit_surface,
        T stiffness=(T)1e4,T separation_parameter=(T)1e-4);
    virtual ~DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES();

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T scale,const T time) const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Update_Position_Based_State_Particle(int p);
    void Update_Penetrating_Particles(int p);
    void Update_Surface_Triangles();
    void Penalty(VECTOR<int,4> nodes,const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,4>&>& X,T& e,VECTOR<TV,4>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he);
    void Penalty(VECTOR<int,3> nodes,const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,4>&>& X,T& e,VECTOR<TV,4>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he);
    void Penalty(VECTOR<int,2> nodes,const INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,4>&>& X,T& e,VECTOR<TV,4>& de,VECTOR<VECTOR<MATRIX<T,TV::m>,4>,4>& he);
    int Estimate_Closest_Undeformed_Surface_Triangle(const TV& X,int p) const;
    void Apply_Friction(ARRAY_VIEW<TV> V,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
