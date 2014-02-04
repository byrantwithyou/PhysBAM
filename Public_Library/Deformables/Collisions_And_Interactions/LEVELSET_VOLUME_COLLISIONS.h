//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_VOLUME_COLLISIONS
//#####################################################################
#ifndef __LEVELSET_VOLUME_COLLISIONS__
#define __LEVELSET_VOLUME_COLLISIONS__

#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Forces/COLLISION_FORCE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class LEVELSET_VOLUME_COLLISIONS:public COLLISION_FORCE<TV>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,TV::m+1> SIMPLEX_NODES;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE SIMPLEX_FACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX SIMPLEX;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT OBJECT;
    typedef INDIRECT_ARRAY<ARRAY_VIEW<TV,int>,VECTOR<int,2*TV::m+2>& > X_ARRAY;
    typedef INDIRECT_ARRAY<ARRAY<T>,VECTOR<int,TV::m+1>& > PHI_ARRAY;
    struct HYPER_PLANE
    {
        unsigned char plane;
        TV normal;
        T s;
        T Signed_Distance(const TV& location) const
        {
            return normal.Dot(location)-s;
        }
    };

public:
    using DEFORMABLES_FORCES<TV>::particles;using COLLISION_FORCE<TV>::coefficient_of_friction;
    OBJECT& collision_body;
    ARRAY<T> undeformed_phi;

    T stiffness;

    ARRAY<VECTOR<int,2*TV::m+2> > overlapping_particles;
    T pe;
    ARRAY<VECTOR<TV,2*TV::m+2> > grad_pe;
    ARRAY<MATRIX<MATRIX<T,TV::m>,2*TV::m+2> > H_pe;

    LEVELSET_VOLUME_COLLISIONS(DEFORMABLE_PARTICLES<TV>& particles,OBJECT& collision_body,IMPLICIT_OBJECT<TV>& implicit_surface,
        T stiffness=(T)1e4);
    virtual ~LEVELSET_VOLUME_COLLISIONS();

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
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
    void Apply_Friction(ARRAY_VIEW<TV> V,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
    void Simplex_Intersection(const SIMPLEX& s,const ARRAY<HYPER_PLANE>& f,ARRAY<unsigned char>& polytope);
    void Integrate_Simplex(VECTOR<unsigned char,TV::m+1> simplex,const X_ARRAY& X,const PHI_ARRAY& nodewise_undeformed_phi,VECTOR<TV,2*TV::m+2>& df,MATRIX<MATRIX<T,TV::m>,2*TV::m+2>& ddf);
};
}
#endif
