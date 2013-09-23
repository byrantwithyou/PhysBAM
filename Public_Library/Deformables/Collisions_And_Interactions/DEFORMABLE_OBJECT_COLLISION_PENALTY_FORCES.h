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
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;

template<class TV>
class DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    typedef DEFORMABLES_FORCES<TV> BASE;
public:
    using BASE::particles;
    TETRAHEDRALIZED_VOLUME<T>* collision_body;
    ARRAY<T> undeformed_phi;
    ARRAY<int> colliding_particles;

    T stiffness;
    T separation_parameter;
    T length_scale;

    ARRAY<VECTOR<int,5>> penetrating_particles;
    T pe;
    ARRAY<VECTOR<TV,5>> grad_pe;
    ARRAY<VECTOR<VECTOR<MATRIX<T,TV::m>,5>,5> > H_pe;

    DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles,TETRAHEDRALIZED_VOLUME<T>* collision_body,ARRAY<T>& undeformed_phi,
        T stiffness=(T)1e4,T separation_parameter=(T)1e-4,T length_scale=1);
    virtual ~DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES();

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Update_Position_Based_State_Particle(int p);
//#####################################################################
};
}
#endif
