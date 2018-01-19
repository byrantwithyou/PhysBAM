//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES
//#####################################################################
#ifndef __IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES__
#define __IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES__

#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Forces/COLLISION_FORCE.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES:public COLLISION_FORCE<TV>
{
    typedef typename TV::SCALAR T;
public:
    using DEFORMABLES_FORCES<TV>::particles;using COLLISION_FORCE<TV>::coefficient_of_friction;
    IMPLICIT_OBJECT<TV>* implicit_object;
    bool own_implicit_object;
    ARRAY<int> colliding_particles;

    T stiffness;
    T separation_parameter;
    T length_scale;

    ARRAY<int> penetrating_particles;
    T pe;
    ARRAY<TV> grad_pe;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > H_pe;

    IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES(DEFORMABLE_PARTICLES<TV>& particles,IMPLICIT_OBJECT<TV>* implicit_object,
        T stiffness=(T)1e4,T separation_parameter=(T)1e-4,T length_scale=1);
    virtual ~IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES();

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    int Velocity_Dependent_Forces_Size() const override;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const override;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA> frequency) override;
    T Potential_Energy(const T time) const override;
    void Update_Position_Based_State_Particle(int p);
    void Apply_Friction(ARRAY_VIEW<TV> V,const T time,const T dt) const override;
//#####################################################################
};
}
#endif
