//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RALEIGH_DAMPING_FORCE
//#####################################################################
#ifndef __RALEIGH_DAMPING_FORCE__
#define __RALEIGH_DAMPING_FORCE__

#include <Deformables/Forces/LAGGED_FORCE.h>
namespace PhysBAM{

template<class TV>
class RALEIGH_DAMPING_FORCE:public LAGGED_FORCE<TV>
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    using LAGGED_FORCE<TV>::particles;using typename LAGGED_FORCE<TV>::FREQUENCY_DATA;
    DEFORMABLES_FORCES<TV>& force;
    T coefficient,dt_dv_over_dx;
    T& dt;
    ARRAY<TV> D_V0;
    mutable ARRAY<TV> tmp;

    RALEIGH_DAMPING_FORCE(DEFORMABLE_PARTICLES<TV>& particles_input,DEFORMABLES_FORCES<TV>* pforce,T coefficient,T dt_dv_over_dx,T& dt)
        :LAGGED_FORCE<TV>(particles_input),force(*pforce),coefficient(coefficient),dt_dv_over_dx(dt_dv_over_dx),dt(dt)
    {force.Enforce_Definiteness(true);}

    virtual ~RALEIGH_DAMPING_FORCE()
    {}

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Lagged_Update_Position_Based_State(const T time) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T Potential_Energy(const T time) const override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override;
//#####################################################################
};
}
#endif
