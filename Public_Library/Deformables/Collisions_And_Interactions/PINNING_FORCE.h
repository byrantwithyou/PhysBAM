//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PINNING_FORCE
//#####################################################################
#ifndef __PINNING_FORCE__
#define __PINNING_FORCE__

#include <Core/Arrays/ARRAY.h>
#include <Core/Vectors/VECTOR.h>
#include <Deformables/Forces/LAGGED_FORCE.h>
#include <functional>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;

template<class TV>
class PINNING_FORCE:public LAGGED_FORCE<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename LAGGED_FORCE<TV>::FREQUENCY_DATA FREQUENCY_DATA;
public:
    using DEFORMABLES_FORCES<TV>::particles;

    T stiffness;
    T damping_coefficient;
    T& dt;

    ARRAY<int> pinned_particles;
    ARRAY<std::function<TV(T)> > targets;

    void Add_Target(int p,std::function<TV(T)> func)
    {pinned_particles.Append(p);targets.Append(func);}

    T E,H_E;
    ARRAY<TV> grad_E;
    ARRAY<TV> X0;

    PINNING_FORCE(DEFORMABLE_PARTICLES<TV>& particles,T& dt,T stiffness,T damping_coefficient);
    virtual ~PINNING_FORCE();

    void Lagged_Update_Position_Based_State(const T time) override;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const override;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) override;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) override;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const override;
    void Enforce_Definiteness(const bool enforce_definiteness_input) override;
    T CFL_Strain_Rate() const override;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) override;
    T Potential_Energy(const T time) const override;
//#####################################################################
};
}
#endif
