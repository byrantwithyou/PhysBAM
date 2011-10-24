//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_AREA_PENALTY_FORCE
//#####################################################################
#ifndef __COLLISION_AREA_PENALTY_FORCE__
#define __COLLISION_AREA_PENALTY_FORCE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class VOLUME_COLLISIONS;
template<class T> class TRIANGULATED_AREA;
template<class TV>
class COLLISION_AREA_PENALTY_FORCE:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    VOLUME_COLLISIONS<TV>& volume_collisions;
    T force_coefficient;

    COLLISION_AREA_PENALTY_FORCE(PARTICLES<TV>& particles);

    virtual ~COLLISION_AREA_PENALTY_FORCE();

//#####################################################################
    void Add_Mesh(TRIANGULATED_AREA<T>& ta);
    virtual void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true);
    virtual void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1);
    virtual void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const;
    virtual void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids);
    virtual void Update_Position_Based_State(const T time,const bool is_position_update);
    virtual void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const;
    virtual void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const;
    virtual int Velocity_Dependent_Forces_Size() const;
    virtual void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const;
    virtual void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const;
    virtual void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const;
    virtual void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const;
    virtual void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const;
    virtual void Enforce_Definiteness(const bool enforce_definiteness_input);
    virtual T CFL_Strain_Rate() const;
    virtual void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency);
    virtual T Potential_Energy(const T time) const;
    virtual void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const;
//#####################################################################
};
}
#endif
