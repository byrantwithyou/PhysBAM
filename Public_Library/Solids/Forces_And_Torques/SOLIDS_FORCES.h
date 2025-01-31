//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FORCES
//#####################################################################
#ifndef __SOLIDS_FORCES__
#define __SOLIDS_FORCES__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class GENERALIZED_VELOCITY;
class SEGMENT_MESH;


template<class TV>
class SOLIDS_FORCES
{
    typedef typename TV::SCALAR T;
public:
    typedef typename RIGIDS_FORCES<TV>::FREQUENCY_DATA RIGID_FREQUENCY_DATA;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;
    DEFORMABLE_PARTICLES<TV>& particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
protected:
    T cfl_number;
    bool allow_external_cfl_number;
public:
    bool cfl_initialized;

    bool use_rest_state_for_strain_rate;
    bool limit_time_step_by_strain_rate;
    T max_strain_per_time_step; // for limiting the timestep in the CFL calculation
    int unique_id;
    bool compute_half_forces;

    SOLIDS_FORCES(DEFORMABLE_PARTICLES<TV>& particles,RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    SOLIDS_FORCES(const SOLIDS_FORCES&) = delete;
    void operator=(const SOLIDS_FORCES&) = delete;
    virtual ~SOLIDS_FORCES();

    static int Get_Unique_Id()
    {static int next_unique_id=0;return ++next_unique_id;}

    void Set_CFL_Number(const T cfl_number_input)
    {if(allow_external_cfl_number) cfl_number=cfl_number_input;}

    bool CFL_Valid() const
    {return cfl_initialized;}

    void Invalidate_CFL()
    {cfl_initialized=false;}

    void Validate_CFL()
    {cfl_initialized=true;}

//#####################################################################
    virtual void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true);
    virtual void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1);
    virtual void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const=0;
    virtual void Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)=0;
    virtual void Update_Position_Based_State(const T time);
    virtual void Add_Velocity_Independent_Forces(GENERALIZED_VELOCITY<TV>& F,const T time) const=0;
    virtual void Add_Velocity_Dependent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time) const=0;
    virtual int Velocity_Dependent_Forces_Size() const;
    virtual void Add_Velocity_Dependent_Forces_First_Half(const GENERALIZED_VELOCITY<TV>& V,ARRAY_VIEW<T> aggregate,const T time) const;
    virtual void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,GENERALIZED_VELOCITY<TV>& F,const T time) const;
    virtual void Add_Implicit_Velocity_Independent_Forces(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,const T time,bool transpose=false) const;
    virtual void Enforce_Definiteness(const bool enforce_definiteness_input);
    virtual T CFL_Strain_Rate() const=0;
    virtual void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<RIGID_FREQUENCY_DATA> rigid_frequency)=0;
    virtual T Potential_Energy(const T time) const;
    virtual void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const;
//#####################################################################
};
}
#endif
