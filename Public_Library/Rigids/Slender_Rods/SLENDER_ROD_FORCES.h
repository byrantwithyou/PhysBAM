//#####################################################################
// Copyright 2002-2009, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SLENDER_ROD_FORCES
//#####################################################################
#ifndef __SLENDER_ROD_FORCES__
#define __SLENDER_ROD_FORCES__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class TV> class RIGID_SLENDER_ROD_COLLECTION;
class SEGMENT_MESH;

template<class TV>
class SLENDER_ROD_FORCES
{
    typedef typename TV::SCALAR T;
public:
    RIGID_SLENDER_ROD_COLLECTION<TV>& rigid_slender_rod_collection;
protected:
    T cfl_number;
    bool allow_external_cfl_number;
public:
    struct FREQUENCY_DATA
    {
        FREQUENCY_DATA()
            :elastic_squared(0),damping(0)
        {}

        T elastic_squared,damping;
    };
    bool cfl_initialized;

    bool use_rest_state_for_strain_rate;
    bool limit_time_step_by_strain_rate;
    T max_strain_per_time_step; // for limiting the timestep in the CFL calculation
    int unique_id;
    bool compute_half_forces;

    SLENDER_ROD_FORCES(RIGID_SLENDER_ROD_COLLECTION<TV>& rigid_slender_rod_collection);
    SLENDER_ROD_FORCES(const SLENDER_ROD_FORCES&) = delete;
    void operator=(const SLENDER_ROD_FORCES&) = delete;
    virtual ~SLENDER_ROD_FORCES();

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
    virtual void Update_Mpi(const ARRAY<bool>& particle_is_simulated)=0;
    virtual void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true);
    virtual void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1);
    virtual void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const=0;
    virtual void Update_Position_Based_State(const T time){}
    virtual void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const=0;
    virtual void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const=0;
    virtual int Velocity_Dependent_Forces_Size() const;
    virtual void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const;
    virtual void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const;
    virtual void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time,bool transpose=false) const;
    virtual void Enforce_Definiteness(const bool enforce_definiteness_input);
    virtual T CFL_Strain_Rate() const=0;
    virtual void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> rigid_frequency)=0;
    virtual T Potential_Energy(const T time) const;
//#####################################################################
};
}
#endif
