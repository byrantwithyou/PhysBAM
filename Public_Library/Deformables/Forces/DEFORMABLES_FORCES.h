//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_FORCES
//#####################################################################
#ifndef __DEFORMABLES_FORCES__
#define __DEFORMABLES_FORCES__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;
template<class TV> class DEFORMABLE_PARTICLES;
class SEGMENT_MESH;

template<class TV>
struct FORCE_DATA{
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef typename TV::SCALAR T;

    FORCE_DATA():state(0)
    {}

    std::string name;
    T state; // holds signed strech for springs
    TV first_action_point,second_action_point;

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,name,state,first_action_point,second_action_point);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary<RW>(output,name,state,first_action_point,second_action_point);}
};

template<class TV>
class DEFORMABLES_FORCES
{
    typedef typename TV::SCALAR T;
public:
    DEFORMABLE_PARTICLES<TV>& particles;
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

    DEFORMABLES_FORCES(DEFORMABLE_PARTICLES<TV>& particles);
    DEFORMABLES_FORCES(const DEFORMABLES_FORCES&) = delete;
    void operator=(const DEFORMABLES_FORCES&) = delete;
    virtual ~DEFORMABLES_FORCES();

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
    virtual void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)=0;
    virtual void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian);
    virtual void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const=0;
    virtual void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const=0;
    virtual int Velocity_Dependent_Forces_Size() const;
    virtual void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const;
    virtual void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const;
    virtual void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const;
    virtual void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time,bool transpose=false) const=0;
    virtual void Enforce_Definiteness(const bool enforce_definiteness_input);
    virtual T CFL_Strain_Rate() const=0;
    virtual void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)=0;
    virtual T Potential_Energy(const T time) const;
    virtual void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const;
    void Test_Diff(const T time);
//#####################################################################
};
// defined in RIGIDS_FORCES
template<class T_ARRAY> void Update_Force_Particles(ARRAY<int>& force_particles,
    const ARRAY_BASE<int,T_ARRAY>& particles,
    const ARRAY<bool>& particle_is_simulated,bool check_dups);
template<int d,class T_ARRAY> void Update_Force_Elements(ARRAY<int>& force_elements,
    const ARRAY_BASE<VECTOR<int,d>,T_ARRAY>& elements,
    const ARRAY<bool>& particle_is_simulated);

// tag classes to allow iteration over certain types of forces
struct SPRINGS_TAG{};
struct FINITE_VOLUME_TAG{};

}
#endif
