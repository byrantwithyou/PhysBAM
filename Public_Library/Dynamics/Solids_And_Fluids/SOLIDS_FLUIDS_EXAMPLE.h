//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_EXAMPLE
//#####################################################################
#ifndef __SOLIDS_FLUIDS_EXAMPLE__
#define __SOLIDS_FLUIDS_EXAMPLE__

#include <Core/Log/DEBUG_SUBSTEPS.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Solids/Solids_Evolution/SOLIDS_EVOLUTION_CALLBACKS.h>
#include <Incompressible/Incompressible_Flows/EXTERNAL_STRAIN_ADJUSTMENT.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class SOLIDS_PARAMETERS;
template<class TV> class SOLIDS_EVOLUTION;
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class SOLIDS_FLUIDS_PARAMETERS;

template<class TV>
class SOLIDS_FLUIDS_EXAMPLE:public EXAMPLE<TV>,public EXTERNAL_STRAIN_ADJUSTMENT<typename TV::SCALAR>,public EXAMPLE_FORCES_AND_VELOCITIES<TV>,
                            public SOLIDS_EVOLUTION_CALLBACKS<TV>,public SOLIDS_FLUIDS_CALLBACKS<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef EXAMPLE<TV> BASE;
public:
    using EXAMPLE_FORCES_AND_VELOCITIES<TV>::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::output_directory;using BASE::frame_title;using BASE::stream_type;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::write_last_frame;using BASE::write_time;using BASE::write_substeps_level;using BASE::Set_Write_Substeps_Level;//using BASE::data_directory;
    using BASE::restart;

protected:
    T minimum_collision_thickness; // needed for ray tracing, etc.
public:
    bool use_melting;
    SOLIDS_PARAMETERS<TV>& solids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_EVOLUTION<TV>* solids_evolution; // defaults to newmark
    bool opt_solidssymmqmr,opt_solidscr,opt_solidscg;

    SOLIDS_FLUIDS_EXAMPLE(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args);
    virtual ~SOLIDS_FLUIDS_EXAMPLE();

    void Set_Minimum_Collision_Thickness(const T minimum_collision_thickness_input=1e-6)
    {minimum_collision_thickness=minimum_collision_thickness_input;}

//#####################################################################
    virtual void Post_Initialization();
    virtual void Preprocess_Frame(const int frame);
    virtual void Postprocess_Frame(const int frame);
    virtual void Preprocess_Substep(const T dt,const T time);
    virtual void Postprocess_Substep(const T dt,const T time);
    virtual void Read_Output_Files_Fluids(const int frame);
    void Log_Parameters() const override;
    // solids
    virtual void Initialize_Bodies();
    virtual void Read_Output_Files_Solids(const int frame);
    // fluids
    virtual void Extrapolate_Phi_Into_Objects(const T time);
    virtual void Postprocess_Phi(const T time);
    virtual bool Adjust_Phi_With_Sources(const T time);
    virtual void Adjust_Phi_With_Objects(const T time);
    virtual void Add_SPH_Particles_For_Sources(const T dt,const T time);
    virtual void Initialize_SPH_Particles();
    virtual void Initialize_Velocities();
    virtual void Initialize_Euler_State();
    virtual void Setup_Initial_Refinement();
    virtual void Initialize_Advection();
    virtual void Clamp_Velocities(const T time);
    // melting
    virtual void Melting_Substep(const T dt,const T time);
    virtual void Modify_Fluid_For_Melting(const T dt,const T time);
    virtual void Update_Melting_Substep_Parameters(const T dt,const T time);
    void After_Initialization() override;
    template<class T_MPI> void Adjust_Output_Directory_For_MPI(const T_MPI mpi);
//#####################################################################
};
}
#endif
