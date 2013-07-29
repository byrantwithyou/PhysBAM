//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_UNIFORM
//#####################################################################
#ifndef __SOLIDS_FLUIDS_DRIVER_UNIFORM__
#define __SOLIDS_FLUIDS_DRIVER_UNIFORM__    

#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class TV>
class SOLIDS_FLUIDS_DRIVER_UNIFORM:public SOLIDS_FLUIDS_DRIVER<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef typename ARRAY<T,TV_INT>::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;

    typedef SOLIDS_FLUIDS_DRIVER<TV> BASE;
public:
    using BASE::project_at_frame_boundaries;using BASE::current_frame;using BASE::next_dt;using BASE::next_done;using BASE::Write_Time;using BASE::Write_First_Frame;
    using BASE::Write_Last_Frame;using BASE::Write_Substep;using BASE::time;
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>& example;
    T last_dt;
    T restart_dt;
    bool reset_with_restart;

    SOLIDS_FLUIDS_DRIVER_UNIFORM(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>& example_input);
    virtual ~SOLIDS_FLUIDS_DRIVER_UNIFORM();

    bool Simulate_Fluids() const
    {const FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    return (solids_fluids_parameters.mpi_solid_fluid || fluids_parameters.simulate) && (fluids_parameters.smoke || fluids_parameters.fire || fluids_parameters.water || fluids_parameters.sph || fluids_parameters.compressible);}

    bool Simulate_Solids() const
    {SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=example.solid_body_collection.deformable_body_collection;
    return (deformable_body_collection.simulate && deformable_body_collection.particles.Size()) || (solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies && example.solid_body_collection.rigid_body_collection.rigid_body_particles.Size());}

    bool Simulate_Incompressible_Fluids() const
    {const FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters=example.fluids_parameters;
    return fluids_parameters.simulate && (fluids_parameters.smoke || fluids_parameters.fire || fluids_parameters.water || fluids_parameters.sph);}

    bool Two_Way_Coupled() const
    {return example.fluids_parameters.fluid_affects_solid && example.fluids_parameters.solid_affects_fluid;}

//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Initialize_Fluids_Grids();
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    T Compute_Dt(const T time,const T target_time,bool& done);
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Integrate_Fluid_Non_Advection_Forces(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const int substep);
    void Setup_Solids(const T time,const int substep);
    void Setup_Fluids(const T time);
    void Solid_Position_Update(const T dt,const int substep);
    void Rigid_Cluster_Fracture(const T dt_full_advance,const T dt_cfl,const int substep);
    void Project_Fluid(const T dt_projection,const T time_projection,const int substep);
    void Advance_Fluid_One_Time_Step_Implicit_Part_For_Object_Compatibility(const T dt_projection,const T time_projection,const int substep);
    void Calculate_Maximum_Allowable_dt(const T dt,T& min_dt,const int substep,RUNGEKUTTA<T_ARRAYS_DIMENSION_SCALAR>& rungekutta_u);
    void Advect_Fluid(const T dt,const int substep);
    void Solid_Velocity_Update(const T dt,const int substep,const bool done);
    void Advance_Fluid_One_Time_Step_Implicit_Part(const bool done,const T dt,const int substep);
//#####################################################################
};
}
#endif
