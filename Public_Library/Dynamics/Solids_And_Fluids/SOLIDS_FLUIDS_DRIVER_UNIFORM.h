//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_UNIFORM
//#####################################################################
#ifndef __SOLIDS_FLUIDS_DRIVER_UNIFORM__
#define __SOLIDS_FLUIDS_DRIVER_UNIFORM__    

#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class TV>
class SOLIDS_FLUIDS_DRIVER_UNIFORM:public SOLIDS_FLUIDS_DRIVER<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;typedef VECTOR<T,TV::m+2> TV_DIMENSION;
    typedef ARRAY<TV_DIMENSION,TV_INT> T_ARRAYS_DIMENSION_SCALAR;

    typedef SOLIDS_FLUIDS_DRIVER<TV> BASE;
public:
    using BASE::current_frame;using BASE::next_dt;using BASE::next_done;using BASE::Write_Time;
    using BASE::Write_Substep;using BASE::time;
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
    void Initialize() override;
    void Initialize_Fluids_Grids();
    void Advance_To_Target_Time(const T target_time) override;
    void Postprocess_Frame(const int frame) override;
    T Compute_Dt(const T time,const T target_time,bool& done);
    T Compute_Solids_Dt(const T time);
    T Compute_Fluids_Dt(const T time);
    void Write_Output_Files() override;
    void Integrate_Fluid_Non_Advection_Forces(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt);
    void Setup_Solids(const T time);
    void Setup_Fluids(const T time);
    void Solid_Position_Update(const T dt);
    void Rigid_Cluster_Fracture(const T dt_cfl);
    void Project_Fluid(const T dt_projection,const T time_projection);
    void Advance_Fluid_One_Time_Step_Implicit_Part_For_Object_Compatibility(const T dt_projection,const T time_projection);
    void Advect_Fluid(const T dt);
    void Solid_Velocity_Update(const T dt);
    void Advance_Fluid_One_Time_Step_Implicit_Part(const T dt);
    void Advance_One_Time_Step(const T dt);
//#####################################################################
};
}
#endif
