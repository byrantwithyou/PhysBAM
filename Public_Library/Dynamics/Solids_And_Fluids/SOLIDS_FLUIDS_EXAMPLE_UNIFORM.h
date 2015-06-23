//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_EXAMPLE_UNIFORM
//#####################################################################
#ifndef __SOLIDS_FLUIDS_EXAMPLE_UNIFORM__
#define __SOLIDS_FLUIDS_EXAMPLE_UNIFORM__

#include <Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <Solids/Collisions/SOLIDS_COLLISIONS_FORWARD.h>
#include <Fluids/Fluids/FLUID_COLLECTION.h>
#include <Incompressible/Collisions_And_Interactions/INCOMPRESSIBLE_COLLISIONS_FORWARD.h>
#include <Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
#include <Dynamics/Solids_And_Fluids/SPH_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class DEBUG_PARTICLES;

template<class TV>
class SOLIDS_FLUIDS_EXAMPLE_UNIFORM:public SOLIDS_FLUIDS_EXAMPLE<TV>,public LEVELSET_CALLBACKS<TV>,public SPH_CALLBACKS<TV>,
    public FLUIDS_PARAMETERS_CALLBACKS<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<int,TV::m> T_VECTOR_INT;
    typedef ARRAY<char,TV_INT> T_ARRAYS_CHAR;
    typedef typename MATRIX_POLICY<TV>::TRANSFORMATION_MATRIX T_TRANSFORMATION_MATRIX;
    typedef typename TV::SPIN T_ANGULAR_VELOCITY;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<TV> T_FACE_LOOKUP_COLLIDABLE;
    typedef AVERAGING_UNIFORM<TV> T_AVERAGING;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE<TV> BASE;
    using BASE::output_directory;using BASE::first_frame;using BASE::restart;using BASE::Write_Frame_Title;using BASE::solids_parameters;using BASE::stream_type;
    using BASE::solids_fluids_parameters;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::Adjust_Phi_With_Sources;using BASE::minimum_collision_thickness;using FLUIDS_PARAMETERS_CALLBACKS<TV>::Adjust_Density_And_Temperature_With_Sources;
    using FLUIDS_PARAMETERS_CALLBACKS<TV>::Get_Source_Reseed_Mask;using FLUIDS_PARAMETERS_CALLBACKS<TV>::Get_Analytic_Velocities;
    using FLUIDS_PARAMETERS_CALLBACKS<TV>::Get_Source_Velocities;using FLUIDS_PARAMETERS_CALLBACKS<TV>::Get_Object_Velocities; // silence -Woverloaded-virtual

    FLUIDS_PARAMETERS_UNIFORM<TV> fluids_parameters;
    FLUID_COLLECTION<TV> fluid_collection;
    int resolution;
    DEBUG_PARTICLES<TV>& debug_particles;
    bool opt_skip_debug_data;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM(const STREAM_TYPE stream_type,PARSE_ARGS& parse_args,const int number_of_regions,const typename FLUIDS_PARAMETERS<TV>::TYPE type);
    virtual ~SOLIDS_FLUIDS_EXAMPLE_UNIFORM();

    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,ARRAY<T,FACE_INDEX<TV::m> >& V_levelset,const T time) const override
    {if(fluids_parameters.analytic_test) Get_Analytic_Velocities(time);V_levelset=fluid_collection.incompressible_fluid_collection.face_velocities;}

    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) override
    {fluids_parameters.Adjust_Particle_For_Domain_Boundaries(particles,index,V,particle_type,dt,time);}

    void Get_Body_Force(ARRAY<T,FACE_INDEX<TV::m> >& force,const T dt,const T time) override
    {if(fluids_parameters.fire||fluids_parameters.smoke) fluids_parameters.Get_Body_Force(force,dt,time);else PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Update_Fluid_Parameters(const T dt,const T time)
    {fluids_parameters.Update_Fluid_Parameters(dt,time);}

    virtual void Evolve_Density_And_Temperature(const T dt,const T time)
    {fluids_parameters.Evolve_Density_And_Temperature(dt,time);}

    virtual void Apply_Isobaric_Fix(const T dt,const T time)
    {fluids_parameters.Apply_Isobaric_Fix(dt,time);}

//#####################################################################
    void Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision=true,bool add_coupling=true);
    void Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,bool add_collision=true,bool add_coupling=true);
    virtual void Initialize_MPI();
    virtual void Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization(); // called before grids are initialized
    virtual void Initialize_Solid_Fluid_Coupling_After_Grid_Initialization(); // called after grids are initialized
    virtual void Initialize_Compressible_Incompressible_Coupling();
    virtual void Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
    void Set_Dirichlet_Boundary_Conditions(const T time) override;
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity);
    template<class GEOMETRY> void Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity,const ARRAY<bool,FACE_INDEX<TV::m> >& invalid_mask);
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source);
    template<class GEOMETRY> void Adjust_Phi_With_Source(const GEOMETRY& source,const int region,const T_TRANSFORMATION_MATRIX& world_to_source);
    template<class GEOMETRY> void Get_Source_Reseed_Mask(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool,TV_INT>*& cell_centered_mask,const bool reset_mask);
    template<class GEOMETRY> void Adjust_Density_And_Temperature_With_Sources(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const T source_density,
        const T source_temperature);
    void Revalidate_Fluid_Scalars();
    void Revalidate_Phi_After_Modify_Levelset();
    void Revalidate_Fluid_Velocity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
    void Post_Velocity_Advection_Callback(const T dt,const T time){}
    void Get_Object_Velocities(LAPLACE_UNIFORM<TV>* elliptic_solver,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time) override;
    void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<TV>& levelset_multiple,ARRAY<T,FACE_INDEX<TV::m> >& V_levelset,const T time) const override;
    void Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,T maximum_fluid_velocity,const bool include_removed_particle_velocities);
    void Read_Output_Files_Fluids(const int frame) override;
    virtual void Write_Output_Files(const int frame) const override;
    void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time) override;
    void Log_Parameters() const override;
    void After_Construction() override;
//#####################################################################
};
}
#endif
