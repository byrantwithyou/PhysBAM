//#####################################################################
// Copyright 2006, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHEETING
//#####################################################################
#ifndef __SHEETING__
#define __SHEETING__

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Particles/SPH_PARTICLES.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class SHEETING:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef VECTOR<int,3> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::viewer_dir;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::solid_body_collection;using BASE::stream_type;using BASE::test_number;using BASE::resolution;using BASE::Adjust_Phi_With_Source;
    using BASE::user_last_frame;
    
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    FLUID_COLLISION_BODY_INACCURATE_UNION<TV> inaccurate_union;
    bool use_variable_density_for_sph,use_two_way_coupling_for_sph,convert_sph_particles_to_fluid,use_analytic_divergence;
    T ballistic_particles_as_percentage_of_target;
    RANDOM_NUMBERS<T> random;
    int bowl;
    VECTOR<T,3> source_center;
    T source_radius;

    ARRAY<CYLINDER<T> > sources;
    ARRAY<MATRIX<T,4> > world_to_source;
    ARRAY<VECTOR<T,3> > source_velocity;

    SHEETING(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,1,fluids_parameters.WATER),rigid_body_collection(solid_body_collection.rigid_body_collection),inaccurate_union(*fluids_parameters.grid),
        use_variable_density_for_sph(false),use_two_way_coupling_for_sph(false),convert_sph_particles_to_fluid(false),ballistic_particles_as_percentage_of_target((T).03),bowl(0)
    {
        parse_args.Parse();
        // set up the standard fluid environment
        if(!this->user_frame_rate) frame_rate=96;
        restart=0;
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_debug_data=true;
        fluids_parameters.write_ghost_values=true;fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.incompressible_iterations=40;
        fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
        fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_vorticity_confinement_fuel=false;
        fluids_parameters.cfl=(T)1.9;
        
        // SPH parameters
        fluids_parameters.use_sph_for_removed_negative_particles=test_number==2;
        use_variable_density_for_sph=true;
        use_two_way_coupling_for_sph=true;
        convert_sph_particles_to_fluid=false;
        use_analytic_divergence=true;
        ballistic_particles_as_percentage_of_target=(T).1;

        // set up the domain
        if(!user_last_frame) last_frame=5000;
        int cells=1*resolution;
        fluids_parameters.grid->Initialize(TV_INT(3*cells+1,4*cells+1,4*cells+1),RANGE<TV>(TV(1,0,0),TV(4,4,4)));

        if(test_number>2){LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Sheeting/Test_Sheeting_%d_Resolution_%d_%d_%d",test_number,
            (fluids_parameters.grid->counts.x-1),(fluids_parameters.grid->counts.y-1),(fluids_parameters.grid->counts.z-1));

        // set up sources for each test case
        source_center=VECTOR<T,3>((T)2.5,(T)3.7,1);source_radius=(T).10;
        sources.Resize(1);source_velocity.Resize(1);world_to_source.Resize(1);
        sources(1)=CYLINDER<T>(VECTOR<T,3>(0,-(T).25,0),VECTOR<T,3>(0,(T).25,0),source_radius);
        ARRAY<MATRIX<T,4> > source_to_world(1);
        MATRIX<T,4> translation=MATRIX<T,4>::Translation_Matrix(source_center);
        source_to_world(1)=translation;
        source_velocity(1)=VECTOR<T,3>(0,-1,0);world_to_source(1)=source_to_world(1).Inverse();
    }

    virtual ~SHEETING()
    {}

//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() override
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=0;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();TV X=iterator.Location();
        fluids_parameters.particle_levelset_evolution->phi(cell)=sources(1).Signed_Distance(world_to_source(1).Homogeneous_Times(X));}
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const VECTOR<T,3>& X) const
{
    return rigid_body_collection.Rigid_Body(bowl).Implicit_Geometry_Extended_Value(X);
}
//#####################################################################
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Bodies() override
{        
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    bowl=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies/anisotropic_bowl",(T).1,true,true,false);
    rigid_body_collection.rigid_body_particles.frame(bowl).t=VECTOR<T,3>((T)2.5,2,1);
    rigid_body_collection.Rigid_Body(bowl).Frame().r=ROTATION<VECTOR<T,3> >::From_Euler_Angles((T)pi/4,0,0);
    rigid_body_collection.Rigid_Body(bowl).Is_Kinematic()=true;
    inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() override
{
    SPH_EVOLUTION_UNIFORM<TV>& sph_evolution=*fluids_parameters.sph_evolution;
    sph_evolution.use_two_way_coupling=use_two_way_coupling_for_sph;
    sph_evolution.use_variable_density_solve=use_variable_density_for_sph;
    sph_evolution.convert_particles_to_fluid=convert_sph_particles_to_fluid;
    sph_evolution.use_analytic_divergence=use_analytic_divergence;
    sph_evolution.target_particles_per_unit_volume=(T).25*fluids_parameters.number_particles_per_cell/fluids_parameters.particle_levelset_evolution->grid.Cell_Size();
    sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override
{
    Get_Source_Velocities(sources(1),world_to_source(1),source_velocity(1));
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    BASE::Set_Kinematic_Positions(frame,time,id);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override
{
    return BASE::Set_Kinematic_Velocities(twist,time,id);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
    Update_Sources(time);
}
//#####################################################################
// Function Update_Sources
//#####################################################################
virtual void Update_Sources(const T time)
{
    T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& removed_positive_particles=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_positive_particles;
    GRID<TV>& grid=*fluids_parameters.grid;
    VECTOR<T,3> source_radius_vector=source_radius*VECTOR<T,3>::All_Ones_Vector();
    RANGE<TV_INT> source_nodes(grid.Clamp_To_Cell(source_center-source_radius_vector),grid.Clamp_To_Cell(source_center+source_radius_vector)+TV_INT::All_Ones_Vector());
    for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(!removed_positive_particles(block)) continue;
        for(int p=0;p<removed_positive_particles(block)->Size();p++)
            if(sources(1).Inside(world_to_source(1).Homogeneous_Times(removed_positive_particles(block)->X(p)),(T)1e-4)) removed_positive_particles(block)->Add_To_Deletion_List(p);
        removed_positive_particles(block)->Delete_Elements_On_Deletion_List();}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) override
{
    Get_Source_Reseed_Mask(sources(1),world_to_source(1),cell_centered_mask,true);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    Adjust_Phi_With_Source(sources(1),world_to_source(1));
    return false;
}
//#####################################################################
};
}
#endif
