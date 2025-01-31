//#####################################################################
// Copyright 2007, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWO_PHASE
//#####################################################################
#ifndef __TWO_PHASE__
#define __TWO_PHASE__

#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_MULTIPHASE_2D.h>

/***************
1. Rising bubble (same as standard test 14)
2. Rising bubble with a rigid sphere at top of it
3. Rising bubble with a rigid sphere at top of it and a velocity inlet from the right.
***************/

namespace PhysBAM{

template<class T_input>
class TWO_PHASE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::viewer_dir;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;using BASE::fluid_collection;
    using BASE::solid_body_collection;using BASE::stream_type;using BASE::test_number;using BASE::resolution;
    using BASE::user_last_frame;
    
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    FLUID_COLLISION_BODY_INACCURATE_UNION<TV> inaccurate_union;
    bool use_inaccurate_body_collisions;
    int sphere;

    TWO_PHASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,2,fluids_parameters.WATER),
        rigid_body_collection(solid_body_collection.rigid_body_collection),inaccurate_union(*fluids_parameters.grid),use_inaccurate_body_collisions(true),sphere(0)
    {
        parse_args.Parse();
        int cells=1*resolution;
        // set up the standard fluid environment
        if(!this->user_frame_rate) frame_rate=24;
        restart=0;
        if(!user_last_frame) last_frame=300;
        if(!this->user_frame_rate) frame_rate=400;
        fluids_parameters.domain_walls[0][0]=true;
        fluids_parameters.domain_walls[0][1]=true;
        fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.viscosity=(T)0;
        fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;
        fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_debug_data=true;
        fluids_parameters.write_ghost_values=true;
        fluids_parameters.write_removed_positive_particles=true;fluids_parameters.write_removed_negative_particles=true;

        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.incompressible_iterations=100;
        fluids_parameters.use_removed_positive_particles=false;
        fluids_parameters.use_removed_negative_particles=false;
        fluids_parameters.write_removed_positive_particles=false;
        fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.second_order_cut_cell_method=true;
        fluids_parameters.store_particle_ids=true;
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_vorticity_confinement_fuel=false;
 
        // set up the example parameters
        if(test_number!=1 && test_number!=2 && test_number!=3){LOG::cout<<"Unrecognized example: "<<test_number<<std::endl;PHYSBAM_FATAL_ERROR();}
        fluids_parameters.densities(0)=(T)1.226;fluids_parameters.densities(1)=1000;
        fluids_parameters.surface_tensions(0,1)=fluids_parameters.surface_tensions(1,0)=(T).0728;
        fluids_parameters.viscosities(0)=(T).0000178;
        fluids_parameters.viscosities(1)=(T).001137;
        fluids_parameters.implicit_viscosity=false;
        fluids_parameters.incompressible_iterations=200;
        fluids_parameters.implicit_viscosity_iterations=200;
        fluids_parameters.solve_neumann_regions=true;

        // set up the domain
        GRID<TV>& grid=*fluids_parameters.grid;
        grid.Initialize(TV_INT(10*cells+1,20*cells+1),RANGE<TV>(TV((T)-.01,(T)-.01),TV((T).01,(T).02)));

        if(!this->user_output_directory)
            viewer_dir.output_directory=LOG::sprintf("Two_Phase/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));
        LOG::cout<<"output directory="<<viewer_dir.output_directory<<std::endl;

        // set example-specific parameters
        fluids_parameters.object_friction=(test_number==2||test_number==3)?(T)1:0;
        inaccurate_union.Initialize_Grid_Structures();
    }

    ~TWO_PHASE()
    {}

//#####################################################################
// Function Initial_Phi
//#####################################################################
void Initialize_Advection()    override
{
    if(test_number==2||test_number==3) fluids_parameters.Use_Fluid_Coupling_Defaults();
    else fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    T radius=(T)1/(T)300;
    for(CELL_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();
        fluids_parameters.particle_levelset_evolution_multiple->phis(0)(iterator.Cell_Index())=X.Magnitude()-radius;// center is at 0,0,0
        fluids_parameters.particle_levelset_evolution_multiple->phis(1)(iterator.Cell_Index())=radius-X.Magnitude();}
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const TV& X) const
{
    if(test_number==2||test_number==3) return rigid_body_collection.Rigid_Body(sphere).Implicit_Geometry_Extended_Value(X);
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
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    if(test_number==2||test_number==3){
        sphere=rigid_body_collection.Add_Rigid_Body(data_directory+"/Rigid_Bodies_2D/circle",(T).001,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T).0,(T).005);
        rigid_body_collection.Rigid_Body(sphere).Is_Kinematic()=true;}
    if(use_inaccurate_body_collisions){
        inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);}
    else fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Set_Dirichlet_Boundary_Conditions(time);
    if(test_number!=3) return;
    //TODO: Is there a good place to set Neuman conditions?
    ARRAY<bool,FACE_INDEX<TV::m> >& psi_N=fluids_parameters.incompressible_multiphase->projection.elliptic_solver->psi_N;
    ARRAY<T,FACE_INDEX<TV::m> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;

    RANGE<TV_INT> right_grid_cells=RANGE<TV_INT>(TV_INT(fluids_parameters.grid->counts.x-2,1),fluids_parameters.grid->Numbers_Of_Cells());
    for(int axis=0;axis<TV::m;axis++){
        RANGE<TV_INT> right_grid_faces=right_grid_cells+RANGE<TV_INT>(TV_INT(),TV_INT::Axis_Vector(axis));
        for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid,right_grid_faces,axis);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();
            psi_N.Component(axis)(face)=true;face_velocities.Component(axis)(face)=axis==1?(T)-1:(T)0;}}
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Construct_Levelsets_For_Objects(time);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
};
}
#endif
