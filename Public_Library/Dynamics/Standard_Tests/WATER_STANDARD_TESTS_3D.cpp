//#####################################################################
// Copyright 2005-2008, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_3D
//#####################################################################
#include <Core/Matrices/MATRIX_4X4.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Basic_Geometry/TORUS.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
using namespace PhysBAM;
//#####################################################################
// Function WATER_STANDARD_TESTS_3D
//#####################################################################
template<class TV> WATER_STANDARD_TESTS_3D<TV>::
WATER_STANDARD_TESTS_3D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS<TV>& fluids_parameters,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :example(example),fluids_parameters(fluids_parameters),rigid_body_collection(rigid_body_collection_input),inaccurate_union(*fluids_parameters.grid),
    use_inaccurate_body_collisions(false),use_variable_density_for_sph(false),use_two_way_coupling_for_sph(false),convert_sph_particles_to_fluid(false),use_analytic_divergence(false),
    use_analytic_divergence_for_expansion_only(false),adjust_cell_weights_on_neumann_boundaries(false),enforce_density_near_interface(true),flip_ratio(1),target_particles_per_unit_volume(1),
    neumann_boundary_slip_multiplier(1),ballistic_particles_as_percentage_of_target((T).03),particle_targeting_time(1),test_number(0),sphere(0)
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Initialize(const int test_number_input,const int resolution)
{
    test_number=test_number_input;
    LOG::cout<<"Running Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;

    // set up the standard fluid environment
    if(!example.user_frame_rate) example.frame_rate=24;
    example.restart=false;example.restart_frame=0;
    fluids_parameters.domain_walls(0)(0)=true;fluids_parameters.domain_walls(0)(1)=true;fluids_parameters.domain_walls(1)(0)=true;
    fluids_parameters.domain_walls(1)(1)=false;fluids_parameters.domain_walls(2)(0)=true;fluids_parameters.domain_walls(2)(1)=true;
    fluids_parameters.number_particles_per_cell=16;
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

    auto lf=[this](int f)
        {
            if(!example.user_last_frame)
                example.last_frame=500;
        };
    // set up the domain
    int cells=1*resolution;
    if(test_number==1){
        lf(10);
        grid.Initialize(TV_INT(15*cells+1,10*cells+1,20*cells+1),RANGE<TV>(TV(),TV((T)1.5,1,2)));}
    else if(test_number==2){
        lf(100);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==3){
        lf(150);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==4){
        lf(2000);
        grid.Initialize(TV_INT(15*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(),TV((T)1.5,1,1)));
        motion_curve.Add_Control_Point(0,TV((T)1.25,(T).55,(T).5));
        motion_curve.Add_Control_Point((T).075,TV((T).8,(T).1,(T).5)); // .03 was old
        motion_curve.Add_Control_Point(3,TV((T).8,(T).1,(T).5));}
    else if(test_number==5){
        lf(150);
        grid.Initialize(TV_INT(10*cells+1,25*cells+1,10*cells+1),RANGE<TV>(TV(),TV((T).1,(T).25,(T).1)));}
    else if(test_number==6){
        lf(50);
        grid.Initialize(TV_INT(25*cells+1,10*cells+1,25*cells+1),RANGE<TV>(TV(),TV(1,(T).4,1)));}
    else if(test_number==7){
        fluids_parameters.gravity=TV();
        lf(100);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==8){
        lf(200);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==9){
        lf(2000);
        grid.Initialize(TV_INT(16*cells+1,6*cells+1,8*cells+1),RANGE<TV>(TV(),TV(4,(T)1.5,2)));}
    else if(test_number==10){
        lf(2000);
        grid.Initialize(TV_INT(6*cells+1,12*cells+1,2*cells+1),RANGE<TV>(TV(),TV((T)1.5,3,(T).5)));}
    else if(test_number==11){
        fluids_parameters.gravity=TV();
        lf(100);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==12){
        fluids_parameters.viscosity=(T)1000;fluids_parameters.implicit_viscosity=true;
        lf(100);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==13){
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=true;
        fluids_parameters.use_explicit_part_of_implicit_viscosity=true;
        lf(100);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==15){
        fluids_parameters.analytic_test=true; //TODO: analytic test
        fluids_parameters.viscosity=(T)1000;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=false;
        fluids_parameters.use_explicit_part_of_implicit_viscosity=false;
        lf(100);
        grid.Initialize(TV_INT(2*cells+1,3*cells+1,2*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==16){
        fluids_parameters.analytic_test=true; //TODO: analytic test
        lf(100);
        grid.Initialize(TV_INT(2*cells+1,3*cells+1,2*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else if(test_number==19){
        lf(1000);
        grid.Initialize(TV_INT(201,81,201),RANGE<TV>(TV(),TV((T)1.5,1,2)));
        fluids_parameters.cfl=(T)0.9;
        fluids_parameters.incompressible_iterations=100;
        //grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));
    }
    else if(test_number==20){ //MIKE example
        lf(2000);
        //grid.Initialize(TV_INT(15*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(),TV((T)1.5,1,1)));
        grid.Initialize(TV_INT(201,81,201),RANGE<TV>(TV(),TV((T)1.5,1,2)));
        fluids_parameters.cfl=(T)0.9;
        fluids_parameters.incompressible_iterations=100;
        //fluids_parameters.cg_restart_iterations=20;
        //fluids_parameters.evolution_solver_type=krylov_solver_cr;
        motion_curve.Add_Control_Point(0,TV((T)1,(T).1,(T)1));
        motion_curve.Add_Control_Point((T).075,TV((T)1,(T).1,(T)1));
        motion_curve.Add_Control_Point(3,TV((T)1,(T).55,(T)1));}
    else if(test_number==21){
        lf(1000);
        fluids_parameters.cfl=(T)1.9;
        fluids_parameters.incompressible_iterations=100;
        grid.Initialize(TV_INT(201,81,201),RANGE<TV>(TV(),TV(2,(T).8,2)));}
    else if(test_number==22){
        lf(1000);
        fluids_parameters.cfl=(T)1.9;
        fluids_parameters.incompressible_iterations=20;
        grid.Initialize(TV_INT(201,41,161),RANGE<TV>(TV(),TV((T)2,(T).4,(T)1.6)));}
    else if(test_number==23){
        lf(100);
        grid.Initialize(TV_INT(10*cells+1,15*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5,1)));}
    else{
        LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

    if(!example.user_output_directory)
        example.viewer_dir.output_directory=LOG::sprintf("Standard_Tests/Test_%d__Resolution_%d_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1),(grid.counts.z-1));

    // set up sources for each test case
    if(test_number==3){
        TV domain_center=grid.domain.Center();
        sources.Resize(4);source_velocity.Resize(4);world_to_source.Resize(4);
        sources(0)=sources(1)=sources(2)=sources(3)=CYLINDER<T>(TV((T)-.0103432,0,0),TV((T).15324,0,0),(T).10787);
        MATRIX<T,4> rot=MATRIX<T,4>::Translation_Matrix(domain_center)*MATRIX<T,4>::Rotation_Matrix_Y_Axis((T)pi/2)*MATRIX<T,4>::Translation_Matrix(-domain_center);
        ARRAY<MATRIX<T,4> > source_to_world(4);
        source_to_world(0)=MATRIX<T,4>::Translation_Matrix(TV((T).38,1,(T).5))*MATRIX<T,4>::Rotation_Matrix_Z_Axis((T)pi);
        source_to_world(1)=rot*rot*source_to_world(0);
        source_to_world(2)=rot*source_to_world(0);
        source_to_world(3)=rot*rot*rot*source_to_world(0);
        for(int i=0;i<4;i++){source_velocity(i)=source_to_world(i).Extract_Rotation()*TV((T)1.2,0,0);world_to_source(i)=source_to_world(i).Inverse();}}
    if(test_number==5){
        world_to_source.Append(MATRIX<T,4>::Identity_Matrix());
        sources.Append(CYLINDER<T>(TV((T).03,(T).18,(T).05),TV((T).06,(T).2,(T).05),(T).015));
        source_velocity.Append((T).8*TV((T)-.03,(T)-.02,0).Normalized());}
    if(test_number==8){
        TV domain_center=grid.domain.Center();
        sources.Resize(1);source_velocity.Resize(1);world_to_source.Resize(1);
        sources(0)=CYLINDER<T>(TV(),TV((T).2,0,0),(T).05);
        ARRAY<MATRIX<T,4> > source_to_world(1);
        MATRIX<T,4> translation=MATRIX<T,4>::Translation_Matrix(domain_center);
        source_to_world(0)=translation*MATRIX<T,4>::Rotation_Matrix_Z_Axis((T).5*(T)pi);world_to_source(0)=source_to_world(0).Inverse();
        source_velocity(0)=TV(0,(T)-1,0);}
    
    // set example-specific parameters
    fluids_parameters.object_friction=(test_number==4 || test_number==20)?(T)1:0;
    if(test_number==8){
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        use_two_way_coupling_for_sph=false;
        convert_sph_particles_to_fluid=false;
        use_variable_density_for_sph=true;}
    if(test_number==9){
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        convert_sph_particles_to_fluid=false;
        use_analytic_divergence=true;
        use_two_way_coupling_for_sph=true;
        use_variable_density_for_sph=true;}
    if(test_number==10){
        fluids_parameters.cfl=(T)1.9;
        flip_ratio=1;
        target_particles_per_unit_volume=(T).5*fluids_parameters.number_particles_per_cell/grid.Cell_Size();
        neumann_boundary_slip_multiplier=0;
        adjust_cell_weights_on_neumann_boundaries=true;
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        convert_sph_particles_to_fluid=false;
        use_analytic_divergence=true;
        use_analytic_divergence_for_expansion_only=false;
        particle_targeting_time=(T).25;
        use_two_way_coupling_for_sph=false;
        use_variable_density_for_sph=false;
        enforce_density_near_interface=true;}

    if(test_number==4||test_number==5||test_number==6||test_number==20){fluids_parameters.solid_affects_fluid=true;fluids_parameters.fluid_affects_solid=false;}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Initialize_Advection(const bool always_use_objects)
{
    PHYSBAM_ASSERT(!always_use_objects,"fluids_parameters.solid_affects_fluids overrides this functionality");  // Deprecated parameter
    if(always_use_objects||test_number==4||test_number==5||test_number==6||test_number==20) fluids_parameters.Use_Fluid_Coupling_Defaults();
    else fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class TV> TV WATER_STANDARD_TESTS_3D<TV>::
Initial_Velocity(const TV& X) const
{
    if(test_number==7 || test_number==11) return TV(0,1,0);
    if(test_number==12 || test_number==13) return TV(-1,0,0);
    if(test_number==22)
    {
        TV u;
        height_noise_random.Fill_Uniform(u,-(T).16,(T).16);
        return u;
    }
    return TV();
}
//#####################################################################
// Class Get_Variable_Viscosity
//#####################################################################
template<class T>
void Get_Variable_Viscosity_Helper(GRID<VECTOR<T,3> >& grid,const int test_number,ARRAY<T,VECTOR<int,3> >& variable_viscosity,const T time)
{
    if(test_number==13){
        for(CELL_ITERATOR<VECTOR<T,3> > iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.2)?(T)0:(T)1000;}
    else if(test_number==15){
        for(CELL_ITERATOR<VECTOR<T,3> > iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.5)?(T)0:(T)1000;}
}
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Get_Variable_Viscosity(ARRAY<T,VECTOR<int,3> >& variable_viscosity,const T time) const
{
    Get_Variable_Viscosity_Helper(*fluids_parameters.grid,test_number,variable_viscosity,time);
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
template<class TV> typename TV::SCALAR WATER_STANDARD_TESTS_3D<TV>::
Initial_Phi(const TV& X) const
{
    T phi=1;
    if(test_number==1) phi=X.y-(T).32424;
    else if(test_number==2){
        static SPHERE<TV> sphere((TV((T).5,(T).75,(T).5)),(T).2);
        phi=min(sphere.Signed_Distance(X),X.y-(T).412134);}
    else if(test_number==4) phi=X.y-(T).400235234;
    else if(test_number==6){
        static SPHERE<TV> sphere((TV((T).7,(T)((T).2*.7),(T).5)),(T).2);
        phi=sphere.Signed_Distance(X);}
    else if(test_number==7){
        static SPHERE<TV> sphere((TV((T).5,(T).5,(T).5)),(T).2);
        phi=sphere.Signed_Distance(X);}
    else if(test_number==9){
        static RANGE<TV> box((TV((T)3.4,1,(T).9)),(TV((T)3.6,(T)1.3,(T)1.1)));
        phi=min(box.Signed_Distance(X),X.y-(T).4);} 
    else if(test_number==11){
        static TV normal((T)1,(T)1,0),location((T).5,(T).75,(T).5);
        static PLANE<T> plane(normal,location);
        phi=plane.Signed_Distance(X);return phi;}
    else if(test_number==12 || test_number==13){
        static SPHERE<TV> sphere((TV((T).5,(T).75,(T).5)),(T).2);
        phi=sphere.Signed_Distance(X);}
    else if(test_number==15){
        static SPHERE<TV> sphere((TV((T).5,(T)1,(T).5)),(T).2);
        phi=sphere.Signed_Distance(X);}
    else if(test_number==16){
        static SPHERE<TV> sphere((TV((T).5,(T).75,(T).5)),(T).2);
        phi=sphere.Signed_Distance(X);}
    else if(test_number==19){
        static SPHERE<TV> sphere((TV((T)1.25,(T).55,(T).5)),(T).1);
        phi=min(sphere.Signed_Distance(X),X.y-(T).412134);}
    else if(test_number==20){
        static SPHERE<TV> sphere((TV((T)1,(T).1,(T)1)),(T).1);
        phi=-sphere.Signed_Distance(X)+(T)1e-5;
        if(phi<0) phi=X.y-(T).400235234;}
    else if(test_number==21) phi=X.y-(T).32424;
    else if(test_number==22){
        T height=height_noise_random.Get_Uniform_Number(grid.X(TV_INT(0,20,0)).y,grid.X(TV_INT(0,21,0)).y);
        phi=X.y-height;}
    else if(test_number==23){
        static SPHERE<TV> sphere((TV((T).7,(T).7,(T).7)),(T).1);
        static TORUS<T> torus((TV((T).3,(T).5,(T).7)),TV(1,0,0),(T).05,(T).2);
        static CYLINDER<T> cylinder((TV((T).5,(T).8,(T).3)),(TV((T).5,(T)1,(T).3)),(T).1);
        phi=min(sphere.Signed_Distance(X),X.y-(T).412134);
        phi=min(torus.Signed_Distance(X),phi);
        phi=min(cylinder.Signed_Distance(X),phi);
    }
    for(int s=0;s<sources.m;s++)phi=min(phi,sources(s).Signed_Distance(world_to_source(s).Homogeneous_Times(X)));
    return phi;
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
template<class TV> typename TV::SCALAR WATER_STANDARD_TESTS_3D<TV>::
Initial_Phi_Object(const TV& X) const
{
    if(test_number==4 || test_number==20) return rigid_body_collection.Rigid_Body(sphere).Implicit_Geometry_Extended_Value(X);
    else if(test_number==5){
        static CYLINDER<T> glass((TV((T).05,(T)-.1,(T).05)),TV((T).05,1,(T).05),(T).045);
        return -glass.Signed_Distance(X);}
    else if(test_number==6){
        static PLANE<T> plane((TV((T)-.2,1,0)).Normalized(),TV());
        return plane.Signed_Distance(X);}
    else return 1;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Initialize_Bodies()
{
    if(test_number==4){
        sphere=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/sphere",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T)1.25,(T).55,(T).5);
        rigid_body_collection.rigid_body_particles.kinematic(sphere)=true;}
    else if(test_number==20){
        sphere=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies/sphere",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T).8,(T).1,(T).5);
        rigid_body_collection.rigid_body_particles.kinematic(sphere)=true;}
    if(use_inaccurate_body_collisions){
        inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);}
    else fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
}
//#####################################################################
// Function Update_Sources
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Update_Sources(const T time)
{
    if(test_number==3)
        if(time>4) sources.Clean_Memory();
    if(test_number==5)
        if(time>4) sources.Clean_Memory();
    if(test_number==8)
        if(time>2) sources.Clean_Memory();
    if(test_number==10) {//TODO : move this to a better callback function
        FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<TV>&>(fluids_parameters);
        SPH_EVOLUTION_UNIFORM<TV>& sph_evolution=*fluids_parameters_uniform.sph_evolution;
        if(time>7){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume/2;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}
        else if(time>5){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}
        else if(time>3){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume*2;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}}
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if((test_number==4 || test_number==20) && id==sphere) frame.t=motion_curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> bool WATER_STANDARD_TESTS_3D<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if((test_number==4 || test_number==20) && id==sphere){twist.linear=motion_curve.Derivative(time);return true;}
    return false;
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Limit_Dt(T& dt,const T time)
{
    if(test_number==4 || test_number==20){
        TV velocity=rigid_body_collection.rigid_body_particles.twist(sphere).linear;
        T rigid_dt_denominator=abs(velocity.x)/grid.dX.x+abs(velocity.y)/grid.dX.y+abs(velocity.z)/grid.dX.z;
        if(rigid_dt_denominator>1e-8) dt=min(dt,1/rigid_dt_denominator);}
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
template<class T> void
Initialize_SPH_Particles_Helper(int test_number,WATER_STANDARD_TESTS_3D<VECTOR<T,3> >& tests,FLUIDS_PARAMETERS<VECTOR<T,3> >& fluids_parameters)
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;

    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<TV>&>(fluids_parameters);
    PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=fluids_parameters_uniform.particle_levelset_evolution->Particle_Levelset(0);
    SPH_EVOLUTION_UNIFORM<TV>& sph_evolution=*fluids_parameters_uniform.sph_evolution;
    GRID<TV>& grid=fluids_parameters_uniform.particle_levelset_evolution->grid;

    sph_evolution.use_two_way_coupling=tests.use_two_way_coupling_for_sph;
    sph_evolution.use_variable_density_solve=tests.use_variable_density_for_sph;
    sph_evolution.convert_particles_to_fluid=tests.convert_sph_particles_to_fluid;
    sph_evolution.use_analytic_divergence=tests.use_analytic_divergence;
    sph_evolution.use_analytic_divergence_for_expansion_only=tests.use_analytic_divergence_for_expansion_only;
    sph_evolution.flip_ratio=tests.flip_ratio;
    sph_evolution.neumann_boundary_slip_multiplier=tests.neumann_boundary_slip_multiplier;
    sph_evolution.adjust_cell_weights_on_neumann_boundaries=tests.adjust_cell_weights_on_neumann_boundaries;
    sph_evolution.enforce_density_near_interface=tests.enforce_density_near_interface;
    sph_evolution.particle_targeting_time=tests.particle_targeting_time;
    sph_evolution.target_particles_per_unit_volume=fluids_parameters.number_particles_per_cell/grid.Cell_Size();
    sph_evolution.ballistic_particles_per_unit_volume=tests.ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;

    int particle_id=0,number_of_sph_particles=0;

    // TODO:  Make examples 8-9 work with MPI
    if(test_number==8){
        RANGE<TV> particle_region((TV((T).2,1,(T).2)),TV((T).3,(T)1.2,(T).3));
        number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume);
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        for(int i=0;i<number_of_sph_particles;i++){
            TV X;
            random.Fill_Uniform(X,TV((T).2,1,(T).2),TV((T).3,(T)1.2,(T).3));
            TV_INT block=grid.Block_Index(X,3);
            if(!removed_negative_particles(block)) removed_negative_particles(block)=particle_levelset.template_removed_particles.Clone();
            int id=removed_negative_particles(block)->Add_Element();
            (*removed_negative_particles(block)->template Get_Array<int>("id"))(id)=particle_id++;
            removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}}
    else if(test_number==9){
        sph_evolution.target_particles_per_unit_volume=200000;
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        RANGE<TV> particle_region((TV((T).4,1,(T).9)),TV((T).6,(T)1.3,(T)1.1));
        T particle_multiplier=3;
        number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume/particle_multiplier);
        for(int region=0;region<3;region++){
            for(int i=0;i<number_of_sph_particles;i++){
                TV X;
                random.Fill_Uniform(X,particle_region);
                TV_INT block=grid.Block_Index(X,3);
                if(!removed_negative_particles(block)) removed_negative_particles(block)=particle_levelset.template_removed_particles.Clone();
                int id=removed_negative_particles(block)->Add_Element();
                (*removed_negative_particles(block)->template Get_Array<int>("id"))(id)=particle_id++;
                removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}
            particle_region.min_corner.x++;
            particle_region.max_corner.x++;
            number_of_sph_particles=int(particle_multiplier*number_of_sph_particles);}}
    else if(test_number==10){
        sph_evolution.target_particles_per_unit_volume=tests.target_particles_per_unit_volume;
        sph_evolution.ballistic_particles_per_unit_volume=tests.ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        GRID<TV>  *grid_global;
        if(fluids_parameters_uniform.mpi_grid) grid_global=&fluids_parameters_uniform.mpi_grid->global_grid;
        else grid_global=&grid;
        //RANGE<TV> particle_region(TV(fluids_parameters_uniform.mpi_grid->global_grid.xmin,fluids_parameters_uniform.mpi_grid->global_grid.ymin,fluids_parameters_uniform.mpi_grid->global_grid.zmin),TV(fluids_parameters_uniform.mpi_grid->global_grid.xmax,fluids_parameters_uniform.mpi_grid->global_grid.ymax/3,fluids_parameters_uniform.mpi_grid->global_grid.zmax));
        RANGE<TV> particle_region(RANGE<TV>::Intersect(grid_global->domain,grid.domain));
        T particle_multiplier=(T).5;//to start compressed
        number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume/particle_multiplier);
        for(int i=0;i<number_of_sph_particles;i++){
            TV X;
            random.Fill_Uniform(X,particle_region);
            TV_INT block=grid.Block_Index(X,3);
            if(!removed_negative_particles(block)) removed_negative_particles(block)=particle_levelset.template_removed_particles.Clone();
            int id=removed_negative_particles(block)->Add_Element();
            (*removed_negative_particles(block)->template Get_Array<int>("id"))(id)=particle_id++;
            removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}}
}
template<class TV> void WATER_STANDARD_TESTS_3D<TV>::
Initialize_SPH_Particles()
{
    Initialize_SPH_Particles_Helper(test_number,*this,fluids_parameters);
}
//#####################################################################
// Function Analytic_Velocity
//#####################################################################
template<class TV> TV  WATER_STANDARD_TESTS_3D<TV>::
Analytic_Velocity(const T time,const TV& location) const
{
    if(test_number==15) return TV(0,(T)-1,0);
    else if(test_number==16) return (location-TV((T).5,(T).75,(T).5)).Normalized();
    return TV();
}
//#####################################################################
namespace PhysBAM{
template class WATER_STANDARD_TESTS_3D<VECTOR<float,3> >;
template class WATER_STANDARD_TESTS_3D<VECTOR<double,3> >;
}
