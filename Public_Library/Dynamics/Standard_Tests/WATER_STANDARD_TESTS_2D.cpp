//#####################################################################
// Copyright 2005-2008, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_STANDARD_TESTS_2D
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Matrices/FRAME.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> WATER_STANDARD_TESTS_2D<TV>::
WATER_STANDARD_TESTS_2D(SOLIDS_FLUIDS_EXAMPLE<TV>& example,FLUIDS_PARAMETERS<TV>& fluids_parameters,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :example(example),fluids_parameters(fluids_parameters),rigid_body_collection(rigid_body_collection_input),inaccurate_union(*fluids_parameters.grid),
    use_inaccurate_body_collisions(false),use_variable_density_for_sph(false),use_two_way_coupling_for_sph(false),convert_sph_particles_to_fluid(false),use_analytic_divergence(false),
    use_analytic_divergence_for_expansion_only(false),adjust_cell_weights_on_neumann_boundaries(false),enforce_density_near_interface(true),flip_ratio(1),target_particles_per_unit_volume(1),
    neumann_boundary_slip_multiplier(0),ballistic_particles_as_percentage_of_target((T).03),particle_targeting_time(1),test_number(0),sphere(0)
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Initialize(const int test_number_input,const int resolution)
{
    test_number=test_number_input;
    LOG::cout<<"Running Standard Test Number "<<test_number<<" at resolution "<<resolution<<std::endl;

    // set up the standard fluid environment
    if(!example.user_frame_rate) example.frame_rate=24;
    example.restart=false;example.restart_frame=0;
    fluids_parameters.domain_walls(0)(0)=true;fluids_parameters.domain_walls(0)(1)=true;fluids_parameters.domain_walls(1)(0)=true;fluids_parameters.domain_walls(1)(1)=false;
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

    // set up the domain
    int cells=1*resolution;
    if(test_number==1){
        if(!example.user_last_frame) example.last_frame=10;
        grid.Initialize(TV_INT(15*cells+1,10*cells+1),RANGE<TV>(TV(),TV((T)1.5,1)));}
    else if(test_number==2){
        if(!example.user_last_frame) example.last_frame=100;
        grid.Initialize(TV_INT(10*cells+1,15*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5)));}
    else if(test_number==3){
        if(!example.user_last_frame) example.last_frame=150;
        grid.Initialize(TV_INT(10*cells+1,15*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5)));}
    else if(test_number==4){
        if(!example.user_last_frame) example.last_frame=200;
        grid.Initialize(TV_INT(15*cells+1,10*cells+1),RANGE<TV>(TV(),TV((T)1.5,1)));
        motion_curve.Add_Control_Point(0,TV((T)1.25,(T).55));
        motion_curve.Add_Control_Point((T).2,TV((T).8,(T).1)); // .03 was old
        motion_curve.Add_Control_Point(3,TV((T).8,(T).1));}
    else if(test_number==5){
        if(!example.user_last_frame) example.last_frame=150;
        grid.Initialize(TV_INT(10*cells+1,25*cells+1),RANGE<TV>(TV(),TV((T).1,(T).25)));}
    else if(test_number==6){
        if(!example.user_last_frame) example.last_frame=50;
        grid.Initialize(TV_INT(25*cells+1,10*cells+1),RANGE<TV>(TV(),TV(1,(T).4)));}
    else if(test_number==7){
        fluids_parameters.gravity=TV();
        if(!example.user_last_frame) example.last_frame=100;
        grid.Initialize(TV_INT(10*cells+1,15*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5)));}
    else if(test_number==8){
        if(!example.user_last_frame) example.last_frame=1000;
        grid.Initialize(TV_INT(10*cells+1,15*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5)));}
    else if(test_number==9){
        if(!example.user_last_frame) example.last_frame=1000;
        grid.Initialize(TV_INT(16*cells+1,6*cells+1),RANGE<TV>(TV(),TV(4,(T)1.5)));}
    else if(test_number==10){
        if(!example.user_last_frame) example.last_frame=1000;
        grid.Initialize(TV_INT(6*cells+1,12*cells+1),RANGE<TV>(TV(),TV((T)1.5,3)));}
    else if(test_number==11){
        fluids_parameters.gravity=TV();
        if(!example.user_last_frame) example.last_frame=100;
        grid.Initialize(TV_INT(4*cells+1,6*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5)));}
    else if(test_number==12){
        fluids_parameters.viscosity=(T)1000;fluids_parameters.implicit_viscosity=true;
        if(!example.user_last_frame) example.last_frame=100;
        grid.Initialize(TV_INT(10*cells+1,15*cells+1),RANGE<TV>(TV(),TV(1,(T)1.5)));}
    else if(test_number==13){
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=true;
        fluids_parameters.use_explicit_part_of_implicit_viscosity=true;
        if(!example.user_last_frame) example.last_frame=100;
        grid.Initialize(TV_INT(10*cells,15*cells),RANGE<TV>(TV(0,0),TV(1,(T)1.5)),true);}
    else if(test_number==14){
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=true;
        fluids_parameters.use_explicit_part_of_implicit_viscosity=true;
        if(!example.user_last_frame) example.last_frame=100;
        grid.Initialize(TV_INT(10*cells,15*cells),RANGE<TV>(TV(0,0),TV(1,(T)1.5)),true);}
    else if(test_number==15){
        fluids_parameters.analytic_test=true;
        //fluids_parameters.viscosity=(T)1000;fluids_parameters.implicit_viscosity=true;fluids_parameters.variable_viscosity=false;
        //fluids_parameters.use_explicit_part_of_implicit_viscosity=false;
        if(!example.user_last_frame) example.last_frame=100;
        grid.Initialize(TV_INT(2*cells,3*cells),RANGE<TV>(TV(0,0),TV(1,(T)1.5)),true);}
    else if(test_number==20){ //MIKE example
        if(!example.user_last_frame) example.last_frame=2000;
        //grid.Initialize(TV_INT(15*cells+1,10*cells+1,10*cells+1),RANGE<TV>(TV(0,0,0),TV((T)1.5,1,1)));
        grid.Initialize(TV_INT(201,81),RANGE<TV>(TV(),TV((T)2,(T).8)));
        fluids_parameters.cfl=(T)1.9;
        fluids_parameters.incompressible_iterations=100;
        //fluids_parameters.cg_restart_iterations=20;
        //fluids_parameters.evolution_solver_type=krylov_solver_cr;
        motion_curve.Add_Control_Point(0,TV((T)1,(T).1));
        motion_curve.Add_Control_Point((T).075,TV((T)1,(T).1));
        motion_curve.Add_Control_Point(3,TV((T)1,(T).55));}
    else if(test_number==21){
        if(!example.user_last_frame) example.last_frame=1000;
        fluids_parameters.cfl=(T)0.9;
        fluids_parameters.incompressible_iterations=100;
        grid.Initialize(TV_INT(201,81),RANGE<TV>(TV(),TV((T)1,(T).8)));}
    else{
        LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

    if(!example.user_output_directory)
        example.output_directory=LOG::sprintf("Standard_Tests/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));

    // set up sources for each test case
    if(test_number==3){
        TV domain_center=grid.domain.Center();domain_center.y=(T)1;
        sources.Resize(2);source_velocity.Resize(2);world_to_source.Resize(2);
        sources(0)=sources(1)=RANGE<TV>(TV((T).10787,-(T).10787),TV((T).2714532,(T).10787));
        ARRAY<MATRIX<T,3> > source_to_world(2);
        MATRIX<T,3> rotation=MATRIX<T,3>::Rotation_Matrix_Z_Axis((T)pi);
        MATRIX<T,3> translation=MATRIX<T,3>::Translation_Matrix(domain_center);
        source_to_world(0)=translation;
        source_to_world(1)=translation*rotation;
        for(int i=0;i<2;i++){source_velocity(i)=source_to_world(i).Extract_Rotation()*TV((T).5,0);world_to_source(i)=source_to_world(i).Inverse();}}
    if(test_number==5){
        world_to_source.Append(MATRIX<T,3>::Rotation_Matrix_Z_Axis((T)pi/(T)4)*MATRIX<T,3>::Translation_Matrix(TV((T)-.03,(T)-.23)));
        sources.Append(RANGE<TV>(TV(-(T).010,-(T).010),TV((T).010,(T).010)));
        source_velocity.Append((T).3*TV(-1,-1).Normalized());}
    if(test_number==8){
        TV domain_center=grid.domain.Center();domain_center.y=(T)1;
        sources.Resize(1);source_velocity.Resize(1);world_to_source.Resize(1);
        sources(0)=RANGE<TV>(TV(),TV((T).1,(T).2));
        ARRAY<MATRIX<T,3> > source_to_world(1);
        MATRIX<T,3> translation=MATRIX<T,3>::Translation_Matrix(domain_center);
        source_to_world(0)=translation;world_to_source(0)=source_to_world(0).Inverse();
        source_velocity(0)=TV(0,-.5);}

    // set example-specific parameters
    fluids_parameters.object_friction=(test_number==4||test_number==20)?(T)1:0;
    if(test_number==8){
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        use_two_way_coupling_for_sph=false;
        convert_sph_particles_to_fluid=false;
        use_variable_density_for_sph=true;}
    else if(test_number==9){ 
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        convert_sph_particles_to_fluid=false;
        use_analytic_divergence=true;
        use_two_way_coupling_for_sph=true;
        use_variable_density_for_sph=true;}
    else if(test_number==10){ 
        flip_ratio=1;
        target_particles_per_unit_volume=fluids_parameters.number_particles_per_cell/grid.Cell_Size();
        neumann_boundary_slip_multiplier=0;
        adjust_cell_weights_on_neumann_boundaries=true;
        fluids_parameters.use_sph_for_removed_negative_particles=true;
        convert_sph_particles_to_fluid=false;
        use_analytic_divergence=true;
        use_analytic_divergence_for_expansion_only=false;
        particle_targeting_time=.25;
        use_two_way_coupling_for_sph=false;
        use_variable_density_for_sph=true;
        enforce_density_near_interface=true;}

    if(test_number==4||test_number==5||test_number==6||test_number==20){fluids_parameters.solid_affects_fluid=true;fluids_parameters.fluid_affects_solid=false;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> WATER_STANDARD_TESTS_2D<TV>::
~WATER_STANDARD_TESTS_2D()
{}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Initialize_Advection(const bool always_use_objects)
{
    PHYSBAM_ASSERT(!always_use_objects,"fluids_parameters.solid_affects_fluids overrides this functionality");  // Deprecated parameter
    if(always_use_objects||test_number==4||test_number==5||test_number==6||test_number==20) fluids_parameters.Use_Fluid_Coupling_Defaults();
    else fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class TV> TV WATER_STANDARD_TESTS_2D<TV>::
Initial_Velocity(const TV& X) const
{
    if(test_number==7 || test_number==11) return TV(0,1);
    if(test_number==12 || test_number==13) return TV(-1,0);
    return TV();
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Get_Variable_Viscosity(ARRAY<T,VECTOR<int,2> >& variable_viscosity,const T time) const
{
    if(test_number==13){
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.2)?(T)0:(T)1000;}
    else if(test_number==14){
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.5)?(T)0:(T)1000;}
    else if(test_number==15){
        for(CELL_ITERATOR<TV> iterator(grid,1);iterator.Valid();iterator.Next())
            variable_viscosity(iterator.Cell_Index())=(iterator.Location().x<.5)?(T)0:(T)1000;}
}
//#####################################################################
// Function Initial_Phi
//#####################################################################
template<class TV> typename TV::SCALAR WATER_STANDARD_TESTS_2D<TV>::
Initial_Phi(const TV& X) const
{
    T phi=1;
    if(test_number==1) phi=X.y-(T).32424;
    else if(test_number==2){
        static SPHERE<TV> circle((TV((T).5,(T).75)),(T).2);
        phi=min(circle.Signed_Distance(X),X.y-(T).412134);}
    else if(test_number==4||test_number==20) phi=X.y-(T).400235234;
    else if(test_number==6){
        static SPHERE<TV> circle((TV((T).7,(T).14)),(T).1);
        phi=circle.Signed_Distance(X);}
    else if(test_number==7){
        static SPHERE<TV> circle((TV((T).5,(T).5)),(T).2);
        phi=circle.Signed_Distance(X);}
    else if(test_number==9){
        static RANGE<TV> box((TV((T)3.4,1)),(TV((T)3.6,(T)1.3)));
        phi=min(box.Signed_Distance(X),X.y-(T).4);}
    else if(test_number==11){
        static TV normal((T)1,(T)1),location((T).5,(T).75);
        static LINE_2D<T> line(normal,location);
        phi=line.Signed_Distance(X);}
    else if(test_number==12 || test_number==13){
        static SPHERE<TV> circle((TV((T).5,(T).75)),(T).2);
        phi=circle.Signed_Distance(X);}
    else if(test_number==14){
        static RANGE<TV> box((TV((T).25,(T).75)),(TV((T)0,(T).5)));
        phi=box.Signed_Distance(X);}
    else if(test_number==15){
        static SPHERE<TV> circle((TV((T).5,(T)1)),(T).2);
        phi=circle.Signed_Distance(X);}
    else if(test_number==21) phi=X.y-(T).32424;
    for(int s=0;s<sources.m;s++) phi=min(phi,sources(s).Signed_Distance(world_to_source(s).Homogeneous_Times(X)));
    return phi;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Initialize_Bodies()
{
    if(test_number==4){
        sphere=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies_2D/circle",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T)1.25,(T).55);
        rigid_body_collection.rigid_body_particles.kinematic(sphere)=true;}
    else if(test_number==20){
        sphere=rigid_body_collection.Add_Rigid_Body(example.data_directory+"/Rigid_Bodies_2D/circle",(T).1,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T).8,(T).1);
        rigid_body_collection.rigid_body_particles.kinematic(sphere)=true;}
    if(use_inaccurate_body_collisions){
        inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);}
    else fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
}
//#####################################################################
// Function Update_Sources
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Update_Sources(const T time)
{
    if(test_number==3)
        if(time>4) sources.Clean_Memory();
    if(test_number==5)
        if(time>(T)1.6) sources.Clean_Memory();
    if(test_number==8)
        if(time>2) sources.Clean_Memory();
    if(test_number==10){ //TODO : move this to a better callback function
        FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<TV>&>(fluids_parameters);
        SPH_EVOLUTION_UNIFORM<TV>& sph_evolution=*fluids_parameters_uniform.sph_evolution;
        if(time>12){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume/2;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}
        else if(time>8){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}
        else if(time>4){
            sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume*2;
            sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;}}
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if((test_number==4||test_number==20) && id==sphere) frame.t=motion_curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> bool WATER_STANDARD_TESTS_2D<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if((test_number==4||test_number==20) && id==sphere){twist.linear=motion_curve.Derivative(time);return true;}
    return false;
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Initialize_SPH_Particles()
{
    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<TV>&>(fluids_parameters);
    PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=fluids_parameters_uniform.particle_levelset_evolution->Particle_Levelset(0);
    SPH_EVOLUTION_UNIFORM<TV>& sph_evolution=*fluids_parameters_uniform.sph_evolution;
    GRID<TV>& grid=fluids_parameters_uniform.particle_levelset_evolution->grid;

    sph_evolution.use_two_way_coupling=use_two_way_coupling_for_sph;
    sph_evolution.use_variable_density_solve=use_variable_density_for_sph;
    sph_evolution.convert_particles_to_fluid=convert_sph_particles_to_fluid;
    sph_evolution.use_analytic_divergence=use_analytic_divergence;
    sph_evolution.use_analytic_divergence_for_expansion_only=use_analytic_divergence_for_expansion_only;
    sph_evolution.flip_ratio=flip_ratio;
    sph_evolution.neumann_boundary_slip_multiplier=neumann_boundary_slip_multiplier;
    sph_evolution.adjust_cell_weights_on_neumann_boundaries=adjust_cell_weights_on_neumann_boundaries;
    sph_evolution.enforce_density_near_interface=enforce_density_near_interface;
    sph_evolution.particle_targeting_time=particle_targeting_time;
    sph_evolution.target_particles_per_unit_volume=fluids_parameters.number_particles_per_cell/grid.Cell_Size();
    sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;

    int particle_id=0,number_of_sph_particles=0;

    if(test_number==8){
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        RANGE<TV> particle_region((TV((T).2,1)),TV((T).3,(T)1.2));
        number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume);
        for(int i=0;i<number_of_sph_particles;i++){
            TV X;
            random.Fill_Uniform(X,particle_region);
            TV_INT block=grid.Block_Index(X,3);
            if(!removed_negative_particles(block)) removed_negative_particles(block)=particle_levelset.template_removed_particles.Clone();
            int id=removed_negative_particles(block)->Add_Element();
            (*removed_negative_particles(block)->template Get_Array<int>("id"))(id)=particle_id++;
            removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}}
    else if(test_number==9){
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        RANGE<TV> particle_region((TV((T).4,1)),TV((T).6,(T)1.3));
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
        sph_evolution.target_particles_per_unit_volume=target_particles_per_unit_volume;
        sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*sph_evolution.target_particles_per_unit_volume;
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=particle_levelset.removed_negative_particles;
        RANDOM_NUMBERS<T> random;
        random.Set_Seed(1);
        RANGE<TV> particle_region(TV(grid.domain.min_corner.x,grid.domain.min_corner.y),TV(grid.domain.max_corner.x,grid.domain.max_corner.y/3));
        T particle_multiplier=1;
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
//#####################################################################
// Function Limit_Dt
//#####################################################################
template<class TV> void WATER_STANDARD_TESTS_2D<TV>::
Limit_Dt(T& dt,const T time)
{
    if(test_number==4||test_number==20){
        TV velocity=rigid_body_collection.rigid_body_particles.twist(sphere).linear;
        T rigid_dt_denominator=abs(velocity.x)/grid.dX.x+abs(velocity.y)/grid.dX.y;
        if(rigid_dt_denominator>1e-8) dt=min(dt,1/rigid_dt_denominator);}
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
template<class TV> typename TV::SCALAR WATER_STANDARD_TESTS_2D<TV>::
Initial_Phi_Object(const TV& X) const
{
    if(test_number==4||test_number==20) return rigid_body_collection.Rigid_Body(sphere).Implicit_Geometry_Extended_Value(X);
    else if(test_number==5){
        static RANGE<TV> glass(TV((T).005,(T)-.1),TV((T).095,1));
        return -glass.Signed_Distance(X);}
    else if(test_number==6){
        static LINE_2D<T> line((TV((T)-.2,1)).Normalized(),TV());
        return line.Signed_Distance(X);}
    return 1;
}
//#####################################################################
namespace PhysBAM{
template class WATER_STANDARD_TESTS_2D<VECTOR<float,2> >;
template class WATER_STANDARD_TESTS_2D<VECTOR<double,2> >;
}
