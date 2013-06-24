//#####################################################################
// Copyright 2006-2007, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DENSITY_TARGETTING
//#####################################################################
#ifndef __DENSITY_TARGETTING__
#define __DENSITY_TARGETTING__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Images/IMAGE.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SPH_CALLBACKS.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
namespace PhysBAM{

template<class T_input>
class DENSITY_TARGETTING:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
public:
    typedef VECTOR<T,2> TV;typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef VECTOR<int,2> TV_INT;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;using BASE::fluid_collection;
    using BASE::solid_body_collection;using BASE::stream_type;using BASE::parse_args;using BASE::test_number;using BASE::resolution;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    FLUID_COLLISION_BODY_INACCURATE_UNION<GRID<TV> > inaccurate_union;
    int sphere;
    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > target_image,target_image2;
    GRID<TV> grid_image,grid_image2;

    RANDOM_NUMBERS<T> random;
    ARRAY<ORIENTED_BOX<TV> > sph_sources,initial_sources;
    ARRAY<TV> sph_sources_velocity,initial_sources_velocity;
    ARRAY<RANGE<TV> > sph_sources_bounding_box,initial_sources_bounding_box;

    ORIENTED_BOX<TV> source1,initial_source1;
    TV source1_velocity,initial_source1_velocity;
    RANGE<TV> source1_bounding_box,initial_source1_bounding_box;

    int particle_id;
    T fraction_of_particles_for_sph;

    T time_pour_start,time_pour_end;
    T targetting_time_start,targetting_falloff_time,targetting_time_end,targetting_switch_time;

    T target_density_factor_outside_image,target_density_factor_low,target_density_factor_high;
    
    TV target_translation;

    /***************
    example explanation:
    1. Targetting SIGGRAPH logo 
    2. Low density region surrounded by high density region in zero gravity and zero initial velocity
    3. Same as 2, but with a non-zero horizontal initial velocity
    4. Same as 1, followed by a blending into Eurographics logo.
    5. Adaptivity: A rigid sphere in the middle and low density targetting below it.
    6. Adaptivity: low density targetting closer to floor.
    ***************/

    DENSITY_TARGETTING(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,1,fluids_parameters.WATER),
        rigid_body_collection(solid_body_collection.rigid_body_collection),inaccurate_union(*fluids_parameters.grid)
    {
    }

    ~DENSITY_TARGETTING()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    int cells=1*resolution;
    frame_rate=24;
    restart=false;restart_frame=0;
    fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
    fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
    fluids_parameters.number_particles_per_cell=4;
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
    fluids_parameters.write_ghost_values=true;
    fluids_parameters.store_particle_ids=true;
    first_frame=0;last_frame=1000;

    if(test_number==1 || test_number==4 || test_number==5 || test_number==6 || test_number==7) fluids_parameters.grid->Initialize(TV_INT(6*cells+1,12*cells+1),RANGE<TV>(TV(0,0),TV(1.5,3)));
    else if(test_number==2 || test_number==3){
        frame_rate=100;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.grid->Initialize(TV_INT(12*cells+1,6*cells+1),RANGE<TV>(TV(0,0),TV(6,3)));}
    else{
        LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

    GRID<TV>& grid=*fluids_parameters.grid;

    output_directory=STRING_UTILITIES::string_sprintf("Density_Targetting/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));
        
    fluids_parameters.use_sph_for_removed_negative_particles=true;
    particle_id=0;
    fraction_of_particles_for_sph=1;

    random.Set_Seed(1);
    source1=ORIENTED_BOX<TV>(RANGE<TV>(TV((T).02,(T)2.25),TV((T).12,(T)2.75)),ROTATION<TV>());
    source1_velocity=source1.edges.Column(2).Normalized()*(T)-.05;
    source1_bounding_box=source1.Axis_Aligned_Bounding_Box();

    sph_sources.Append(source1);
    sph_sources_velocity.Append(source1_velocity);
    sph_sources_bounding_box.Append(source1_bounding_box);
        
    initial_source1=ORIENTED_BOX<TV>(RANGE<TV>(TV((T)1.5,1),TV((T)2.5,2)),ROTATION<TV>());
    if(test_number==6 || test_number==7) initial_source1=ORIENTED_BOX<TV>(RANGE<TV>(TV(),TV(grid.domain.max_corner.x,grid.domain.max_corner.y/(T)1.5)),ROTATION<TV>());
    if(test_number==2) initial_source1_velocity=TV();
    else if (test_number==3) initial_source1_velocity=TV(1,0);
    else if (test_number==6 || test_number==7) initial_source1_velocity=TV();
    initial_source1_bounding_box=initial_source1.Axis_Aligned_Bounding_Box();

    if(test_number==1){
        time_pour_start=2;
        time_pour_end=(T)12.5;
        IMAGE<T>::Read("Density_Targetting/siggraphlogobw.jpg",target_image);
        IMAGE<T>::Invert(target_image);
        IMAGE<T>::Threshold(target_image,(T).1,VECTOR<T,3>(0,0,0),VECTOR<T,3>(1,1,1));
        grid_image.Initialize(target_image.Size(),grid.domain);

        target_translation=TV();
        targetting_time_start=time_pour_start;
        targetting_falloff_time=2;
        targetting_time_end=time_pour_end+4;
        targetting_switch_time=targetting_time_end+1000;

        target_density_factor_outside_image=1;
        target_density_factor_low=1;
        target_density_factor_high=5;}
    else if(test_number==2 || test_number==3){
        initial_sources.Append(initial_source1);
        initial_sources_velocity.Append(initial_source1_velocity);
        initial_sources_bounding_box.Append(initial_source1_bounding_box);

        fluids_parameters.gravity=(T)0;
        time_pour_end=-1;
        IMAGE<T>::Read("Density_Targetting/circle.jpg",target_image);
        grid_image.Initialize(target_image.Size(),RANGE<TV>(TV((T)1.5,1),TV((T)2.5,2)));
            
        target_translation=TV();
        targetting_time_start=time_pour_start;
        targetting_falloff_time=2;
        targetting_time_end=time_pour_end+4;
        targetting_switch_time=targetting_time_end+1000;

        target_density_factor_outside_image=1;
        target_density_factor_low=(T)0.1;
        target_density_factor_high=1;}
    else if(test_number==4){
        time_pour_start=2;
        time_pour_end=16;
        IMAGE<T>::Read("Density_Targetting/siggraphlogobw.jpg",target_image);
        IMAGE<T>::Invert(target_image);
        IMAGE<T>::Threshold(target_image,(T).1,VECTOR<T,3>(0,0,0),VECTOR<T,3>(1,1,1));
        grid_image.Initialize(target_image.Size(),grid.domain);

        IMAGE<T>::Read("Density_Targetting/eg.jpg",target_image2);
        IMAGE<T>::Invert(target_image2);
        IMAGE<T>::Threshold(target_image2,(T).1,VECTOR<T,3>(0,0,0),VECTOR<T,3>(1,1,1));
        grid_image2.Initialize(target_image2.Size(),grid.domain);

        target_translation=TV();
        targetting_time_start=time_pour_start;
        targetting_falloff_time=2;
        targetting_switch_time=time_pour_end+1;
        targetting_time_end=targetting_switch_time+6;

        target_density_factor_outside_image=1;
        target_density_factor_low=1;
        target_density_factor_high=5;}
    else if(test_number==5){
        time_pour_start=2;
        time_pour_end=16;
        IMAGE<T>::Read("Density_Targetting/ieee.jpg",target_image);
        IMAGE<T>::Invert(target_image);
        IMAGE<T>::Threshold(target_image,(T).1,VECTOR<T,3>(0,0,0),VECTOR<T,3>(1,1,1));
        grid_image.Initialize(target_image.Size(),grid.domain);

        IMAGE<T>::Read("Density_Targetting/tvcg.jpg",target_image2);
        IMAGE<T>::Invert(target_image2);
        IMAGE<T>::Threshold(target_image2,(T).1,VECTOR<T,3>(0,0,0),VECTOR<T,3>(1,1,1));
        grid_image2.Initialize(target_image2.Size(),grid.domain);

        target_translation=TV();
        targetting_time_start=time_pour_start;
        targetting_falloff_time=2;
        targetting_switch_time=time_pour_end+1;
        targetting_time_end=targetting_switch_time+6;

        target_density_factor_outside_image=1;
        target_density_factor_low=1;
        target_density_factor_high=5;}
    else if(test_number==6){
        initial_sources.Append(initial_source1);
        initial_sources_velocity.Append(initial_source1_velocity);
        initial_sources_bounding_box.Append(initial_source1_bounding_box);

        time_pour_start=2;
        time_pour_end=(T)12.5;

        target_translation=TV();
        targetting_time_start=0;
        targetting_falloff_time=0;
        targetting_time_end=time_pour_end+1000;
        targetting_switch_time=targetting_time_end+1000;

        target_density_factor_outside_image=1;
        target_density_factor_low=(T).2;
        target_density_factor_high=1;}
    else if(test_number==7){
        initial_sources.Append(initial_source1);
        initial_sources_velocity.Append(initial_source1_velocity);
        initial_sources_bounding_box.Append(initial_source1_bounding_box);

        time_pour_start=2;
        time_pour_end=(T)12.5;

        target_translation=TV();
        targetting_time_start=0;
        targetting_falloff_time=0;
        targetting_time_end=time_pour_end+1000;
        targetting_switch_time=targetting_time_end+1000;

        target_density_factor_outside_image=1;
        target_density_factor_low=(T).5;
        target_density_factor_high=1;}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    if(test_number==6) fluids_parameters.Use_Fluid_Coupling_Defaults();
    else fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR<TV> iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        fluid_collection.incompressible_fluid_collection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=T();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=1;
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const TV& X) const
{
    if(test_number==6) return rigid_body_collection.Rigid_Body(sphere).Implicit_Geometry_Extended_Value(X);
    else return 1;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    if(test_number==6){
        sphere=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies_2D/circle",(T).3,true,true,false);
        rigid_body_collection.rigid_body_particles.frame(sphere).t=TV((T).75,(T)1.05);
        rigid_body_collection.Rigid_Body(sphere).Is_Kinematic()=true;

        inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);}
    else fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Construct_Levelsets_For_Objects(time);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    BASE::Set_Kinematic_Positions(frame,time,id);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    return BASE::Set_Kinematic_Velocities(twist,time,id);
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() PHYSBAM_OVERRIDE
{
    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >&>(fluids_parameters);
    SPH_EVOLUTION_UNIFORM<GRID<TV> >& sph_evolution=*fluids_parameters_uniform.sph_evolution;
    GRID<TV> & grid=fluids_parameters_uniform.particle_levelset_evolution->grid;
    sph_evolution.target_particles_per_unit_volume=fluids_parameters.number_particles_per_cell/fluids_parameters.grid->Cell_Size();
    sph_evolution.flip_ratio=1;
    sph_evolution.convert_particles_to_fluid=false;
    sph_evolution.use_analytic_divergence=true;
    sph_evolution.use_analytic_divergence_for_expansion_only=false;
    sph_evolution.particle_targeting_time=(T).25;
    sph_evolution.use_two_way_coupling=false;
    sph_evolution.use_variable_density_solve=false;
    sph_evolution.adjust_cell_weights_on_neumann_boundaries=true;
    sph_evolution.neumann_boundary_slip_multiplier=0;
    sph_evolution.enforce_density_near_interface=true;
    sph_evolution.Set_SPH_Callbacks(*this);

    if(test_number==2 || test_number==3){
        Add_SPH_Particles_For_Sources(initial_sources,initial_sources_bounding_box,initial_sources_velocity,true);
        targetting_time_end=-1;
        return;}
    if(test_number==6 || test_number==7){
        Add_SPH_Particles_For_Sources(initial_sources,initial_sources_bounding_box,initial_sources_velocity,true);
        return;}

    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=fluids_parameters_uniform.particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles;
    T height_scaling=8;
    RANGE<TV> particle_region(TV(grid.domain.min_corner.x,grid.domain.min_corner.y),TV(grid.domain.max_corner.x,grid.domain.max_corner.y/height_scaling));
    int number_of_sph_particles=0;
    T particle_multiplier=1;
    number_of_sph_particles=int(particle_region.Size()*sph_evolution.target_particles_per_unit_volume/particle_multiplier);
    for(int i=0;i<number_of_sph_particles;i++){
        TV X=random.Get_Uniform_Vector(particle_region);
        TV_INT block=grid.Block_Index(X,3);
        if(!removed_negative_particles(block)) removed_negative_particles(block)=new PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>();
        int id=removed_negative_particles(block)->Add_Element();
        (*removed_negative_particles(block)->template Get_Array<int>(ATTRIBUTE_ID_ID))(id)=particle_id++;
        removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}
}
//#####################################################################
// Function Add_SPH_Particles_For_Sources
//#####################################################################
void Add_SPH_Particles_For_Sources(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(time<time_pour_start || time>time_pour_end) return;
    Add_SPH_Particles_For_Sources(sph_sources,sph_sources_bounding_box,sph_sources_velocity);
}
//#####################################################################
// Function Add_SPH_Particles_For_Sources
//#####################################################################
void Add_SPH_Particles_For_Sources(const ARRAY<ORIENTED_BOX<TV> > &sources,const ARRAY<RANGE<TV> > sources_bounding_box,const ARRAY<TV> sources_velocity,
        bool use_initial_targetting=false)
{
    GRID<TV> & grid=fluids_parameters.particle_levelset_evolution->grid;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles;
    for(int s=0;s<sources.m;s++){
        RANGE<TV> source_bounding_box=sources_bounding_box(s);
        RANGE<TV_INT> source_bounding_box_int=RANGE<TV_INT>(grid.Block_Index(TV(source_bounding_box.min_corner.x,source_bounding_box.min_corner.y),1),
                grid.Block_Index(TV(source_bounding_box.max_corner.x,source_bounding_box.max_corner.y),1));
        for(NODE_ITERATOR<TV> node_iterator(grid,source_bounding_box_int);node_iterator.Valid();node_iterator.Next()){
            TV location=node_iterator.Location();TV_INT block=grid.Block_Index(location,1);
            T target_density_factor=1;
            if(use_initial_targetting) target_density_factor=Target_Density_Factor(location,targetting_time_start+targetting_falloff_time);
            BLOCK_UNIFORM<GRID<TV> > block_uniform(grid,block);
            RANGE<TV> block_bounding_box=block_uniform.Bounding_Box();
            if(sources(s).Lazy_Inside(location)){
                if(!removed_negative_particles(block)) removed_negative_particles(block)=new PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>();
                while(removed_negative_particles(block)->Size()<target_density_factor*fraction_of_particles_for_sph*fluids_parameters.number_particles_per_cell){
                    TV X=random.Get_Uniform_Vector(block_bounding_box);
                    int id=removed_negative_particles(block)->Add_Element();
                    (*removed_negative_particles(block)->template Get_Array<int>(ATTRIBUTE_ID_ID))(id)=particle_id++;
                    removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->V(id)=sources_velocity(s);
                    removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}}}}
}
//#####################################################################
// Function Create_Levelset_From_Particles
//#####################################################################
void Create_Levelset_For_SPH(LEVELSET<TV> &levelset) const
{
    GRID<TV>& grid=*fluids_parameters.grid;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) levelset.phi(iterator.Cell_Index())=fluids_parameters.sph_evolution->cell_weight(iterator.Cell_Index())>0?1:-1;
    levelset.Fast_Marching_Method();
}
//#####################################################################
// Function Target_Density_Factor
//#####################################################################
T Target_Density_Factor(const TV& location,const T time) const PHYSBAM_OVERRIDE
{
    if(test_number==7){
        GRID<TV>& grid=*fluids_parameters.grid;
        if(location.y>grid.domain.max_corner.y/(T)3.5) return target_density_factor_high;
        return target_density_factor_low;
    }
    else if(test_number==6){
        if(location.x<(T)1.05&&location.x>(T).45&&location.y<(T)1.05) return target_density_factor_low;
        else return target_density_factor_high;}

    TV location_transformed=location-target_translation;

    if(time<targetting_time_start || time>targetting_time_end) return 1;
    if(!grid_image.Domain().Lazy_Inside(location_transformed)) return target_density_factor_outside_image;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> > interpolation;
    VECTOR<T,3> color=interpolation.Clamped_To_Array(grid_image,target_image,location_transformed);
    if(time>targetting_switch_time){
        LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,VECTOR<T,3> > interpolation2;
        VECTOR<T,3> color2=interpolation2.Clamped_To_Array(grid_image2,target_image2,location_transformed);
        if(time<targetting_switch_time+targetting_falloff_time) color=(1-(time-targetting_switch_time)/targetting_falloff_time)*color+((time-targetting_switch_time)/targetting_falloff_time)*color2;
        else color=color2;}

    T target_density_alpha,target_density;
    if(time<targetting_time_start+targetting_falloff_time) target_density_alpha=(time-targetting_time_start)/targetting_falloff_time;
    else target_density_alpha=1;

    target_density=(1-color.Magnitude()*(T)root_three)*target_density_factor_low+color.Magnitude()*(T)root_three*target_density_factor_high;
    return 1+target_density_alpha*(target_density-1);
}
//#####################################################################
};
}
#endif
