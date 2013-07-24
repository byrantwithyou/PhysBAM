//#####################################################################
// Copyright 2006, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GLASS
//#####################################################################
#ifndef __GLASS__
#define __GLASS__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Images/IMAGE.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SPH_CALLBACKS.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
namespace PhysBAM{

template<class T_input>
class GLASS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
public:
    typedef VECTOR<T,2> TV;typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef VECTOR<int,2> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;using BASE::fluid_collection;
    using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;using BASE::Adjust_Phi_With_Source;

    ARRAY<VECTOR<T,3> ,VECTOR<int,2> > target_image;
    GRID<TV> grid_image;

    RANDOM_NUMBERS<T> random;
    ARRAY<RANGE<TV> > sources;
    ARRAY<ORIENTED_BOX<TV> > sph_sources;
    ARRAY<MATRIX<T,3> > world_to_source;
    ARRAY<VECTOR<T,2> > sph_sources_velocity,sources_velocity;
    ARRAY<RANGE<TV> > sph_sources_bounding_box,sources_bounding_box;

    RANGE<TV> source1,source4;
    ORIENTED_BOX<TV> source2,source3;
    VECTOR<T,2> source1_velocity,source2_velocity,source3_velocity,source4_velocity;
    RANGE<TV> source1_bounding_box,source2_bounding_box,source3_bounding_box,source4_bounding_box;
    MATRIX<T,3> world_to_source1,world_to_source4;

    T time_pour_start_sph,time_pour_end_sph,time_pour_end_pls;
    int particle_id;

    bool use_variable_density_for_sph,use_two_way_coupling_for_sph,convert_sph_particles_to_fluid;
    bool use_analytic_divergence,use_analytic_divergence_for_expansion_only;
    bool adjust_cell_weights_on_neumann_boundaries;
    T ballistic_particles_as_percentage_of_target,particle_targeting_time,fraction_of_particles_for_sph,flip_ratio;

    /***************
    example explanation:
    1. SPH source
    2. Levelset source with OWC VD SPH for NRPs
    3. SPH source and levelset source with TWC VD SPH
    ***************/

    GLASS(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,1,fluids_parameters.WATER)
    {
    }

    ~GLASS()
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
    fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[3][0]=true;fluids_parameters.domain_walls[3][1]=true;
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
    GRID<TV>& grid=*fluids_parameters.grid;

    grid.Initialize(TV_INT(10*cells+1,20*cells+1),RANGE<TV>(TV((T)-.05,0),TV((T).05,(T).2)));
    //grid.Initialize(TV_INT(15*cells+1,20*cells+1),RANGE<TV>(TV((T)-.075,0),TV((T).075,(T).2)));

    output_directory=STRING_UTILITIES::string_sprintf("Glass/Test_%d__Resolution_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1));
        
    random.Set_Seed(1);
    particle_id=0;
    // SPH paramters
    fluids_parameters.use_sph_for_removed_negative_particles=true;
    use_variable_density_for_sph=true;
    use_two_way_coupling_for_sph=true;
    convert_sph_particles_to_fluid=false;
    use_analytic_divergence=false;
    use_analytic_divergence_for_expansion_only=false;
    adjust_cell_weights_on_neumann_boundaries=false;
    ballistic_particles_as_percentage_of_target=(T).1;
    particle_targeting_time=(T).25;
    flip_ratio=1;
    fraction_of_particles_for_sph=1;

    if(test_number==2){
        convert_sph_particles_to_fluid=true;
        use_two_way_coupling_for_sph=true;
        use_variable_density_for_sph=false;
        flip_ratio=(T).2;
        use_analytic_divergence=true;
        use_analytic_divergence_for_expansion_only=false;}
    if(test_number==3){
        use_analytic_divergence=true;
        use_analytic_divergence_for_expansion_only=false;
        use_two_way_coupling_for_sph=true;
        fraction_of_particles_for_sph=(T).8;}


    time_pour_start_sph=(T).7;
    time_pour_end_sph=1000;
    time_pour_end_pls=1000;
    source1=RANGE<TV>(TV((T)-.01,(T)-.0075),TV((T).01,(T).0075));
    MATRIX<T,3> rotation=MATRIX<T,3>::Rotation_Matrix_Z_Axis(-(T)pi/4);
    MATRIX<T,3> translation=MATRIX<T,3>::Translation_Matrix(TV((T)-.01,(T).175));
    MATRIX<T,3> source1_to_world=translation*rotation;
    world_to_source1=source1_to_world.Inverse();
    source1_velocity=source1_to_world.Extract_Rotation()*TV(0,-1)*(T).5;
    source1_bounding_box=source1;

    rotation=MATRIX<T,3>::Rotation_Matrix_Z_Axis(-(T)pi/4);
    translation=MATRIX<T,3>::Translation_Matrix(TV((T)-.01,(T).175));
    MATRIX<T,3> source2_to_world=translation*rotation;
    source2=ORIENTED_BOX<TV>(source1,source2_to_world);
    source2_velocity=source2.edges.Column(2).Normalized()*(T)-.5;
    source2_bounding_box=source2.Axis_Aligned_Bounding_Box();

    rotation=MATRIX<T,3>::Rotation_Matrix_Z_Axis(-(T)pi/3);
    translation=MATRIX<T,3>::Translation_Matrix(TV((T)-.03,(T).175));
    MATRIX<T,3> source3_to_world=translation*rotation;
    source3=ORIENTED_BOX<TV>(source1,source3_to_world);
    source3_velocity=source3.edges.Column(2).Normalized()*(T)-.2;
    source3_bounding_box=source3.Axis_Aligned_Bounding_Box();

    source4=source1;
    rotation=MATRIX<T,3>::Rotation_Matrix_Z_Axis((T)pi/3);
    translation=MATRIX<T,3>::Translation_Matrix(TV((T).03,(T).175));
    MATRIX<T,3> source4_to_world=translation*rotation;
    world_to_source4=source4_to_world.Inverse();
    source4_velocity=source4_to_world.Extract_Rotation()*TV(0,-1)*(T).3;
    source4_bounding_box=source4;

    if(test_number==1){
        sources.Append(source1);
        world_to_source.Append(world_to_source1);
        sources_velocity.Append(source1_velocity);
        sources_bounding_box.Append(source1_bounding_box);}
    else if(test_number==2){
        sph_sources.Append(source2);
        sph_sources_velocity.Append(source2_velocity);
        sph_sources_bounding_box.Append(source2_bounding_box);}
    else if(test_number==3){
        frame_rate=192;
        time_pour_end_pls=(T)1.6;
        time_pour_end_sph=(T)1.6;
        time_pour_end_pls=(T)1.12;
        time_pour_end_sph=(T)1.12;
        sph_sources.Append(source3);
        sph_sources_velocity.Append(source3_velocity);
        sph_sources_bounding_box.Append(source3_bounding_box);

        sources.Append(source4);
        world_to_source.Append(world_to_source4);
        sources_velocity.Append(source4_velocity);
        sources_bounding_box.Append(source4_bounding_box);}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
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
    T initial_phi;
    LOG::cout<<"setting phi values"<<std::endl;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();
        initial_phi=1;
        for(int s=0;s<sources.m;s++) initial_phi=min(initial_phi,sources(s).Signed_Distance(world_to_source(s).Homogeneous_Times(X)));
        phi(iterator.Cell_Index())=initial_phi;}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_collection);
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
    if(time>time_pour_end_pls) sources.Clean_Memory();
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    for(int s=0;s<sources.m;s++)Adjust_Phi_With_Source(sources(s),world_to_source(s));
    T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& removed_positive_particles=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_positive_particles;
    GRID<TV>& grid=*fluids_parameters.grid;

    for(int s=0;s<sources.m;s++){
        for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
            if(!removed_positive_particles(block)) continue;
            for(int p=0;p<removed_positive_particles(block)->Size();p++)
                if(sources(s).Inside(world_to_source(s).Homogeneous_Times(removed_positive_particles(block)->X(p)),(T)-1e-4)) removed_positive_particles(block)->Add_To_Deletion_List(p);
            removed_positive_particles(block)->Delete_Elements_On_Deletion_List();}}
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
{
    bool first=true;
    for(int s=0;s<sources.m;s++){Get_Source_Reseed_Mask(sources(s),world_to_source(s),cell_centered_mask,first);first=false;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    for(int s=0;s<sources.m;s++)Get_Source_Velocities(sources(s),world_to_source(s),sources_velocity(s));
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
    for(int s=0;s<sources.m;s++) dt=min(dt,fluids_parameters.cfl*fluids_parameters.grid->dX.Min()/sources_velocity(s).Magnitude());
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() PHYSBAM_OVERRIDE
{
    GRID<TV> & grid=*fluids_parameters.grid;
    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_uniform=dynamic_cast<FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >&>(fluids_parameters);
    SPH_EVOLUTION_UNIFORM<GRID<TV> >& sph_evolution=*fluids_parameters_uniform.sph_evolution;
    sph_evolution.use_analytic_divergence_for_expansion_only=use_analytic_divergence_for_expansion_only;
    sph_evolution.use_two_way_coupling=use_two_way_coupling_for_sph;
    sph_evolution.use_variable_density_solve=use_variable_density_for_sph;
    sph_evolution.convert_particles_to_fluid=convert_sph_particles_to_fluid;
    sph_evolution.use_analytic_divergence=use_analytic_divergence;
    sph_evolution.adjust_cell_weights_on_neumann_boundaries=adjust_cell_weights_on_neumann_boundaries;
    sph_evolution.particle_targeting_time=particle_targeting_time;
    sph_evolution.flip_ratio=flip_ratio;
    sph_evolution.target_particles_per_unit_volume=fraction_of_particles_for_sph*fluids_parameters.number_particles_per_cell/grid.Cell_Size();
    sph_evolution.ballistic_particles_per_unit_volume=ballistic_particles_as_percentage_of_target*fluids_parameters.sph_evolution->target_particles_per_unit_volume;
}
//#####################################################################
// Function Add_SPH_Particles_For_Sources
//#####################################################################
void Add_SPH_Particles_For_Sources(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==1) return; 
    if(time<time_pour_start_sph||time>time_pour_end_sph) return;
    Add_SPH_Particles_For_Sources(sph_sources,sph_sources_bounding_box,sph_sources_velocity);
}
//#####################################################################
// Function Add_SPH_Particles_For_Sources
//#####################################################################
void Add_SPH_Particles_For_Sources(const ARRAY<ORIENTED_BOX<TV> > &sph_sources,const ARRAY<RANGE<TV> > sph_sources_bounding_box,const ARRAY<VECTOR<T,2> > sources_velocity,
        bool use_initial_targetting=false)
{
    GRID<TV> & grid=fluids_parameters.particle_levelset_evolution->grid;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,2> >& removed_negative_particles=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles;
    for(int s=0;s<sph_sources.m;s++){
        RANGE<TV> source_bounding_box=sph_sources_bounding_box(s);
        RANGE<TV_INT> source_bounding_box_int=RANGE<TV_INT>(grid.Block_Index(TV(source_bounding_box.min_corner.x,source_bounding_box.min_corner.y),1),
                grid.Block_Index(TV(source_bounding_box.max_corner.x,source_bounding_box.max_corner.y),1));
        for(NODE_ITERATOR<TV> node_iterator(grid,source_bounding_box_int);node_iterator.Valid();node_iterator.Next()){
            TV location=node_iterator.Location();TV_INT block=grid.Block_Index(location,1);
            T target_density_factor=1;
            BLOCK_UNIFORM<GRID<TV> > block_uniform(grid,block);
            RANGE<TV> block_bounding_box=block_uniform.Bounding_Box();
            if(sph_sources(s).Lazy_Inside(location)){
                if(!removed_negative_particles(block)) removed_negative_particles(block)=new PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>();
                while(removed_negative_particles(block)->Size()<target_density_factor*fraction_of_particles_for_sph*fluids_parameters.number_particles_per_cell){
                    TV X=random.Get_Uniform_Vector(block_bounding_box);
                    int id=removed_negative_particles(block)->Add_Element();
                    (*removed_negative_particles(block)->template Get_Array<int>(ATTRIBUTE_ID_ID))(id)=particle_id++;
                    removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->V(id)=sources_velocity(s);
                    removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}}}}
}
//#####################################################################
};
}
#endif
