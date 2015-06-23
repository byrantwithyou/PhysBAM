//#####################################################################
// Copyright 2006-2007, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GLASS
//#####################################################################
#ifndef __GLASS__
#define __GLASS__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Log/LOG.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Dynamics/Standard_Tests/WATER_STANDARD_TESTS_3D.h>
namespace PhysBAM{

template<class T_input>
class GLASS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;
public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,TV_INT> T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;using BASE::solid_body_collection;
    using BASE::stream_type;using BASE::test_number;using BASE::resolution;using BASE::Adjust_Phi_With_Source;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    FLUID_COLLISION_BODY_INACCURATE_UNION<TV> inaccurate_union;
    int glass;

    int particle_id;
    ARRAY<CYLINDER<T> > sph_sources;
    ARRAY<VECTOR<T,3> > sph_sources_velocity;
    ARRAY<RANGE<TV> > sph_sources_bounding_box;
    RANDOM_NUMBERS<T> random;
    T time_pour;

    ARRAY<CYLINDER<T> > sources;
    ARRAY<MATRIX<T,4> > world_to_sources;
    ARRAY<VECTOR<T,3> > sources_velocity;

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

    GLASS(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>(stream_type_input,parse_args,1,fluids_parameters.WATER),
        rigid_body_collection(solid_body_collection.rigid_body_collection),inaccurate_union(*fluids_parameters.grid),particle_id(0)
    {
        parse_args.Parse();
        random.Set_Seed(1);
        // set up the standard fluid environment
        frame_rate=96;
        restart=false;restart_frame=0;
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
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
        fluids_parameters.cfl=(T)1.9;

        // MacCormack parameters
        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;
        fluids_parameters.use_maccormack_compute_mask=true;
        fluids_parameters.use_maccormack_for_incompressible=true;
        fluids_parameters.bandwidth_without_maccormack_near_interface=1;

        int cells=1*resolution;
        first_frame=0;last_frame=10;
        GRID<TV>& grid=*fluids_parameters.grid;
        grid.Initialize(TV_INT(10*cells+1,20*cells+1,10*cells+1),RANGE<TV>(TV(-(T).05,0,-(T).05),TV((T).05,(T).2,(T).05)));

        last_frame=1000;
        
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
        convert_sph_particles_to_fluid=true;
        flip_ratio=1;
        fraction_of_particles_for_sph=1;

        if(test_number==2){
            fluids_parameters.second_order_cut_cell_method=false;
            convert_sph_particles_to_fluid=true;
            use_two_way_coupling_for_sph=true;
            use_variable_density_for_sph=true;
            flip_ratio=(T)0;
            use_analytic_divergence=true;
            use_analytic_divergence_for_expansion_only=false;}
        if(test_number==3){
            use_analytic_divergence_for_expansion_only=false;
            use_two_way_coupling_for_sph=true;
            fraction_of_particles_for_sph=(T).25;}

        time_pour=4;
        
        output_directory=LOG::sprintf("Glass/Glass_%d__Resolution_%d_%d_%d",test_number,(grid.counts.x-1),(grid.counts.y-1),(grid.counts.z-1));
        LOG::cout<<"Running SPH simulation to "<<output_directory<<std::endl;

        CYLINDER<T> source1;
        source1.Set_Endpoints(TV(-(T).0025,(T).174,0),TV((T).0077,(T).1875,0));
        source1.radius=(T).01;
        VECTOR<T,3> source1_velocity=(T)1.2*source1.plane1.normal;

        CYLINDER<T> source2;
        source2.Set_Endpoints(TV(-(T).04,(T).1875,-(T).01),TV(-(T).02,(T).1725,-(T).01));
        source2.radius=(T).01;
        VECTOR<T,3> source2_velocity=-(T).35*source2.plane1.normal;

        CYLINDER<T> source3;
        source3.Set_Endpoints(TV((T).02,(T).1725,(T).01),TV((T).04,(T).1875,(T).01));
        source3.radius=(T).01;
        VECTOR<T,3> source3_velocity=(T).35*source3.plane1.normal;

        world_to_sources.Append(MATRIX<T,4>::Identity_Matrix());
        if(test_number==1){
            sph_sources.Append(source1);
            sph_sources_velocity.Append(source1_velocity);}
        else if(test_number==2){
            sources.Append(source1);
            sources_velocity.Append(source1_velocity);}
        else if(test_number==3){
            sph_sources.Append(source2);
            sph_sources_velocity.Append(source2_velocity);
            sources.Append(source3);
            sources_velocity.Append(source3_velocity);}
    }

    ~GLASS()
    {}

    // Unused callbacks
    void Initialize_Velocities() override {}
    void Construct_Levelsets_For_Objects(const T time){}
    void Preprocess_Frame(const int frame) override {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Solids_Substep(const T time,const int substep) override {}
    void Postprocess_Frame(const int frame) override {}
    void Postprocess_Phi(const T time) override {}

void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() override
{
    fluids_parameters.Use_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    T initial_phi;
    LOG::cout<<"setting phi values"<<std::endl;
    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV X=iterator.Location();
        initial_phi=1;
        for(int s=0;s<sources.m;s++)initial_phi=min(initial_phi,sources(s).Signed_Distance(world_to_sources(s).Homogeneous_Times(X)));
        phi(iterator.Cell_Index())=initial_phi;}
}
//#####################################################################
// Function Initialize_SPH_Particles
//#####################################################################
void Initialize_SPH_Particles() override
{
    GRID<TV>& grid=*fluids_parameters.grid;
    SPH_EVOLUTION_UNIFORM<TV>& sph_evolution=*fluids_parameters.sph_evolution;
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
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const VECTOR<T,3>& X) const
{
    return rigid_body_collection.Rigid_Body(glass).Implicit_Geometry_Extended_Value(X);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    glass=rigid_body_collection.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/glass_cylinder",(T).03,true,true,false);
    rigid_body_collection.rigid_body_particles.frame(glass).t=VECTOR<T,3>(0,(T)0.06,0);
    rigid_body_collection.Rigid_Body(glass).Is_Kinematic()=true;

    inaccurate_union.collision_bodies.Add_Bodies(rigid_body_collection);
    fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&inaccurate_union,0,false);
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) override
{
    for(int s=0;s<sources.m;s++) dt=min(dt,fluids_parameters.cfl*fluids_parameters.grid->dX.Min()/sources_velocity(s).Magnitude());
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) override
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV>::Update_Fluid_Parameters(dt,time);
    if(time>1.85) sources.Clean_Memory();
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) override
{
    for(int s=0;s<sources.m;s++)Adjust_Phi_With_Source(sources(s),world_to_sources(s));
    
    T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& removed_positive_particles=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_positive_particles;
    GRID<TV>& grid=*fluids_parameters.grid; 

    for(int s=0;s<sources.m;s++){
        RANGE<TV_INT> source_nodes(grid.Clamp_To_Cell(sources(s).Bounding_Box().Minimum_Corner()),grid.Clamp_To_Cell(sources(s).Bounding_Box().Maximum_Corner())+TV_INT::All_Ones_Vector());
        for(NODE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
            if(!removed_positive_particles(block)) continue;
            for(int p=0;p<removed_positive_particles(block)->Size();p++)
                if(sources(s).Inside(removed_positive_particles(block)->X(p),-(T)1e-4)) removed_positive_particles(block)->Add_To_Deletion_List(p);
            removed_positive_particles(block)->Delete_Elements_On_Deletion_List();}}
    return false;
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) override
{
    bool first=true;
    for(int s=0;s<sources.m;s++){Get_Source_Reseed_Mask(sources(s),world_to_sources(s),cell_centered_mask,first);first=false;}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,const T time) override
{
    for(int s=0;s<sources.m;s++)Get_Source_Velocities(sources(s),world_to_sources(s),sources_velocity(s));
}
//#####################################################################
// Function Add_SPH_Particles_For_Sources
//#####################################################################
void Add_SPH_Particles_For_Sources(const T dt,const T time) override
{
    if(time>time_pour) return;
    GRID<TV>& grid=fluids_parameters.particle_levelset_evolution->grid;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> >& removed_negative_particles=fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).removed_negative_particles;
    for(int s=0;s<sph_sources.m;s++){
        RANGE<TV> source_bounding_box=sph_sources(s).Bounding_Box();
        RANGE<TV_INT> source_bounding_box_int=RANGE<TV_INT>(grid.Block_Index(TV(source_bounding_box.min_corner.x,source_bounding_box.min_corner.y,source_bounding_box.min_corner.z),1),
             grid.Block_Index(TV(source_bounding_box.max_corner.x,source_bounding_box.max_corner.y,source_bounding_box.max_corner.z),1));
        for(NODE_ITERATOR<TV> node_iterator(grid,source_bounding_box_int);node_iterator.Valid();node_iterator.Next()){
            TV location=node_iterator.Location();TV_INT block=grid.Block_Index(location,1);
            BLOCK_UNIFORM<TV> block_uniform(grid,block);
            RANGE<TV> block_bounding_box=block_uniform.Bounding_Box();
            if(sph_sources(s).Lazy_Inside(location)){
                if(!removed_negative_particles(block)) removed_negative_particles(block)=new PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>();
                while(removed_negative_particles(block)->Size()<.5*fraction_of_particles_for_sph*fluids_parameters.number_particles_per_cell){
                    TV X=random.Get_Uniform_Vector(block_bounding_box);
                    ARRAY_VIEW<int>& id_attr=*removed_negative_particles(block)->template Get_Array<int>(ATTRIBUTE_ID_ID);
                    int id=removed_negative_particles(block)->Add_Element();
                    id_attr(id)=particle_id++;
                    removed_negative_particles(block)->X(id)=X;removed_negative_particles(block)->V(id)=sph_sources_velocity(s);
                    removed_negative_particles(block)->radius(id)=(T).1*grid.dX.Min();}}}}
}
//#####################################################################
};
}
#endif
