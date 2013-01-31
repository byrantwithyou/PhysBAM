//#####################################################################
// Copyright 2006-2007, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MASS_CONSERVATION 
// Examples:
// 1. Box in analytic cardinal velocity field
// 2. Box in analytic non-axis-aligned velocity field
// 3. Zalesak disk
// 4. Zalesak not at center of grid.
// 5. Vortex test
// 6. Deformation field
// 7. Milk crown
// 8. Disappearing opposed drops
// 9. Multiple drops intercolliding zero-g
// 10. Crossing streams
// 11. Zalesak's migrating disk
// 12. Zalesak's inverted migrating disk
//##################################################################### 
#ifndef __MASS_CONSERVATION__
#define __MASS_CONSERVATION__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_ELEMENTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class RW=T_input>
class MASS_CONSERVATION:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef GRID<TV> T_GRID;typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;typedef VECTOR<int,2> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;typedef UNIFORM_GRID_ITERATOR_FACE<TV> FACE_ITERATOR;typedef UNIFORM_GRID_ITERATOR_NODE<TV> NODE_ITERATOR;
public:
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::fluids_parameters;
    using BASE::fluid_collection;using BASE::parse_args;
    using BASE::solids_parameters;using BASE::write_time;using BASE::write_frame_title;using BASE::data_directory;using BASE::abort_when_dt_below;using BASE::Set_External_Velocities;
    using BASE::resolution;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::test_number; // silence -Woverloaded-virtual

    SPHERE<TV> source,other_source;
    RANGE<TV> box;
    MATRIX<T,3> world_to_source;
    TV source_velocity;
    ARRAY<T,FACE_INDEX<2> > velocity;
    T initial_water_level;
    TV zalesak_center,zalesak_velocity_center;
    RANDOM_NUMBERS<T> random;
    T pi_over_314;
    T root_two_over_two;

    PARAMETER_LIST parameter_list;
    T period;

    MASS_CONSERVATION(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID >(stream_type,1,fluids_parameters.WATER),pi_over_314((T)pi/314),root_two_over_two((T).5*(T)root_two)
    {
    }

    virtual ~MASS_CONSERVATION() 
    {}

    // Unused callbacks
    void Initialize_Bodies() PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Update_Rigid_Bodies(const T time){}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}
    void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}

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
    LOG::cout<<"Example number "<<test_number<<std::endl;
    std::string filename=STRING_UTILITIES::string_sprintf("Mass_Conservation/example_%d.param",test_number);
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading parameter file '"<<filename<<"'"<<std::endl;parameter_list.Read(filename);}
    random.Set_Seed(7411);

    // Common parameters
    output_directory=STRING_UTILITIES::string_sprintf("Mass_Conservation/example_%d_resolution_%d",test_number%resolution);
    first_frame=0;restart=false;restart_frame=0;frame_rate=300;

    // Fluids parameters
    fluids_parameters.number_particles_per_cell=16;
    fluids_parameters.reseeding_frame_rate=last_frame/5;
    fluids_parameters.incompressible_iterations=50;
    fluids_parameters.cfl=(T).5;
    fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][0]=fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.bias_towards_negative_particles=false;
    fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=true;
    fluids_parameters.write_flattened_particles=true;        
    fluids_parameters.write_removed_negative_particles=true;
    fluids_parameters.store_particle_ids=true;
    fluids_parameters.delete_fluid_inside_objects=false;
    fluids_parameters.gravity=0;
    fluids_parameters.enforce_divergence_free_extrapolation=false;
    fluids_parameters.incompressible_tolerance=(T)1e-12;
    fluids_parameters.use_maccormack_for_level_set=false;
    fluids_parameters.use_maccormack_for_incompressible=false;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;

    period=8;
    T final_time=period;
    last_frame=628;
    world_to_source=MATRIX<T,3>::Identity_Matrix();

    if(test_number==1 || test_number==2){
        initial_water_level=-1; // no initial height
        box=RANGE<TV>((T).4,(T).6,(T).4,(T).6);
        fluids_parameters.grid->Initialize(TV_INT()+10*resolution+1,RANGE<TV>::Unit_Box());}
    else if(test_number==3){
        final_time=15;frame_rate=(T).5;
        zalesak_center=TV(50,75);
        zalesak_velocity_center=TV(50,50);
        fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,10*resolution+1),RANGE<TV>(TV(0,0),TV(100,100)));}
    else if(test_number==4){
        final_time=628;frame_rate=(T).5;
        zalesak_center=TV(25,25);
        zalesak_velocity_center=TV(50,75);
        fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,10*resolution+1),RANGE<TV>(TV(0,0),TV(100,100)));}
    else if(test_number==5 || test_number==13){
        fluids_parameters.grid->Initialize(TV_INT()+(1<<resolution)+1,RANGE<TV>::Unit_Box());
        initial_water_level=-1;
        fluids_parameters.cfl=(T).5;
        last_frame=100;
        frame_rate=(T)last_frame/final_time;
        source=SPHERE<TV>(TV((T).5,(T).75),(T).15);}
    else if(test_number==6){
        fluids_parameters.number_particles_per_cell=128;
        period=final_time=2;
        fluids_parameters.grid->Initialize(TV_INT((1<<resolution)+1,3*(1<<resolution)+1),RANGE<TV>(TV((T)0,(T)-1),TV((T)1,(T)2)));
        initial_water_level=-2;
        fluids_parameters.cfl=(T).5;
        last_frame=100;
        frame_rate=(T)last_frame/final_time;
        source=SPHERE<TV>(TV((T).5,(T).75),(T).15);}
    else if(test_number==7){
        frame_rate=500;
        fluids_parameters.grid->Initialize(TV_INT(10*resolution,(int)(7.5*resolution)),RANGE<TV>(TV((T)0,(T)0),TV((T)2,(T)1.5)));
        source_velocity=TV(0,-40);
        initial_water_level=(T).02;
        source=SPHERE<TV>(TV(1,initial_water_level+(T).15),(T).05);
        fluids_parameters.surface_tension=(T)5e-6;}
    else if(test_number==8){
        frame_rate=100;
        fluids_parameters.grid->Initialize(TV_INT()+10*resolution+1,RANGE<TV>::Unit_Box());
        fluids_parameters.gravity=0;
        source=SPHERE<TV>(TV((T).1,(T).5),(T).08);
        other_source=SPHERE<TV>(TV((T).9,(T).5),(T).08);
        source_velocity=TV((T)2,(T)0);
        initial_water_level=-1;}
    else if(test_number==9){
        fluids_parameters.grid->Initialize(TV_INT()+10*resolution+1,RANGE<TV>::Unit_Box());
        fluids_parameters.gravity=0;
        fluids_parameters.surface_tension=(T)1e-5;}
    else if(test_number==10){
        fluids_parameters.grid->Initialize(TV_INT()+10*resolution+1,RANGE<TV>::Unit_Box());}
    else if(test_number==11 || test_number==12){
        final_time=628;frame_rate=(T).5;
        zalesak_center=TV(25,25);
        zalesak_velocity_center=TV(50,50);
        fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,10*resolution+1),RANGE<TV>(TV(0,0),TV(100,100)));}

    // Debugging
    fluids_parameters.write_debug_data=true;
    abort_when_dt_below=(T)1e-7;
    write_time=true;
    write_frame_title=true;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
    fluids_parameters.particle_levelset_evolution->particle_levelset.reincorporate_removed_particles_everywhere=false;
    if(test_number == 1 || test_number == 2 || test_number==3 || test_number==4 || test_number==5 || test_number==6 || test_number==11 || test_number==12 || test_number==13){
        fluids_parameters.particle_levelset_evolution->levelset_advection.Use_Local_WENO_For_Advection();
        fluids_parameters.particle_levelset_evolution->Set_Runge_Kutta_Order_Levelset(3);
        fluids_parameters.particle_levelset_evolution->Set_Runge_Kutta_Order_Particles(3);
        fluids_parameters.analytic_test=true;
        fluids_parameters.reseeding_frame_rate=200000;
        fluids_parameters.reinitialize_geometry_frame_rate=200000;
        fluids_parameters.particle_levelset_evolution->Use_Reinitialization();
        fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).reincorporate_removed_particles_everywhere=true;}
    if(test_number==5 || test_number==6 || test_number==13) fluids_parameters.particle_levelset_evolution->Use_Frozen_Velocity(false);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if((test_number==11 || test_number==12) && frame && frame%40==0){
        T_GRID& grid=*fluids_parameters.grid;
        T_ARRAYS_SCALAR copy_phis=fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi;
        T_ARRAYS_SCALAR copy_volumes;
        fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi.Fill((T)1);
        TV position_offset;
        if(test_number==11){
            TV_INT offset=TV_INT(0,50);
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(cell_index.y<grid.counts.x-50){
                    fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi(cell_index)=copy_phis(cell_index+offset);}}
            position_offset=TV((T)0,(T)-50);
            for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT block=iterator.Node_Index();
                if(fluids_parameters.particle_levelset_evolution->particle_levelset.negative_particles(block)){
                    PARTICLE_LEVELSET_PARTICLES<TV>& block_particles=*fluids_parameters.particle_levelset_evolution->particle_levelset.negative_particles(block);
                    for(int k=0;k<block_particles.Size();k++) block_particles.X(k)+=position_offset;}
                if(fluids_parameters.particle_levelset_evolution->particle_levelset.positive_particles(block)){
                    PARTICLE_LEVELSET_PARTICLES<TV>& block_particles=*fluids_parameters.particle_levelset_evolution->particle_levelset.positive_particles(block);
                    for(int k=0;k<block_particles.Size();k++) block_particles.X(k)+=position_offset;}}}
        else if(test_number==12){
            TV_INT offset=TV_INT(0,-50);
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(cell_index.y>50){
                    fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi(cell_index)=copy_phis(cell_index+offset);}}
            position_offset=TV((T)0,(T)50);
            for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT block=iterator.Node_Index();
                if(fluids_parameters.particle_levelset_evolution->particle_levelset.negative_particles(block)){
                    PARTICLE_LEVELSET_PARTICLES<TV>& block_particles=*fluids_parameters.particle_levelset_evolution->particle_levelset.negative_particles(block);
                    for(int k=0;k<block_particles.Size();k++) block_particles.X(k)+=position_offset;}
                if(fluids_parameters.particle_levelset_evolution->particle_levelset.positive_particles(block)){
                    PARTICLE_LEVELSET_PARTICLES<TV>& block_particles=*fluids_parameters.particle_levelset_evolution->particle_levelset.positive_particles(block);
                    for(int k=0;k<block_particles.Size();k++) block_particles.X(k)+=position_offset;}}}
        
        fluids_parameters.particle_levelset_evolution->particle_levelset.Update_Particle_Cells(fluids_parameters.particle_levelset_evolution->particle_levelset.negative_particles);
        fluids_parameters.particle_levelset_evolution->particle_levelset.Update_Particle_Cells(fluids_parameters.particle_levelset_evolution->particle_levelset.positive_particles);}
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID >::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    // Not so good to set up a heaviside function here because then the interface will be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    T_GRID& grid=*fluids_parameters.grid;

    if(test_number==1 || test_number==2){
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=box.Signed_Distance(location);}}
    else if(test_number==3 || test_number==4 || test_number==11 || test_number==12){
        T radius=15,slot_width=5,slot_depth=25;
        RANGE<TV> rect(zalesak_center.x-slot_width/2,zalesak_center.x+slot_width/2,0,zalesak_center.y-radius+slot_depth);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=max((location-zalesak_center).Magnitude()-radius,-rect.Signed_Distance(location));}}
    else if(test_number==5 || test_number==6 || test_number==7 || test_number==8 || test_number==13)
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(location.y-initial_water_level,source.Signed_Distance(location));}

    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    if(test_number==7)
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(iterator.Location().y>initial_water_level+grid.dX.Min()*2) face_velocities(axis,face)=source_velocity[axis];}
    else if(test_number==8){
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),other_source.Signed_Distance(location));}
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(iterator.Location().x<(T).5) face_velocities(axis,face)=source_velocity[axis];
            else face_velocities(axis,face)=-source_velocity[axis];}}
    else if(test_number==9){
        // list of source circles and velocities with random locations, random radii within a range, random velocities within a range.  Overlap is fine.
        int drops=40;
        T radius_min=(T).02,radius_max=(T).03;
        ARRAY<SPHERE<TV> > sources(drops);
        ARRAY<TV> velocities(drops);
        RANGE<TV> velocity_box((T)-.2,(T).2,(T)-.2,(T).2);
        for(int i=0;i<drops;i++){
            sources(i).center=random.Get_Uniform_Vector(grid.Domain());
            sources(i).radius=random.Get_Uniform_Number(radius_min,radius_max);
            velocities(i)=random.Get_Uniform_Vector(velocity_box);}
         for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();TV_INT cell=iterator.Cell_Index();
            T& phi=fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index());
            phi=FLT_MAX;
            for(int i=0;i<drops;i++) phi=min(phi,sources(i).Signed_Distance(location));}
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
            T phi=FLT_MAX;
            for(int i=0;i<drops;i++){
                T phi_source=sources(i).Signed_Distance(location);
                if(phi_source<phi){
                    phi=phi_source;
                    face_velocities(axis,face)=velocities(i)(axis);}}}}
    else if(test_number==10){
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=1;}}
    
    fluids_parameters.particle_levelset_evolution->Reseed_Particles(0);
    fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.refine_fmm_initialization_with_iterative_solver=true;
    if(test_number==12){
        // invert everything around y=50
        T_ARRAYS_SCALAR copy_phis=fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi;
        // fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi.Fill((T)1);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            TV_INT mirror_index=cell_index;
            mirror_index.y=101-mirror_index.y;
            fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi(cell_index)=copy_phis(mirror_index);}
        for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT block=iterator.Node_Index();
            if(fluids_parameters.particle_levelset_evolution->particle_levelset.negative_particles(block)){
                PARTICLE_LEVELSET_PARTICLES<TV>& block_particles=*fluids_parameters.particle_levelset_evolution->particle_levelset.negative_particles(block);
                for(int k=0;k<block_particles.Size();k++) block_particles.X(k).y=(T)100-block_particles.X(k).y;}
            if(fluids_parameters.particle_levelset_evolution->particle_levelset.positive_particles(block)){
                PARTICLE_LEVELSET_PARTICLES<TV>& block_particles=*fluids_parameters.particle_levelset_evolution->particle_levelset.positive_particles(block);
                for(int k=0;k<block_particles.Size();k++) block_particles.X(k).y=(T)100-block_particles.X(k).y;}}}

    PHYSBAM_DEBUG_WRITE_SUBSTEP("initialized phi",0,0);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
bool Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==10){
        // should write something to rasterize source mass using our negative material calc in cell
        const T_GRID& grid=*fluids_parameters.grid;
        RANGE<TV> source_one_box((T).2,(T).3,(T).14,(T).16);ORIENTED_BOX<TV> oriented_one_box(source_one_box,ROTATION<TV>::From_Angle((T)pi/4),source_one_box.Minimum_Corner());
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV source_X=world_to_source.Homogeneous_Times(iterator.Location());T source_phi=oriented_one_box.Signed_Distance(source_X);
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),source_phi);}
        RANGE<TV> source_two_box((T).4,(T).5,(T).14,(T).16);ORIENTED_BOX<TV> oriented_two_box(source_two_box,ROTATION<TV>::From_Angle((T)pi/3),source_two_box.Minimum_Corner());
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV source_X=world_to_source.Homogeneous_Times(iterator.Location());T source_phi=oriented_two_box.Signed_Distance(source_X);
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),source_phi);}}
    return false;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<T,FACE_INDEX<2> >& face_velocities,ARRAY<bool,FACE_INDEX<2> >& psi_N,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==10){
        RANGE<TV> local_box;ORIENTED_BOX<TV> oriented_box;
        TV velocity;
        T theta;

        theta=(T)pi/4;local_box=RANGE<TV>((T).2,(T).3,(T).14,(T).16);oriented_box=ORIENTED_BOX<TV>(local_box,ROTATION<TV>::From_Angle(theta),local_box.Minimum_Corner());
        velocity=TV(cos(theta),sin(theta));
        for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
            if(oriented_box.Lazy_Inside(source_X))
                face_velocities.Component(iterator.Axis())(iterator.Face_Index())=velocity(iterator.Axis());}

        theta=(T)pi/3;local_box=RANGE<TV>((T).4,(T).5,(T).14,(T).16);oriented_box=ORIENTED_BOX<TV>(local_box,ROTATION<TV>::From_Angle(theta),local_box.Minimum_Corner());
        velocity=TV(cos(theta),sin(theta));
        for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
            if(oriented_box.Lazy_Inside(source_X))
            face_velocities.Component(iterator.Axis())(iterator.Face_Index())=velocity(iterator.Axis());}}
}
//#####################################################################
// Function Get_Analytic_Velocities
//#####################################################################
void Get_Analytic_Velocities(const T time) const PHYSBAM_OVERRIDE
{
    // call this here to allow use of velocity divergence in our conservation
    // just overwrite fluids_parameters.face_velocities
    PHYSBAM_FATAL_ERROR("broken");
#if 0
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    if(test_number == 1 || test_number==11)
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(axis==1) face_velocities.Component(0)(face)=(T)0;
            else face_velocities.Component(1)(face)=(T).625;}
    else if(test_number==12)
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(axis==1) face_velocities.Component(0)(face)=(T)0;
            else face_velocities.Component(1)(face)=(T)-.625;}
    else if(test_number==2)
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
            face_velocities.Component(iterator.Axis())(iterator.Face_Index())=root_two_over_two;
    else if(test_number==3 || test_number==4){
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();TV location=iterator.Location();
            if(axis==1) face_velocities.Component(axis)(face_index)=pi_over_314*(zalesak_velocity_center.y-location.y);
            else face_velocities.Component(axis)(face_index)=pi_over_314*(location.x-zalesak_velocity_center.x);}}
    else if(test_number==5 || test_number==13){
        T time_reversal_factor=test_number==5?(T)cos((T)pi*time/period):(T)1;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT face_index=iterator.Face_Index();int axis=iterator.Axis();TV location=iterator.Location();
            T x=location.x,y=location.y;
            if(axis==1) face_velocities.Component(axis)(face_index)=-sqr((T)sin((T)pi*x))*(T)2*(T)sin((T)pi*y)*(T)cos((T)pi*y)*time_reversal_factor;
            else face_velocities.Component(axis)(face_index)=(T)2*sin((T)pi*x)*cos((T)pi*x)*sqr(sin((T)pi*y))*time_reversal_factor;}}
    else if (test_number==6){
        T time_reversal_factor=(T)cos((T)pi*time/period);
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT face_index=iterator.Face_Index();int axis=iterator.Axis();TV location=iterator.Location();
            T x=location.x,y=location.y;
            if(axis==1) face_velocities.Component(axis)(face_index)=sin((T)(4*pi*(.5+x)))*sin((T)(4*pi*(.5+y)))*time_reversal_factor;
            else face_velocities.Component(axis)(face_index)=cos((T)(4*pi*(.5+x)))*cos((T)(4*pi*(.5+y)))*time_reversal_factor;}}
#endif
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==5 || test_number==6 || test_number==13){
        Get_Analytic_Velocities(period/(T)4.0);
        dt=min(dt,(T).5*fluids_parameters.particle_levelset_evolution->CFL());
        Get_Analytic_Velocities(time);
    }
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(TV& X,TV& V,const T r,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE
{
    return true;
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==6){
        T_GRID& grid=*fluids_parameters.grid;
        TV DX=grid.Domain().Edge_Lengths(),X=particles.X(index),X_new=X+dt*V,X_new_original=X_new;
        if(particle_type==PARTICLE_LEVELSET_REMOVED_NEGATIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE){
            if(X_new.x>grid.domain.max_corner.x) X.x-=DX.x;else if(X_new.x<grid.domain.min_corner.x) X.x+=DX.x;
            if(X_new.y>grid.domain.max_corner.y) X.y-=DX.y;else if(X_new.y<grid.domain.min_corner.y) X.y+=DX.y;
            particles.X(index)=X;}
        else{
            if(X_new.x>grid.domain.max_corner.x) X_new.x-=DX.x;else if(X_new.x<grid.domain.min_corner.x) X_new.x+=DX.x;
            if(X_new.y>grid.domain.max_corner.y) X_new.y-=DX.y;else if(X_new.y<grid.domain.min_corner.y) X_new.y+=DX.y;
            if(X_new!=X_new_original){particles.X(index)=X_new;V=TV();}}}
    else SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Adjust_Particle_For_Domain_Boundaries(particles,index,V,particle_type,dt,time);
}
//#####################################################################
// Function Get_Analytic_Velocity
//#####################################################################
VECTOR<T,2> Get_Analytic_Velocity(const VECTOR<T,2>& location,const T time) const PHYSBAM_OVERRIDE
{
    if(test_number==1 || test_number==11) return TV((T)0,(T).625);
    else if (test_number==12) return TV((T)0,(T)-.625);
    else if(test_number==2) return TV::All_Ones_Vector()*root_two_over_two;
    else if(test_number==3 || test_number==4) return TV(zalesak_velocity_center.y-location.y,location.x-zalesak_velocity_center.x)*pi_over_314;
    else if(test_number==5 || test_number==13){
        T time_reversal_factor=test_number==5?(T)cos((T)pi*time/period):(T)1;T x=location.x,y=location.y;
        return TV(-sqr((T)sin((T)pi*x))*(T)2*(T)sin((T)pi*y)*(T)cos((T)pi*y),(T)2*sin((T)pi*x)*cos((T)pi*x)*sqr(sin((T)pi*y)))*time_reversal_factor;}
    else if (test_number==6){
        T time_reversal_factor=(T)cos((T)pi*time/period);T x=location.x,y=location.y;
        return TV(sin((T)(4*pi*(.5+x)))*sin((T)(4*pi*(.5+y))),cos((T)(4*pi*(.5+x)))*cos((T)(4*pi*(.5+y))))*time_reversal_factor;}
    return TV();
}
//#####################################################################
};    
}
#endif
