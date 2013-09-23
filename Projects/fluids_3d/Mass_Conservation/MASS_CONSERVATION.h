//#####################################################################
// Copyright 2006-2007, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
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
//#####################################################################
#ifndef __MASS_CONSERVATION__
#define __MASS_CONSERVATION__

#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/SEGMENT_BENDING_ELEMENTS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input,class RW=T_input>
class MASS_CONSERVATION:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;typedef VECTOR<int,3> TV_INT;
public:
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::fluids_parameters;
    using BASE::fluid_collection;using BASE::solids_parameters;using BASE::write_time;using BASE::write_frame_title;using BASE::data_directory;using BASE::abort_when_dt_below;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
    using BASE::test_number;using BASE::resolution;

    SPHERE<TV> source,other_source;
    RANGE<TV> box;
    MATRIX<T,4> world_to_source;
    TV source_velocity;
    ARRAY<T,FACE_INDEX<TV::m> > velocity;
    T initial_water_level;
    T root_three_over_three;
    T pi_over_314;
    TV zalesak_center,zalesak_velocity_center;
    RANDOM_NUMBERS<T> random;

    PARAMETER_LIST parameter_list;
    T period;

    MASS_CONSERVATION(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV >(stream_type,1,fluids_parameters.WATER)
    {
    }

    virtual ~MASS_CONSERVATION()
    {}

    // Unused callbacks
    void Initialize_Bodies() PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Update_Rigid_Bodies(const T time) {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}
    void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Velocities() PHYSBAM_OVERRIDE {}
    void Initialize_Euler_State() PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}

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
    std::string filename=STRING_UTILITIES::string_sprintf("Mass_Conservation/example_%d.param",test_number);
    if(FILE_UTILITIES::File_Exists(filename)){std::cout<<"Reading parameter file '"<<filename<<"'"<<std::endl;parameter_list.Read(filename);}

    random.Set_Seed(3891);
    root_three_over_three=(T)root_three/(T)3;
    pi_over_314 = (T)pi/314;

    // Common parameters
    output_directory=STRING_UTILITIES::string_sprintf("Mass_Conservation/example_%d_resolution_%d",test_number,resolution);
    first_frame=0;restart=false;restart_frame=0;frame_rate=300;

    // Fluids parameters
    fluids_parameters.number_particles_per_cell=16;
    fluids_parameters.reseeding_frame_rate=last_frame/5;
    fluids_parameters.incompressible_iterations=50;
    fluids_parameters.cfl=(T).5;
    fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][0]=fluids_parameters.domain_walls[1][1]=true;
    fluids_parameters.use_removed_positive_particles=true;fluids_parameters.use_removed_negative_particles=true;
    fluids_parameters.bias_towards_negative_particles=false;
    fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;fluids_parameters.variable_viscosity=false;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=true;
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
    world_to_source=MATRIX<T,4>::Identity_Matrix();

    if(test_number==1 || test_number==2){
        fluids_parameters.analytic_test=true;
        initial_water_level=-1; // no initial height
        box=RANGE<TV>((T).4,(T).6,(T).4,(T).6,(T).4,(T).6);
        fluids_parameters.grid->Initialize(TV_INT()+10*resolution+1,RANGE<TV>::Unit_Box());}
    else if(test_number==3){
        final_time=628;frame_rate=(T).5;
        zalesak_center=TV(50,75,50);
        zalesak_velocity_center=TV(50,50,50);
        fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,10*resolution+1,10*resolution+1),RANGE<TV>(TV(0,0,0),TV(100,100,100)));}
    else if(test_number==4){
        final_time=628;frame_rate=(T).5;
        zalesak_center=TV(50,50,50);
        zalesak_velocity_center=TV(50,50,50);
        fluids_parameters.grid->Initialize(TV_INT(10*resolution+1,10*resolution+1,10*resolution+1),RANGE<TV>(TV(0,0,0),TV(100,100,100)));}
    else if(test_number==5){
        fluids_parameters.particle_levelset_evolution->Use_Frozen_Velocity(false);
        fluids_parameters.grid->Initialize((1<<resolution)+1,(1<<resolution)+1,(1<<resolution)+1,RANGE<TV>::Unit_Box());
        initial_water_level=-1;
        fluids_parameters.cfl=(T).5;
        last_frame=100;
        frame_rate=(T)last_frame/final_time;
        source=SPHERE<TV>(TV((T).5,(T).75,(T).5),(T).15);}
    else if(test_number==6){
        period=final_time=3;
        fluids_parameters.grid->Initialize((1<<resolution)+1,(1<<resolution)+1,(1<<resolution)+1,RANGE<TV>::Unit_Box());
        initial_water_level=-1;
        fluids_parameters.cfl=(T).5;
        last_frame=100;
        frame_rate=(T)last_frame/final_time;
        source=SPHERE<TV>(TV((T).35,(T).35,(T).35),(T).15);}
    else if(test_number==7){
        frame_rate=500;
        fluids_parameters.grid->Initialize(TV_INT(10*resolution,(int)(7.5*resolution),10*resolution),RANGE<TV>(TV((T)0,(T)0,(T)0),TV((T)2,(T)1.5,(T)2)));
        source_velocity=TV(0,-40,0);
        initial_water_level=(T).02;
        source=SPHERE<TV>(TV(1,initial_water_level+(T).15,1),(T).05);
        fluids_parameters.surface_tension=(T)5e-6;}
    else if(test_number==8){
        frame_rate=100;
        fluids_parameters.grid->Initialize(TV_INT()+10*resolution+1,RANGE<TV>::Unit_Box());
        fluids_parameters.gravity=0;
        source=SPHERE<TV>(TV((T).1,(T).5,(T).5),(T).08);
        other_source=SPHERE<TV>(TV((T).9,(T).5,(T).5),(T).08);
        source_velocity=TV((T)2,(T)0,(T)0);
        initial_water_level=-1;}
    else if(test_number==9){
        fluids_parameters.grid->Initialize(TV_INT()+10*resolution+1,RANGE<TV>::Unit_Box());
        fluids_parameters.gravity=0;
        fluids_parameters.surface_tension=(T)1e-5;}

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
    if(test_number==1 || test_number==2 || test_number==3 || test_number==4 || test_number==5 || test_number==6){
        fluids_parameters.particle_levelset_evolution->levelset_advection.Use_Local_WENO_For_Advection();
        fluids_parameters.particle_levelset_evolution->Set_Runge_Kutta_Order_Levelset(3);
        fluids_parameters.particle_levelset_evolution->Set_Runge_Kutta_Order_Particles(3);
        fluids_parameters.analytic_test=true;
        fluids_parameters.reseeding_frame_rate=200000;
        fluids_parameters.reinitialize_geometry_frame_rate=200000;
        fluids_parameters.particle_levelset_evolution->Use_Reinitialization();
        fluids_parameters.particle_levelset_evolution->Particle_Levelset(0).reincorporate_removed_particles_everywhere=true;}
    if(test_number==6){
        fluids_parameters.phi_boundary=new BOUNDARY_MAC_GRID_PERIODIC<TV,T>();
        fluids_parameters.incompressible->Set_Custom_Boundary(*fluids_parameters.phi_boundary);
        fluids_parameters.particle_levelset_evolution->Use_Frozen_Velocity(false);}
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV >::Update_Fluid_Parameters(dt,time);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    // Not so good to set up a heaviside function here because then the interface will be exactly between the two nodes which can lead to roundoff issues when setting dirichlet cells, etc.
    GRID<TV>& grid=*fluids_parameters.grid;

    if(test_number==1 || test_number==2){
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=box.Signed_Distance(location);}}
    else if(test_number==3 || test_number==4){
        T radius=15,slot_width=5,slot_depth=(T)12.5;
        RANGE<TV> rect(zalesak_center.x-slot_width/2,zalesak_center.x+slot_width/2,0,zalesak_center.y-radius+slot_depth,0,zalesak_center.y-radius+slot_depth);
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=max((location-zalesak_center).Magnitude()-radius,-rect.Signed_Distance(location));}}
    else if(test_number==5 || test_number==6 || test_number==7 || test_number==8)
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(location.y-initial_water_level,source.Signed_Distance(location));}

    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    if(test_number==7)
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(iterator.Location().y>initial_water_level+grid.dX.Min()*2) face_velocities(axis,face)=source_velocity[axis];}
    else if(test_number==8){
        for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),other_source.Signed_Distance(location));}
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(iterator.Location().x<.5) face_velocities(axis,face)=source_velocity[axis];
            else face_velocities(axis,face)=-source_velocity[axis];}}
    else if(test_number==9){
        // list of source circles and velocities with random locations, random radii within a range, random velocities within a range.  Overlap is fine.
        int drops=40;
        T radius_min=(T).02,radius_max=(T).03;
        ARRAY<SPHERE<TV> > sources(drops);
        ARRAY<TV> velocities(drops);
        RANGE<TV> velocity_box(-(T).2,(T).2,-(T).2,(T).2,-(T).2,(T).2);
        for(int i=0;i<drops;i++){
            sources(i).center=random.Get_Uniform_Vector(grid.Domain());
            sources(i).radius=random.Get_Uniform_Number(radius_min,radius_max);
            velocities(i)=random.Get_Uniform_Vector(velocity_box);}
         for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();TV_INT cell=iterator.Cell_Index();
            T& phi=fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index());
            phi=FLT_MAX;
            for(int i=0;i<drops;i++) phi=min(phi,sources(i).Signed_Distance(location));}
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV location=iterator.Location();TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
            T phi=FLT_MAX;
            for(int i=0;i<drops;i++){
                T phi_source=sources(i).Signed_Distance(location);
                if(phi_source<phi){
                    phi=phi_source;
                    face_velocities(axis,face)=velocities(i)(axis);}}}}

    fluids_parameters.particle_levelset_evolution->Reseed_Particles(0);
    fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.refine_fmm_initialization_with_iterative_solver=true;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("initialized phi",0,0);
}
//#####################################################################
// Function Get_Analytic_Velocities
//#####################################################################
void Get_Analytic_Velocities(const T time) const PHYSBAM_OVERRIDE
{
    PHYSBAM_FATAL_ERROR("broken");
#if 0
    // call this here to allow use of velocity divergence in our conservation
    // just overwrite fluids_parameters.face_velocities
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities=fluid_collection.incompressible_fluid_collection.face_velocities;
    if(test_number==1)
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(axis==1) face_velocities.Component(0)(face)=(T).5;
            else face_velocities.Component(1)(face)=(T)0;}
    else if(test_number==2)
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next())
            face_velocities.Component(iterator.Axis())(iterator.Face_Index())=root_three_over_three;
    else if(test_number==3 || test_number==4){
        static bool initialized=false;
        static ARRAY<T,FACE_INDEX<TV::m> > stored_velocities;
        if(!initialized){
            stored_velocities.Resize(grid);
            stored_velocities.Component(2).Fill((T)0);
            for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();TV location=iterator.Location();
                if(axis==1) stored_velocities.Component(axis)(face_index)=pi_over_314*(zalesak_velocity_center.y-location.y);
                else if(axis==2) stored_velocities.Component(axis)(face_index)=pi_over_314*(location.x-zalesak_velocity_center.x);}
            initialized=true;
            face_velocities=stored_velocities;}}
    else if (test_number==6){
        T time_reversal_factor=(T)cos((T)pi*time/period);
        for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT face_index=iterator.Face_Index();int axis=iterator.Axis();TV location=iterator.Location();
            T x=location.x,y=location.y,z=location.z;
            if(axis==1) face_velocities.Component(axis)(face_index)=(T)2*sqr(sin((T)pi*x))*sin((T)two_pi*y)*sin((T)two_pi*z)*time_reversal_factor;
            else if(axis==2) face_velocities.Component(axis)(face_index)=-sin((T)two_pi*x)*sqr(sin((T)pi*y))*sin((T)two_pi*z)*time_reversal_factor;
            else face_velocities.Component(axis)(face_index)=-sin((T)two_pi*x)*sin((T)two_pi*y)*sqr(sin((T)pi*z))*time_reversal_factor;}}
#endif
}
//#####################################################################
// Function Get_Analytic_Velocity
//#####################################################################
VECTOR<T,3> Get_Analytic_Velocity(const VECTOR<T,3>& location,const T time) const PHYSBAM_OVERRIDE
{
    if(test_number==1) return TV((T).5,(T)0,(T)0);
    else if (test_number==2) return TV(root_three_over_three,root_three_over_three,root_three_over_three);
    else if(test_number==3 || test_number==4) return TV(zalesak_velocity_center.y-location.y,location.x-zalesak_velocity_center.x,(T)0)*pi_over_314;
    else if (test_number==6){
        T time_reversal_factor=(T)cos((T)pi*time/period);T x=location.x,y=location.y,z=location.z;
        return TV((T)2*sqr(sin((T)pi*x))*sin((T)two_pi*y)*sin((T)two_pi*z),-sin((T)two_pi*x)*sqr(sin((T)pi*y))*sin((T)two_pi*z),-sin((T)two_pi*x)*sin((T)two_pi*y)*sqr(sin((T)pi*z)))*time_reversal_factor;}
    return TV();
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(TV& X,TV& V,const T r, const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time) PHYSBAM_OVERRIDE
{
    return true;
}
//#####################################################################
};
}
#endif
