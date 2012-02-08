//#####################################################################
// Copyright 2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_CONTROL
//#####################################################################
// Provides 3 tests that are consistent across all water code
//   1. Rising drop
//   2. Sine wave
//   3. Deep water
// Also supports a variety of standard resolutions.
//#####################################################################
// TODO: Add deep water
#ifndef __FLUID_CONTROL__
#define __FLUID_CONTROL__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Fourier_Transforms_Calculations/DEEP_WATER_EVOLUTION_HEIGHTS.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/FLUID_CONTROL_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Solids_And_Fluids/FLUID_CONTROL_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Standard_Tests/WATER_STANDARD_TESTS_2D.h>
namespace PhysBAM{

template<class T,class RW>
class FLUID_CONTROL:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>,public FLUID_CONTROL_CALLBACKS<GRID<TV> >
{
public:
    typedef GRID<TV> T_GRID;
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID<TV>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;
    using BASE::parse_args;using BASE::test_number;

    int test_number;
    RIGID_BODY_LIST<T,TV>& rigid_body_list;
    GRID<TV> grid;
    SEGMENTED_CURVE_2D<T>* sine_curve;

    T time_control_start,time_control_end;
    TV circle_center_start,circle_center_end;
    T circle_radius;
    T plane_initial_y,plane_control_y;
    
    typedef PhysBAM::VECTOR<T,1> TV_HORIZONTAL;
    mutable DEEP_WATER_EVOLUTION_HEIGHTS<TV_HORIZONTAL> deep_water;
    mutable T time_old;
        
    FLUID_CONTROL(const int test_number_input,const int resolution)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(1,fluids_parameters.WATER),rigid_body_list(solids_parameters.rigid_body_parameters.list),time_old(0)
    {
        LOG::cout<<"Running Fluid Control Test "<<test_number<<" at resolution "<<resolution<<std::endl;

        // set up the standard fluid environment
        frame_rate=24;
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
        fluids_parameters.use_body_force=true;

        // set up the domain
        int cells=1*resolution;
        if(test_number==1||test_number==2||test_number==3){
            first_frame=0;last_frame=1000;
            grid.Initialize(16*cells,16*cells,0,1.5,0,1.5);
            *fluids_parameters.grid=grid;}
        else{LOG::cerr<<"unrecognized test number "<<test_number<<std::endl;exit(1);}

        // set up the sine curve
        int horizontal_cells=15*cells+1;
        if(test_number==2){
            SEGMENT_MESH* sine_segment=new SEGMENT_MESH();sine_segment->Initialize_Straight_Mesh(horizontal_cells,false);
            PARTICLES<TV>* sine_particles=new PARTICLES<TV>();sine_particles->Increase_Array_Size(horizontal_cells);
            for(int p=0;p<horizontal_cells;p++){T t=(p-1)/(T)(horizontal_cells-1);sine_particles->X(p)=TV(1.5*t,.4+.25*sin(13*pi*t));}
            sine_curve=new SEGMENTED_CURVE_2D<T>(*sine_segment,*sine_particles);}

        output_directory=STRING_UTILITIES::string_sprintf("Fluid_Control/Test_%d__Resolution_%d_%d",test_number,grid.m-1,grid.n-1);
        time_control_start=0;time_control_end=20;
        circle_center_start=TV(.75,.15);circle_center_end=TV(.75,1.1);circle_radius=.2;
        plane_initial_y=.35,plane_control_y=.25;

        //set up for deep water
        if(test_number==3){
            deep_water.grid=grid.Get_Horizontal_Grid();
            deep_water.phillips_spectrum.amplitude=7e-6;
            deep_water.period=0;
            deep_water.phillips_spectrum.wind=VECTOR<T,1>(10);
            deep_water.filter_high_frequencies=true;
            deep_water.Initialize();}
    }
    
    virtual ~FLUID_CONTROL()
    {}
    
    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    
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
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next())
        fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=(T)0;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=iterator.Location().y-plane_initial_y;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_particles);
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    BASE::Construct_Levelsets_For_Objects(time);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
void Update_Fluid_Parameters(const T dt,const T time) PHYSBAM_OVERRIDE
{
    BASE::Update_Fluid_Parameters(dt,time);
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
// Function Phi
//#####################################################################
T Phi(const TV& location,const T time) const PHYSBAM_OVERRIDE
{
    if(test_number==1){
        T clamped_time=clamp(time,time_control_start,time_control_end);
        TV circle_center=LINEAR_INTERPOLATION<T,TV>::Linear(time_control_start,time_control_end,circle_center_start,circle_center_end,clamped_time);
        CIRCLE<T> target_circle=CIRCLE<T>(circle_center,circle_radius);
        return min(target_circle.Signed_Distance(location),location.y-plane_control_y);}
    else if(test_number==2) return sine_curve->Calculate_Signed_Distance(location,0);
    else return (T)0;
}
//#####################################################################
// Function Face_Velocity
//#####################################################################
T Face_Velocity(const TV& location,const int axis,const T time) const PHYSBAM_OVERRIDE
{
    typedef ARRAY<T,VECTOR<int,1> > T_ARRAYS_HORIZONTAL;typedef VECTOR<int,2> TV_INT;
    if(test_number==1){
        if(axis==2){
            T clamped_time=clamp(time,time_control_start,time_control_end);
            TV circle_center=LINEAR_INTERPOLATION<T,TV>::Linear(time_control_start,time_control_end,circle_center_start,circle_center_end,clamped_time);
            if(CIRCLE<T>(circle_center,circle_radius).Signed_Distance(location)<location.y-plane_control_y) return (circle_center_end.y-circle_center_start.y)/(time_control_end-time_control_start);
            else return 0;}
        else return 0;}
    else if(test_number==2) return 0;
    else if(test_number==3){
        if(axis==1) return 0;
        TV_INT cell=fluids_parameters.grid->Clamp_To_Cell(location);
        T dt=time-time_old;if(dt>0) deep_water.Advance_Height(dt);time_old=time;
        deep_water.Intersect_With_Bodies(rigid_body_particles);
        T_ARRAYS_HORIZONTAL v(fluids_parameters.grid->Get_Regular_Grid().Get_Horizontal_Grid(),false,false);
        deep_water.Get_Vertical_Velocity(v);
        return v(cell.x)*5e5;}
    else return 0;
}
//#####################################################################
// Function Normal
//#####################################################################
virtual TV Normal(const TV& location,const T time)const PHYSBAM_OVERRIDE
{
    if(test_number==1){
        T clamped_time=clamp(time,time_control_start,time_control_end);
        TV circle_center=LINEAR_INTERPOLATION<T,TV>::Linear(time_control_start,time_control_end,circle_center_start,circle_center_end,clamped_time);
        if(CIRCLE<T>(circle_center,circle_radius).Signed_Distance(location)<location.y-plane_control_y)
            return CIRCLE<T>(circle_center,circle_radius).Normal(location) PHYSBAM_OVERRIDE;
        else return TV(0,1);}
    else if(test_number==2) return TV(0,1); // TODO: Fix normal for sine curve
    else return TV();
}
//#####################################################################
// Function Control_Falloff
//#####################################################################
T Control_Falloff(const TV& location,const T time) const PHYSBAM_OVERRIDE
{
    if(test_number==1){
        return 1;
        //if(location.x>.55&&location.x<.95) return 1; else return 0;
        T clamped_time=clamp(time,time_control_start,time_control_end);
        TV circle_center=LINEAR_INTERPOLATION<T,TV>::Linear(time_control_start,time_control_end,circle_center_start,circle_center_end,clamped_time);
        if(CIRCLE<T>(circle_center,circle_radius).Signed_Distance(location)<0) return 1; else return 0;}
    else if(test_number==2) return 1;
    else if(test_number==3) return 1;
    else return 0;
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(T_FACE_ARRAYS_SCALAR& force,const T dt,const T time) PHYSBAM_OVERRIDE
{
    BASE::Get_Body_Force(force,dt,time);
    PROJECTION_UNIFORM<T,T_GRID> projection(*fluids_parameters.grid);
    FLUID_CONTROL_UNIFORM<T,T_GRID> control(fluids_parameters.particle_levelset_evolution->particle_levelset.levelset,fluids_parameters.incompressible->projection.face_velocities,
        fluids_parameters.incompressible->projection.elliptic_solver->psi_N,projection,false,false,false,false,false,false);
    control.Set_Fluid_Control_Callbacks(*this);
    control.use_control_falloff=true;
    if(test_number==3){control.alpha=0;control.beta=20;}//dont care about the shape feedback
    else control.alpha=10;control.beta=20;
    control.Get_Fluid_Control_Force(force,time);
}
//#####################################################################
};
}
#endif

