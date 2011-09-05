//#####################################################################
// Copyright 2002-2005, Ron Fedkiw, Duc Nguyen, Andrew Selle, Tamar Shinar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_FLAME 
//#####################################################################
#ifndef __SIMPLE_FLAME__
#define __SIMPLE_FLAME__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T,class RW=T>
class SIMPLE_FLAME:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;typedef VECTOR<T,2> TV;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::first_frame;using BASE::data_directory;
    using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::verbose_dt;

    BOX_2D<T> source_box;
    BOX_2D<T> box;
    T fuel_inject_speed;
    T source_xmin,source_xmax,source_ymin,source_ymax;
    RIGID_BODY_LIST<T,TV>& rigid_body_list;
    bool use_object;
    
    SIMPLE_FLAME(int argc,char* argv[])
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(2,fluids_parameters.FIRE),source_box(T(.35),T(.65),T(-.5),T(0.2)),box(T(.25),T(.75),T(-.5),T(0.2)),
         rigid_body_list(solids_parameters.rigid_body_parameters.list),use_object(false)
    {
        fluids_parameters.write_particles=true;
        PARAMETER_LIST parameters;
        output_directory="Simple_Flame/output";
        fluids_parameters.use_reacting_flow=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][2]=false;fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.grid->Initialize(100,200,T(-.5),T(1.5),0,4);
        last_frame=parameters.Get_Parameter("last_frame,l",int(512),"Last frame of animation",true);
        fluids_parameters.temperature_container.Set_Ambient_Temperature(T(283.15));fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);
        fluids_parameters.density_container.Set_Ambient_Density(0);
        fluids_parameters.normal_flame_speed=parameters.Get_Parameter("normal_flame_speed,fs",T(0.03));
        fuel_inject_speed=parameters.Get_Parameter("fuel_inject_speed,i",T(0.4));
        fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;
        fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=parameters.Get_Parameter("buoyancy_constant,b",T(0));
        fluids_parameters.gravity=0;
        fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=parameters.Get_Parameter("confinement_parameter,epsilon",T(.1));
        fluids_parameters.write_debug_data=fluids_parameters.write_velocity=true;
        fluids_parameters.object_friction=use_object?1:0;
        fluids_parameters.incompressible_iterations=200;
//      fluids_parameters.use_old_velocities_for_boundary_conditions=true;
        fluids_parameters.densities(1)=(T)1;
        fluids_parameters.densities(2)=(T).1;
        fluids_parameters.normal_flame_speeds(1,2)=fluids_parameters.normal_flame_speeds(2,1)=(T).03;
        fluids_parameters.normal_flame_speeds(1,1)=fluids_parameters.normal_flame_speeds(2,2)=(T)0;
        fluids_parameters.curvature_flame_speeds(1,2)=fluids_parameters.curvature_flame_speeds(2,1)=(T)0;
        fluids_parameters.curvature_flame_speeds(1,1)=fluids_parameters.curvature_flame_speeds(2,2)=(T)0;
        fluids_parameters.fuel_region(1)=true;
        fluids_parameters.fuel_region(2)=false;
    }

    virtual ~SIMPLE_FLAME()
    {}

    // unused callbacks
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    if(use_object) fluids_parameters.Use_Fluid_Coupling_Defaults();
    else fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    TV fuel_inject_velocity(0,fuel_inject_speed);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(source_box.Lazy_Inside(iterator.Location())){
        fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=fuel_inject_velocity[iterator.Axis()];}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    TV fuel_inject_velocity(0,fuel_inject_speed);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(box.Lazy_Inside(iterator.Location())){
        fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(iterator.Axis())(iterator.Face_Index())=true;
        fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=fuel_inject_velocity[iterator.Axis()];}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<ARRAY<T,VECTOR<int,2> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        if(source_box.Lazy_Inside(iterator.Location()-TV(T(0),(T)0.3))) 
            phis(1)(iterator.Cell_Index())=-grid.dx;
        else phis(1)(iterator.Cell_Index())=grid.dx;
    ARRAY<T,VECTOR<int,2> >::copy(-1,phis(1),phis(2));
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    Adjust_Phi_With_Source(source_box,1,MATRIX<T,3>::Identity_Matrix());
}
//#####################################################################
// Function Initial_Phi_Object
//#####################################################################
T Initial_Phi_Object(const TV& X) const
{
    if(use_object) return rigid_body_particles.Rigid_Body(1).Implicit_Geometry_Extended_Value(X);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    if(use_object){
        int id=rigid_body_list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies_2D/circle",(T).1,true,true,false);
        rigid_body_particles.Rigid_Body(id).frame.t=TV((T)1.25,(T).55);
        rigid_body_particles.Rigid_Body(id).is_kinematic=true;
        fluids_parameters.collision_bodies_affecting_fluid->Add_Bodies(rigid_body_particles);}
}
//#####################################################################
// Function Update_Kinematic_Rigid_Body_States
//#####################################################################
void Update_Kinematic_Rigid_Body_States(const T dt,const T time) PHYSBAM_OVERRIDE{
    if(use_object){
        INTERPOLATION_CURVE<T,TV> motion_curve;
/*      motion_curve.Add_Control_Point(0,TV((T).75,(T).5));
        motion_curve.Add_Control_Point(.2,TV((T).5,(T).5));
        motion_curve.Add_Control_Point(.4,TV((T).25,(T).5));
        motion_curve.Add_Control_Point(.6,TV((T).5,(T).5));
        motion_curve.Add_Control_Point(.4,TV((T).25,(T).5));*/
        motion_curve.Add_Control_Point(0,TV((T).9,(T).5));
        motion_curve.Add_Control_Point(1.7,TV((T).5,(T).5));
        motion_curve.Add_Control_Point(2.8,TV((T).2,(T).6));
        motion_curve.Add_Control_Point(3.6,TV((T).9,(T)1));
        motion_curve.Add_Control_Point(4.4,TV((T).5,(T)1.2));

        rigid_body_particles.Rigid_Body(1).frame.t=motion_curve.Value(time);
        rigid_body_particles.Rigid_Body(1).twist.linear=motion_curve.Derivative(time);}
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    if(use_object){
        const TV& velocity=rigid_body_particles.Rigid_Body(1).twist.linear;
        T rigid_dt_denominator=(T)fabs(velocity.x)/grid.dx+(T)fabs(velocity.y)/grid.dy;
        if(rigid_dt_denominator>1e-8) dt=min(dt,1/rigid_dt_denominator);}
}
//#####################################################################
};    
}
#endif
