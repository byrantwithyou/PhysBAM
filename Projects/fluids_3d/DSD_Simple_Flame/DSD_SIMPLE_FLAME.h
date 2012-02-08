//#####################################################################
// Copyright 2001-2005, Duc Nguyen, Ronald Fedkiw, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DSD_SIMPLE_FLAME  
//##################################################################### 
#ifndef __DSD_SIMPLE_FLAME__
#define __DSD_SIMPLE_FLAME__


#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;

template<class T_input,class RW=T_input>
class DSD_SIMPLE_FLAME:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> >,RW>
{
    typedef T_input T;
public:
    typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID<TV>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID<TV>::FAST_LEVELSET T_FAST_LEVELSET;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;typedef VECTOR<T,3> TV;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::write_output_files;
    using BASE::output_directory;using BASE::restart;using BASE::restart_frame;using BASE::verbose_dt;using BASE::data_directory;

    CYLINDER<T> source_cylinder;
    TV fuel_inject_velocity;
    RIGID_BODY_LIST<T,TV>& rigid_body_list;
    bool use_object;
    T Dn_initial;

    DSD_SIMPLE_FLAME()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>(2,fluids_parameters.FIRE),
        source_cylinder(TV((T).5,(T)-.25,(T).5),TV((T).5,(T).25,(T).5),(T).15),fuel_inject_velocity(0,5,0),rigid_body_list(solids_parameters.rigid_body_parameters.list),use_object(false)
    {
        fluids_parameters.grid->Initialize(33,65,33,T(0),T(1),0,2,T(0),T(1));
        fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[1][0]=true;
        last_frame=int(T(20)*frame_rate);
        fluids_parameters.incompressible_iterations=200;
        fluids_parameters.gravity=0;
        fluids_parameters.density_container.Set_Ambient_Density(0);
        fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);fluids_parameters.temperature_container.Set_Ambient_Temperature((T)283.15);
        fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;
        fluids_parameters.density_buoyancy_constant=fluids_parameters.temperature_buoyancy_constant=T(.001);
        fluids_parameters.use_vorticity_confinement=false;
        fluids_parameters.use_reacting_flow=true;
        fluids_parameters.use_dsd=true;
        fluids_parameters.densities(1)=(T)1;
        fluids_parameters.densities(2)=(T).1;
        fluids_parameters.normal_flame_speeds(1,2)=fluids_parameters.normal_flame_speeds(2,1)=(T).03;
        fluids_parameters.normal_flame_speeds(1,1)=fluids_parameters.normal_flame_speeds(2,2)=(T)0;
        fluids_parameters.fuel_region(1)=true;
        fluids_parameters.fuel_region(2)=false;
        fluids_parameters.delete_fluid_inside_objects=true;
        write_output_files=true;fluids_parameters.write_debug_data=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.number_particles_per_cell=0;
        output_directory="DSD_Simple_Flame/output";

        // initialize dsd and dsd parameters
        fluids_parameters.incompressible_multiphase->projection.Initialize_Dsd();
        DETONATION_SHOCK_DYNAMICS<T,GRID<TV> >& dsd=*fluids_parameters.incompressible_multiphase->projection.dsd;
        dsd.order=3;
        dsd.Dcj=0.2;
        dsd.Dcj_min_clamp=0.0;
        dsd.Dcj_max_clamp=2*dsd.Dcj;
        dsd.A_coeff=0.2;
        dsd.B_coeff=1;
        dsd.C_coeff=100;
        dsd.D_coeff=0;
        dsd.mutheta=(T)5;
        dsd.dtheta=(T)1;
        dsd.use_log_Lcj=true;
        dsd.nb_width=5;
    }

    virtual ~DSD_SIMPLE_FLAME()
    {}

    // unused callbacks
    void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE {}
    void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE {}
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
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(source_cylinder.Lazy_Inside(iterator.Location())){
        fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=fuel_inject_velocity[iterator.Axis()];}

    // TODO : move this dsd initialization elsewhere (after the levelset has been fast marched)
    DETONATION_SHOCK_DYNAMICS<T,GRID<TV> >& dsd=*fluids_parameters.incompressible_multiphase->projection.dsd;
    T_ARRAYS_SCALAR::copy(Dn_initial,fluids_parameters.incompressible_multiphase->projection.dsd->Dn.array);
    T_FAST_LEVELSET& levelset=*fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.levelsets(1);
    levelset.Compute_Curvature();dsd.curvature_old.array=*levelset.curvature;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    T cylinder_thickness=3*grid.Minimum_Edge_Length();
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(!source_cylinder.Outside(iterator.Location(),cylinder_thickness)){
        fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(iterator.Axis())(iterator.Face_Index())=true;
        fluids_parameters.incompressible->projection.face_velocities.Component(iterator.Axis())(iterator.Face_Index())=fuel_inject_velocity[iterator.Axis()];}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<ARRAY<T,VECTOR<int,3> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        if(source_cylinder.Lazy_Inside(iterator.Location())) phis(1)(iterator.Cell_Index())=-grid.dx;
        else phis(1)(iterator.Cell_Index())=grid.dx;
    ARRAY<T,VECTOR<int,3> >::copy(-1,phis(1),phis(2));
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=*fluids_parameters.grid;
    ARRAY<ARRAY<T,VECTOR<int,3> > >& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        if(source_cylinder.Lazy_Inside(iterator.Location())) phis(1)(iterator.Cell_Index())=-grid.dx;
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
        int id=rigid_body_list.template Add_Rigid_Body<float>(data_directory+"/Rigid_Bodies/sphere",(T).1,true,true,false);
        std::cout<<"Added rigid_body :"<<id<<std::endl;
        rigid_body_particles.Rigid_Body(id).frame.t=TV((T).75,(T).55,(T).75);
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
        motion_curve.Add_Control_Point(0,TV((T).75,(T).5,(T).5));
        motion_curve.Add_Control_Point(1.7,TV((T).5,(T).5,(T).5));
        motion_curve.Add_Control_Point(2.8,TV((T).2,(T).6,(T).5));
        motion_curve.Add_Control_Point(3.6,TV((T).75,(T)1,(T).5));
        motion_curve.Add_Control_Point(4.4,TV((T).5,(T)1.2,(T).5));

        rigid_body_particles.Rigid_Body(1).frame.t=motion_curve.Value(time);
        rigid_body_particles.Rigid_Body(1).twist.linear=motion_curve.Derivative(time);}
}
//#####################################################################
// Function Limit_Dt
//#####################################################################
void Limit_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    if(use_object){
        TV& velocity=rigid_body_particles.Rigid_Body(1).twist.linear;
        T rigid_dt_denominator=TV::Dot_Product(abs(velocity),fluids_parameters.grid->One_Over_DX());
        if(rigid_dt_denominator>1e-8) dt=min(dt,1/rigid_dt_denominator);}
}
//##################################################################### 
};      
#endif


