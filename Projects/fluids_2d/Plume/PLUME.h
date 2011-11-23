//#####################################################################
// Copyright 2002-2004, Ron Fedkiw, Frank Losasso, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLUME
//##################################################################### 
#ifndef __PLUME__
#define __PLUME__

#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>

namespace PhysBAM{


template <class T,class RW=T>
class PLUME:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::last_frame;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::abort_when_dt_below;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::write_time;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::write_frame_title;

    T rho,rho_bottom,rho_top,buoyancy_constant;
    BOX_2D<T> source_domain;
    MATRIX<T,3> world_to_source;

    PLUME()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.SMOKE)
    {
        fluids_parameters.grid.Initialize(101,151,0,1,0,1.5);
        first_frame=0;
        last_frame=3840;
        frame_rate=24;
        fluids_parameters.cfl=.9;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).04;
        fluids_parameters.kolmogorov=(T)0;
        rho=1;rho_bottom=1;rho_top=(T).65;buoyancy_constant=0;fluids_parameters.gravity=0;
        fluids_parameters.temperature_container.Set_Cooling_Constant(0);fluids_parameters.temperature_products=3000;
        fluids_parameters.write_velocity=true;write_frame_title=true;
        fluids_parameters.write_debug_data=true;//solids_parameters.write=false;solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=false;
        source_domain=BOX_2D<T>((T).45,(T).55,(T)0,(T).1);
        world_to_source=MATRIX<T,3>::Identity_Matrix();
        fluids_parameters.Initialize_Domain_Boundary_Conditions();
        output_directory="Plume/output";
        //fluids_parameters.use_back_and_forth_advection=false;
        fluids_parameters.use_vorticity_confinement=true;
        fluids_parameters.use_maccormack_semi_lagrangian_advection=true;
    }

    virtual ~PLUME()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Adjust_Velocity_With_Objects(const T time){}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Phi(const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initial_Velocity
//#####################################################################
virtual VECTOR<T,2>  Initial_Velocity(const VECTOR<T,2>& X) const
{
    return VECTOR<T,2>();
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Adjust_Density_And_Temperature_With_Sources(source_domain,world_to_source,rho,fluids_parameters.temperature_products);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    // not equivalent to commented out code below since code below doesn't set psi_N_u's
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Get_Source_Velocities(source_domain,world_to_source,VECTOR<T,2>(0,(T).5));
#if 0
    GRID<TV> &u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid;
    PROJECTION_2D<T>& projection=fluids_parameters.incompressible.projection;LAPLACE_2D<T>& elliptic_solver=*projection.elliptic_solver;
    for(int i=1;i<=u_grid.m;i++)for(int j=1;j<=u_grid.n;j++)if(source_domain.Lazy_Inside(u_grid.X(i,j))) projection.u(i,j)=0;
    for(int i=1;i<=v_grid.m;i++)for(int j=1;j<=v_grid.n;j++)if(source_domain.Lazy_Inside(v_grid.X(i,j))){projection.v(i,j)=(T).5;elliptic_solver.psi_N_v(i,j)=true;}
#endif
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& force,const T dt,const T time) PHYSBAM_OVERRIDE
{
    return;
#if 0
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++){ // y-direction forces only
        T rho_atm=rho_bottom+(rho_top-rho_bottom)*(grid.y(j)-grid.ymin)/(grid.ymax-grid.ymin);
        if(density(i,j)>.5){T difference=density(i,j)-rho_atm;if(difference>0)force(i,j).y=-buoyancy_constant*difference;}}
#endif
}
//#####################################################################
};    
}
#endif


