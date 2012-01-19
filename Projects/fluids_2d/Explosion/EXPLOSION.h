//#####################################################################
// Copyright 2002-2004, Ron Fedkiw, Frank Losasso, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXPLOSION
//##################################################################### 
#ifndef __EXPLOSION__
#define __EXPLOSION__

#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>

namespace PhysBAM{

template <class T,class RW=T>
class EXPLOSION:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::abort_when_dt_below;

    T rho,rho_bottom,rho_top,buoyancy_constant,thermal_buoyancy_constant;T explosion_divergence;T explosion_end_time;
    BOX_2D<T> source_domain;

    EXPLOSION()
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.SMOKE)
    {
        fluids_parameters.grid.Initialize(101,151,0,1,0,1.5);
        first_frame=0;
        last_frame=3840;
        frame_rate=24;
        fluids_parameters.cfl=3;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).05;
        fluids_parameters.kolmogorov=(T)0;
        rho=1;rho_bottom=1;rho_top=(T).65;buoyancy_constant=0;fluids_parameters.gravity=0;
        fluids_parameters.temperature_container.Set_Cooling_Constant(0);fluids_parameters.temperature_products=3000;
        fluids_parameters.write_velocity=true;fluids_parameters.use_body_force=true;
        fluids_parameters.write_debug_data=true;solids_parameters.deformable_body_parameters.write=false;solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=false;
        source_domain=BOX_2D<T>((T).45,(T).55,(T)0,(T).1);
        fluids_parameters.Initialize_Domain_Boundary_Conditions();
        output_directory="Explosion/output";

        explosion_divergence=20;
        explosion_end_time=(T).5;
        fluids_parameters.use_non_zero_divergence=true;
        thermal_buoyancy_constant=.00005;
    }

    ~EXPLOSION()
    {}

//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    if(time<=explosion_end_time) for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)if(source_domain.Lazy_Inside(grid.X(i,j))){
        fluids_parameters.density_container.density_2d(i,j)=rho;fluids_parameters.temperature_container.temperature_2d(i,j)=fluids_parameters.temperature_products;}
    // keep density >= 0 and T >=0
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++){
        fluids_parameters.density_container.density_2d(i,j)=max((T)0,fluids_parameters.density_container.density_2d(i,j));
        fluids_parameters.temperature_container.temperature_2d(i,j)=max((T)fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature_2d(i,j));}
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Get_Divergence(ARRAY<T,VECTOR<int,2> >& divergence,const T dt,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    T expansion=explosion_divergence*sin(time)/exp(time);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) if(source_domain.Lazy_Inside(grid.X(i,j))) divergence(i,j)=expansion;
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >& force,const T dt,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;DENSITY_CONTAINER<T>& density=fluids_parameters.density_container;
    TEMPERATURE_CONTAINER<T>& temperature=fluids_parameters.temperature_container;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){ // y-direction forces only
        T rho_atm=rho_bottom+(rho_top-rho_bottom)*(grid.y(j)-grid.ymin)/(grid.ymax-grid.ymin);
        if(density.density_2d(i,j)>.5){
            T density_difference=density.density_2d(i,j)-rho_atm;
            T temperature_difference=temperature.temperature_2d(i,j)-fluids_parameters.temperature_container.ambient_temperature;
            if(density_difference>0||temperature_difference>0)force(i,j).y=thermal_buoyancy_constant*temperature_difference-buoyancy_constant*density_difference;}}
}
//#####################################################################
};    
}
#endif


