//#####################################################################
// Copyright 2002-2004, Ron Fedkiw, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLANE_JUMP 
//#####################################################################
#ifndef __PLANE_JUMP__
#define __PLANE_JUMP__
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>

namespace PhysBAM{

template<class T,class RW=T>
class PLANE_JUMP:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::thin_shells_semi_lagrangian_density;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::thin_shells_semi_lagrangian_temperature;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::thin_shells_semi_lagrangian_velocity;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;

    PLANE_JUMP()
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.FIRE)
    {
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.phi_boundary->Set_Constant_Extrapolation(!fluids_parameters.domain_walls[0][0],!fluids_parameters.domain_walls[0][1],!fluids_parameters.domain_walls[1][1],!fluids_parameters.domain_walls[1][0],!fluids_parameters.domain_walls[2][0],!fluids_parameters.domain_walls[2][1]);
        fluids_parameters.phi_boundary->Set_Fixed_Boundary(false);
        fluids_parameters.fluid_boundary->Set_Constant_Extrapolation(!fluids_parameters.domain_walls[0][0],!fluids_parameters.domain_walls[0][1],!fluids_parameters.domain_walls[1][1],!fluids_parameters.domain_walls[1][0],!fluids_parameters.domain_walls[2][0],!fluids_parameters.domain_walls[2][1]);
        frame_rate=240;last_frame=int(T(200)*frame_rate);
        fluids_parameters.grid.Initialize(100,100,0,1,0,1);
        fluids_parameters.temperature_container.Set_Ambient_Temperature(T(283.15));fluids_parameters.temperature_container.Set_Cooling_Constant((T)4000);
        fluids_parameters.density_container.Set_Ambient_Density(0);
        fluids_parameters.normal_flame_speed=T(1);fluids_parameters.temperature_products=3000;fluids_parameters.temperature_fuel=298;
        fluids_parameters.density_fuel=1;fluids_parameters.density=T(.2);fluids_parameters.buoyancy_constant=T(0);fluids_parameters.use_body_force=false;
        fluids_parameters.use_vorticity_confinement_fuel=fluids_parameters.use_vorticity_confinement=false;fluids_parameters.gravity=0;
        fluids_parameters.confinement_parameter_fuel=fluids_parameters.confinement_parameter=(T)0;
        fluids_parameters.write_debug_data=fluids_parameters.write_velocity=true;
        output_directory="Plane_Jump/output";
    }

    virtual ~PLANE_JUMP()
    {}

    bool Is_Fuel(const VECTOR_2D<T>& X)
    {if(X.x>=(T).25&&X.x<=(T).75)return true;else return false;}

    virtual void Initialize_Velocities()
    {for(int i=0;i<fluids_parameters.grid.m;i++)for(int j=0;j<fluids_parameters.grid.n;j++){
        if(fluids_parameters.grid.x(i)<(T).25)fluids_parameters.incompressible.V(i,j)=VECTOR_2D<T>(-4,0);
        else if(fluids_parameters.grid.x(i)>(T).75)fluids_parameters.incompressible.V(i,j)=VECTOR_2D<T>(4,0);
        else fluids_parameters.incompressible.V(i,j)=VECTOR_2D<T>(0,0);}}
    
    void Get_Source_Velocities(const T time)
    {}

    void Initialize_Phi()
    {for(int i=0;i<fluids_parameters.grid.m;i++)for(int j=0;j<fluids_parameters.grid.n;j++)
        if(Is_Fuel(fluids_parameters.grid.X(i,j))) fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j)=-fluids_parameters.grid.dx;
        else fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi(i,j)=fluids_parameters.grid.dx;}

    void Adjust_Phi_With_Sources(const T time)
    {}
//#####################################################################
};    
}
#endif

