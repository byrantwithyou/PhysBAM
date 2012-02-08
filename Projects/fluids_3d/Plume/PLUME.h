//#####################################################################
// Copyright 2001-2004, Ron Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLUME
//#####################################################################
#ifndef __PLUME__
#define __PLUME__

#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW=T>
class PLUME:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;

    T rho,rho_bottom,rho_top,buoyancy_constant;
    BOX_3D<T> source_domain;
    
    PLUME()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.SMOKE)
    {
        fluids_parameters.grid.Initialize(100,100,100,0,1,0,1,0,1);
        first_frame=0;
        last_frame=200;
        frame_rate=24;
        fluids_parameters.cfl=3;
        fluids_parameters.domain_walls[0][0]=fluids_parameters.domain_walls[0][1]=fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[2][0]=fluids_parameters.domain_walls[2][1]=false;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).3;
        fluids_parameters.kolmogorov=(T)0;
        rho=(T)1;rho_bottom=(T)1;rho_top=(T).65;buoyancy_constant=(T)0;fluids_parameters.gravity=(T)0;
        fluids_parameters.temperature_container.Set_Cooling_Constant((T)1);fluids_parameters.temperature_products=(T)3000;
        write_output_files=true;fluids_parameters.write_debug_data=true;
        source_domain=BOX_3D<T>((T).45,(T).55,(T)0,(T).1,(T).45,(T).55);
        fluids_parameters.Initialize_Domain_Boundary_Conditions();
        output_directory="Plume/output";
    }
    
    ~PLUME()
    {}
    
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<fluids_parameters.grid.m;i++)for(int j=0;j<fluids_parameters.grid.n;j++)for(int ij=0;ij<fluids_parameters.grid.mn;ij++)if(source_domain.Lazy_Inside(fluids_parameters.grid.X(i,j,ij))){
        fluids_parameters.density_container.density_3d(i,j,ij)=rho;fluids_parameters.temperature_container.temperature_3d(i,j,ij)=fluids_parameters.temperature_products;}
    // keep density >= 0 and T >=0
    for(int i=0;i<fluids_parameters.grid.m;i++)for(int j=0;j<fluids_parameters.grid.n;j++)for(int ij=0;ij<fluids_parameters.grid.mn;ij++){
        fluids_parameters.density_container.density_3d(i,j,ij)=max((T)0,fluids_parameters.density_container.density_3d(i,j,ij));
        fluids_parameters.temperature_container.temperature_3d(i,j,ij)=max((T)fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature_3d(i,j,ij));}

}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<fluids_parameters.u_grid.m;i++)for(int j=0;j<fluids_parameters.u_grid.n;j++)for(int ij=0;ij<fluids_parameters.u_grid.mn;ij++)
        if(source_domain.Lazy_Inside(fluids_parameters.u_grid.X(i,j,ij))) fluids_parameters.incompressible.projection.u(i,j,ij)=0;
    for(int i=0;i<fluids_parameters.v_grid.m;i++)for(int j=0;j<fluids_parameters.v_grid.n;j++)for(int ij=0;ij<fluids_parameters.v_grid.n;ij++)
        if(source_domain.Lazy_Inside(fluids_parameters.v_grid.X(i,j,ij))){fluids_parameters.incompressible.projection.v(i,j,ij)=(T).5;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v(i,j,ij)=true;}
    for(int i=0;i<fluids_parameters.w_grid.m;i++)for(int j=0;j<fluids_parameters.w_grid.n;j++)for(int ij=0;ij<fluids_parameters.w_grid.mn;ij++)
        if(source_domain.Lazy_Inside(fluids_parameters.w_grid.X(i,j,ij))) fluids_parameters.incompressible.projection.w(i,j,ij)=0;
}
//#####################################################################
// Function Update_Forces
//#####################################################################
void Get_Body_Force(ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& force,const T dt,const T time) PHYSBAM_OVERRIDE
{
    // control buoyancy force
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++) for(int ij=0;ij<fluids_parameters.grid.mn;ij++){ // y-direction forces only
        T rho_atm=rho_bottom+(rho_top-rho_bottom)*(fluids_parameters.grid.y(j)-fluids_parameters.grid.ymin)/(fluids_parameters.grid.ymax-fluids_parameters.grid.ymin);
        if(fluids_parameters.grid.y(j)<=2)rho_atm=1;
        T den=fluids_parameters.density_container.density_3d(i,j,ij);
        if(den>.5){ // only add the buoyancy force inside the smoke
            T difference=den-rho_atm;//if(difference > 0) difference*=3;
            if(difference>0)force(i,j,ij).y=fluids_parameters.buoyancy_constant*difference;}}
}
//#####################################################################
};      
}
#endif


