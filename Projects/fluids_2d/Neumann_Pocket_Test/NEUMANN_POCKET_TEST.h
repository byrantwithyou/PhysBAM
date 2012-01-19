//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEUMANN_POCKET_TEST
//#####################################################################
#ifndef __NEUMANN_POCKET_TEST__
#define __NEUMANN_POCKET_TEST__

#include "../SMOKE_2D_EXAMPLE.h"

namespace PhysBAM{

template <class T,class RW=T>
class NEUMANN_POCKET_TEST:public SMOKE_2D_EXAMPLE<T,RW>
{
public:
    T rho,rho_bottom,rho_top,buoyancy_constant;
    BOX_2D<T> source_domain;

    NEUMANN_POCKET_TEST()
    {
        grid.Initialize(29,29,0,1,0,1);
        Use_Smoke_2D_Defaults();

        first_frame=0;
        last_frame=3840;
        frame_rate=24;
        cfl=3;
        domain_walls[1][1]=domain_walls[1][2]=domain_walls[2][2]=false;
        domain_walls[2][1]=true;
        use_vorticity_confinement=true;confinement_parameter=(T).05;
        kolmogorov=(T)0;
        rho=1;rho_bottom=1;rho_top=(T).65;buoyancy_constant=0;gravity=1;
        temperature_container.Set_Cooling_Constant(0);temperature_products=3000;
        write_output_files=true;write_debug_data=true;
        source_domain=BOX_2D<T>((T).4,(T).6,(T)0,(T).15);
        Initialize_Domain_Boundary_Conditions();
        output_directory="Neumann_Pocket_Test/output";

        solve_neumann_regions=true;

        Use_Smoke_2D_Defaults(); // Need to call after grid is set
    }

    ~NEUMANN_POCKET_TEST()
    {}

//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time)
{
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)if(source_domain.Lazy_Inside(grid.X(i,j))){
        density_container.density_2d(i,j)=rho;temperature_container.temperature_2d(i,j)=temperature_products;}
    // keep density >= 0 and T >=0
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++){
        density_container.density_2d(i,j)=max((T)0,density_container.density_2d(i,j));
        temperature_container.temperature_2d(i,j)=max((T)temperature_container.ambient_temperature,temperature_container.temperature_2d(i,j));}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(ARRAY<bool,VECTOR<int,2> >& psi_N_u,ARRAY<bool,VECTOR<int,2> >& psi_N_v,ARRAY<T,VECTOR<int,2> >& u_fixed,ARRAY<T,VECTOR<int,2> >& v_fixed,const T time)
{
#if 0
    for(int i=0;i<u_grid.m;i++)for(int j=0;j<u_grid.n;j++)if(source_domain.Lazy_Inside(u_grid.X(i,j)))u_fixed(i,j)=0;
    for(int i=0;i<v_grid.m;i++)for(int j=0;j<v_grid.n;j++)if(source_domain.Lazy_Inside(v_grid.X(i,j))){v_fixed(i,j)=(T).5;psi_N_v(i,j)=true;}
#endif
}
//#####################################################################
// Function Add_Neumann_Pocket
//#####################################################################
void Add_Neumann_Pocket(const BOX_2D<int> &pocket,const T u_left,const T u_right,const T v_bottom,const T v_top,
                        ARRAY<bool,VECTOR<int,2> >& psi_N_u,ARRAY<bool,VECTOR<int,2> >& psi_N_v,ARRAY<T,VECTOR<int,2> >& u_fixed,ARRAY<T,VECTOR<int,2> >& v_fixed)
{
    for(int j=pocket.ymin;j<pocket.ymax;j++){
        psi_N_u(pocket.xmin,j)=true;u_fixed(pocket.xmin,j)=u_left;
        psi_N_u(pocket.xmax,j)=true;u_fixed(pocket.xmax,j)=u_right;}
    for(int i=pocket.xmin;i<pocket.xmax;i++){
        psi_N_v(i,pocket.ymin)=true;v_fixed(i,pocket.ymin)=v_bottom;
        psi_N_v(i,pocket.ymax)=true;v_fixed(i,pocket.ymax)=v_top;}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(ARRAY<bool,VECTOR<int,2> >& psi_N_u,ARRAY<bool,VECTOR<int,2> >& psi_N_v,ARRAY<T,VECTOR<int,2> >& u_fixed,ARRAY<T,VECTOR<int,2> >& v_fixed,const T dt,const T time)
{
    Add_Neumann_Pocket(BOX_2D<int>(20,25,13,18),-1,1,-1,1,psi_N_u,psi_N_v,u_fixed,v_fixed);
    Add_Neumann_Pocket(BOX_2D<int>(5,10,13,18),0,0,0,0,psi_N_u,psi_N_v,u_fixed,v_fixed);
    for(int j=0;j<u_grid.n;j++){psi_N_u(15,j)=true;u_fixed(15,j)=0;}
    //Add_Neumann_Pocket(BOX_2D<int>(2,3,2,4),-2,3,-4,2,psi_N_u,psi_N_v,u_fixed,v_fixed);
    //Add_Neumann_Pocket(BOX_2D<int>(3,4,2,4),-1,1,-1,1,psi_N_u,psi_N_v,u_fixed,v_fixed);
}
//#####################################################################
};      
}
#endif


