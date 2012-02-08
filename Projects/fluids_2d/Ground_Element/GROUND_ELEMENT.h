//#####################################################################
// Copyright 2002, 2003, 2004, Ron Fedkiw, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GROUND_ELEMENT
//#####################################################################
#ifndef __GROUND_ELEMENT__
#define __GROUND_ELEMENT__


#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_SOLID_WALL_SLIP.h>
#include "../SMOKE_2D_EXAMPLE.h"

namespace PhysBAM{

template <class T>
class GROUND_ELEMENT:public SMOKE_2D_EXAMPLE<T>
{
public:
    T rho,rho_bottom,rho_top,buoyancy_constant,T_inject;
    BOUNDARY_SOLID_WALL_SLIP<T,VECTOR_2D<T> > slip_boundary;
    BOX_2D<T> source_domain;

    GROUND_ELEMENT()
    {
        start_frame=0;
        end_frame=384;
        frame_rate=24;
        m=500;n=100;
        xmin=(T)-5;xmax=(T)5;ymin=(T)0;ymax=(T)2;
        domain_walls[0][0]=false;domain_walls[0][1]=false;domain_walls[1][0]=true;domain_walls[1][1]=false;
        slip_boundary=BOUNDARY_SOLID_WALL_SLIP<T,VECTOR_2D<T> >(!domain_walls[0][0],!domain_walls[0][1],!domain_walls[1][0],!domain_walls[1][1]);
        use_vorticity_confinement=true;confinement_parameter=(T).05;
        komolgorov=(T)0;
        incompressible_enforce_compatibility=false;
        rho=(T)1;rho_bottom=(T)1;rho_top=(T).65;buoyancy_constant=(T)10;gravity=(T)0;
        cooling_constant=(T)0;T_air=(T)298;T_burnt=(T)3000;T_inject=(T)298;
        boundary_V=&slip_boundary;
        boundary_density=&constant_extrapolation;
        boundary_T=&constant_extrapolation;
        write_matlab_files=true;
        write_output_files=true;
        matlab_directory="Ground_Element/matlab";
        //output_directory="Ground_Element/output";
        source_domain=BOX_2D<T>((T)4,(T)5,(T)0,(T).3);
    }

    ~GROUND_ELEMENT()
    {}

//#####################################################################
// Function Update_Source
//#####################################################################
void Update_Sources(const GRID<TV>& grid,ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >& V,ARRAY<T,VECTOR<int,2> >& density,ARRAY<T,VECTOR<int,2> >& temperature,const T time)
{
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)if(source_domain.Lazy_Inside(grid.X(i,j))){
        density(i,j)=rho;temperature(i,j)=T_burnt;V(i,j).x=T(-1.3);}
    // keep density >= 0 and T >=0
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++){if(density(i,j)<0)density(i,j)=0;if(temperature(i,j)<T_air) temperature(i,j)=T_air;}
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
void Set_Boundary_Conditions(const GRID<TV>& p_grid,const GRID<TV>& u_grid,const GRID<TV>& v_grid,ARRAY<bool,VECTOR<int,2> >& psi_N_u,ARRAY<bool,VECTOR<int,2> >& psi_N_v,ARRAY<T,VECTOR<int,2> >& u_fixed,
                             ARRAY<T,VECTOR<int,2> >& v_fixed,ARRAY<bool,VECTOR<int,2> >& psi_D,ARRAY<T,VECTOR<int,2> >& p,const T time)
{
    // Set up wall boundary conditions
    SMOKE_2D_EXAMPLE<T>::Set_Boundary_Conditions(p_grid,u_grid,v_grid,psi_N_u,psi_N_v,u_fixed,v_fixed,psi_D,p,time);

    for(int i=0;i<u_grid.m;i++)for(int j=0;j<u_grid.n;j++)if(source_domain.Lazy_Inside(u_grid.X(i,j))){
        psi_N_u(i,j)=true;u_fixed(i,j)=(T)-1.3;}
}
//#####################################################################
// Function Update_Forces
//#####################################################################
void Update_Forces(const GRID<TV>& grid, const ARRAY<VECTOR_3D<T> ,VECTOR<int,2> >& V,const ARRAY<T,VECTOR<int,2> >& density_ghost,const ARRAY<T,VECTOR<int,2> >& temperature_ghost,ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >& force,
                   const T time)
{
    if(time >= 0)
        for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++){ // y-direction forces only
            T rho_atm=rho_bottom+(rho_top-rho_bottom)*(grid.y(j)-ymin)/(ymax-ymin);
            if(density_ghost(i,j)>.5){T difference=density_ghost(i,j)-rho_atm;if(difference>0)force(i,j).y=-buoyancy_constant*difference;}}
}    
//#####################################################################
};      
}
#endif


