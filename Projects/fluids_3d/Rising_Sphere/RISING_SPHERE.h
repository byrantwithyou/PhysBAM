//#####################################################################
// Copyright 2002, Doug Enright, Sergey Koltakov, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RISING_SPHERE 
//##################################################################### 
//
//#####################################################################
// Enright - July 1, 2002
// Nguyen  - July 15, 2002
// Fedkiw - September 11, 2002
// Koltakov - July 23, 2003
//#####################################################################
#ifndef __RISING_SPHERE__
#define __RISING_SPHERE__

#include <Tools/Vectors/VECTOR_3D.h>
#include <Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_SPHERE.h>

#include "../WATER_FREE_SURFACE_3D_EXAMPLE.h"

namespace PhysBAM{

template<class T,class RW=T>
class RISING_SPHERE : public WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>
{
public:
    RISING_SPHERE() {
        initial_time=0.;
        frame_rate=3*24;last_frame=12*frame_rate;
        restart_frame=0;
        m=30;n=60;mn=30;//m=84;n=60;mn=84;//105,75,105
        xmin=(T)0;xmax=(T)1;ymin=(T)0;ymax=(T)2;zmin=(T)0;zmax=(T)1;
        domain_walls[0][0]=true;domain_walls[0][1]=true;domain_walls[1][0]=true;domain_walls[1][1]=false;domain_walls[2][0]=true;domain_walls[2][1]=true;
        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
        number_particles_per_cell=32;
        use_removed_positive_particles=false;use_removed_negative_particles=false;
        second_order_pressure=true;
        write_levelset=true;write_velocity=true;write_removed_positive_particles=true;write_removed_negative_particles=true;write_particles=true;
        matlab_directory="RISING_SPHERE/matlab"; output_directory="RISING_SPHERE/output";
    }

    ~RISING_SPHERE() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    int m=grid.m,n=grid.n,mn=grid.mn;
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++) phi(i,j,ij)=grid.y(j)-.5;
}
//#####################################################################
// Function Adjust_Phi_And_Get_Velocities_For_Objects
//#####################################################################
void Adjust_Phi_And_Get_Velocities_For_Objects(const T time)
{    
    int i,j,ij;
    for(i=0;i<m;i++) for(j=0;j<n;j++) for(ij=0;ij<mn;ij++){
        T y_velocity=(T)1.0;
        RENDERING_SPHERE sphere(VECTOR<T,3>(.5,.25+y_velocity*time,.5),.075);
        VECTOR<T,3> point_loc(grid.x(i),grid.y(j),grid.z(ij));
        // adjust the values of flow variables on the grid due to the presence of objects if any exist
        if (sphere.Inside(point_loc)){
            phi(i,j,ij)=-grid.dy;
            psi_N(i,j,ij)=true;
            u_fixed(i,j,ij)=0;v_fixed(i,j,ij)=y_velocity;w_fixed(i,j,ij)=0;}}
    
    T lambda=(T).25,y_velocity=(T)1,x0=(T).5,y0=(T).25+y_velocity*time,z0=(T).5,radius=(T).075;
    for(i=0;i<m+1;i++) for(j=0;j<n;j++) for(ij=0;ij<mn;ij++)
        if(radius > sqrt(sqr(grid->x(i)-x0)+sqr(grid->y(j)-y0)+sqr(grid->z(ij)-z0))) (*u)(i,j,ij)=(*u)(i,j,ij)*lambda+(1-lambda)*0;
    for(i=0;i<m+1;i++) for(j=0;j<n+1;j++) for(ij=0;ij<mn;ij++)
        if(radius > sqrt(sqr(grid->x(i)-x0)+sqr(grid->y(j)-y0)+sqr(grid->z(ij)-z0))) (*v)(i,j,ij)=(*v)(i,j,ij)*lambda+(1-lambda)*y_velocity;
    for(i=0;i<m+1;i++) for(j=0;j<n;j++) for(ij=0;ij<mn+1;ij++)
        if(radius > sqrt(sqr(grid->x(i)-x0)+sqr(grid->y(j)-y0)+sqr(grid->z(ij)-z0))) (*w)(i,j,ij)=(*w)(i,j,ij)*lambda+(1-lambda)*0;
}
//#####################################################################
};      
}
#endif    
