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
#ifndef __SPHERE_INTO_WATER__
#define __SPHERE_INTO_WATER__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_SPHERE.h>

#include "../WATER_FREE_SURFACE_3D_EXAMPLE.h"

namespace PhysBAM{

template<class T,class RW=T>
class SPHERE_INTO_WATER : public WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>
{
public:
    SPHERE_INTO_WATER() {
        initial_time=0.;
        frame_rate=3*24;last_frame=12*frame_rate;
        restart_frame=0;
        m=140;n=110;mn=90;
        xmin=(T)-.4;xmax=(T)1;ymin=(T)0;ymax=(T)1.1;zmin=(T).05;zmax=(T).95;
        domain_walls[1][1]=true;domain_walls[1][2]=true;domain_walls[2][1]=true;domain_walls[2][2]=false;domain_walls[3][1]=true;domain_walls[3][2]=true;
        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
        number_particles_per_cell=32;
        use_removed_positive_particles=false;use_removed_negative_particles=false;
        second_order_pressure=true;
        write_levelset=true;write_velocity=true;write_removed_positive_particles=false;write_removed_negative_particles=false;write_particles=true;
        matlab_directory="SPHERE_INTO_WATER/matlab"; output_directory="SPHERE_INTO_WATER/output";
    }

    ~SPHERE_INTO_WATER() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    int m=grid.m,n=grid.n,mn=grid.mn;
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++) phi(i,j,ij)=grid.y(j)-.4;
}
//#####################################################################
// Function Adjust_Phi_And_Get_Velocities_For_Objects
//#####################################################################
void Adjust_Phi_And_Get_Velocities_For_Objects(const T time)
{    
    int i,j,ij;
    // falling sphere
    SPHERE<T> sphere(VECTOR<T,3>(.85,.55,.5),.1);
    double temp_time=time;
    if(time >= .075) 
        temp_time=.075;
    for(i=1;i<=m;i++) for(j=1;j<=n;j++) for(ij=1;ij<=mn;ij++)
    {
        VECTOR_3D point_loc(grid->x(i),grid->y(j),grid->z(ij));
        if (sphere.Inside(point_loc))
        {
            //always set phi
            VECTOR_3D surface=sphere.Surface(point_loc);
            VECTOR_3D normal=sphere.Normal(surface);
            surface+=epsilon*normal;
            (*phi)(i,j,ij)=interpolate(*grid,*phi,surface.x,surface.y,surface.z);

            if(time <= .075)
            {
                (*psi_N)(i,j,ij)=1;
                (*u_fixed)(i,j,ij)=-6;
                (*v_fixed)(i,j,ij)=-6;
                (*w_fixed)(i,j,ij)=0;
            }
            else
            {
                (*psi_N)(i,j,ij)=1;
                (*u_fixed)(i,j,ij)=0;
                (*v_fixed)(i,j,ij)=0;
                (*w_fixed)(i,j,ij)=0;
            }
        }
    }
}
//#####################################################################
};      
}
#endif    
