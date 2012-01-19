//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNCING_DROP 
//##################################################################### 
#ifndef __BOUNCING_DROP__
#define __BOUNCING_DROP__

#include "../WATER_FREE_SURFACE_3D_EXAMPLE.h"
namespace PhysBAM{

template<class T,class RW=T>
class BOUNCING_DROP:public WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>
{
public:
    BOUNCING_DROP()
    {
        first_frame=0;last_frame=300;
        frame_rate=24;
        restart=false;restart_frame=100;
        grid.Initialize(50,50,50,0,1,0,1,0,1);
        domain_walls[1][1]=true;domain_walls[1][2]=true;domain_walls[2][1]=true;domain_walls[2][2]=false;domain_walls[3][1]=true;domain_walls[3][2]=true;
        number_particles_per_cell=16;
        write_levelset=true;write_velocity=true;write_particles=true;write_removed_positive_particles=false;write_removed_negative_particles=false;
        write_debug_data=true;
        output_directory="Bouncing_Drop/output2";
        delete_fluid_inside_objects=true;
        enforce_divergence_free_extrapolation=false;

        viscosity=(T)100;
        use_strain=true;write_strain=true;
        elastic_modulus=20000;
        plasticity_alpha=0;
        //plasticity_gamma=0;
        cfl/=4;
        frame_rate*=10;
    }

    ~BOUNCING_DROP() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    //VECTOR<T,3> center((T).5,(T).7,(T).5);
    VECTOR<T,3> center((T).5,(T).5,(T).5);
    T radius=(T).2;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++)
        phi(i,j,ij)=(grid.X(i,j,ij)-center).Magnitude()-radius;
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities()
{
    VECTOR<T,3> center((T).5,(T).7,(T).5);
    T scale=20;
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
    //    V(i,j,ij)=scale*VECTOR<T,3>(grid.x(i)-center.x,0,0);
        V(i,j,ij)=VECTOR<T,3>(0,-scale,0);
}
//#####################################################################
};      
}
#endif    


