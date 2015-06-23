//#####################################################################
// Copyright 2003-2004, Doug Enright, Ronald Fedkiw, Frank Losasso, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FALLING_DROP 
//##################################################################### 
#ifndef __FALLING_DROP__
#define __FALLING_DROP__

#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class FALLING_DROP:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;

    FALLING_DROP()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.WATER)
    {
        first_frame=0;last_frame=1000;
        frame_rate=24;
        restart=false;restart_frame=18;
        fluids_parameters.grid.Initialize(TV_INT(101,101,101),RANGE<TV>(TV(0,0,0),TV(1,1,1)));
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;
        output_directory="Falling_Drop/output";
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
    }
    
    ~FALLING_DROP() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() override
{
    GRID<TV>& grid=fluids_parameters.grid;
    VECTOR<T,3> center((T).5,(T).7,(T).5);
    T radius=(T).2;
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
        fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=min((grid.X(i,j,ij)-center).Magnitude()-radius, grid.y(j)-(T).21);
}
//#####################################################################
};      
}
#endif    


