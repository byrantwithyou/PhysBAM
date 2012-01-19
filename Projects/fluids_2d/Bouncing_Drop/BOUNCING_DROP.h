//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNCING_DROP
//#####################################################################
#ifndef __BOUNCING_DROP__
#define __BOUNCING_DROP__

#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class BOUNCING_DROP:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_frame_title;

    BOUNCING_DROP():SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER)
    {
        fluids_parameters.grid.Initialize(64,64,0,1,0,1);
        first_frame=0;last_frame=2000;
        frame_rate=24;
        restart=false;restart_frame=18;
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=false; 
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.particle_half_bandwidth=1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fluids_parameters.reseeding_frame_rate=1;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.incompressible_iterations=20;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;
        output_directory="Bouncing_Drop/output";
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=true;
        
        fluids_parameters.viscosity=(T)100;
        fluids_parameters.use_strain=true;fluids_parameters.write_strain=true;
        fluids_parameters.elastic_modulus=10000;
        fluids_parameters.plasticity_alpha=0;
        fluids_parameters.cfl/=3;
        fluids_parameters.levelset_substeps=1;
        frame_rate*=10;
    }
    
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    VECTOR_2D<T> center((T).5,(T).5);
    T radius=(T).2;
    for(int i=0;i<fluids_parameters.grid.m;i++)for(int j=0;j<fluids_parameters.grid.n;j++)
        fluids_parameters.particle_levelset_evolution.phi(i,j)=(fluids_parameters.grid.X(i,j)-center).Magnitude()-radius;
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities()
{
    VECTOR_2D<T> center((T).5,(T).7);
    T scale=10;
    for(int i=0;i<fluids_parameters.grid.m;i++)for(int j=0;j<fluids_parameters.grid.n;j++)
        //    V(i,j)=scale*VECTOR_2D<T>(grid.x(i)-center.x,0);
        fluids_parameters.particle_levelset_evolution.V(i,j)=VECTOR_2D<T>(0,-scale);
}
//#####################################################################
};      
}
#endif
