//#####################################################################
// Copyright 2002-2004, Doug Enright, Ronald Fedkiw, Duc Nguyen, Sergey Koltakov, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FALLING_DROP
//#####################################################################
#ifndef __FALLING_DROP__
#define __FALLING_DROP__

#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class FALLING_DROP:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_frame_title;

    FALLING_DROP()
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER)
    {
        fluids_parameters.grid.Initialize(64,64,0,1,0,1);
        first_frame=0;last_frame=2000;
        frame_rate=24;
        restart=false;restart_frame=18;
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=true; 
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.particle_half_bandwidth=1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fluids_parameters.reseeding_frame_rate=1;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.viscosity=0;//(T)(.001137*2e7);
        fluids_parameters.incompressible_iterations=20;
        write_frame_title=true;;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;
        output_directory="Falling_Drop/output";
        //gravity=0;
        //gravity_direction.x=.5;gravity_direction.Normalize();
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.scalar_substeps=(T)2.9;
        fluids_parameters.second_order_cut_cell_method=true;
    }

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    VECTOR_2D<T> center((T).5,(T).7);
    T radius=(T).2;
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++)
        fluids_parameters.particle_levelset_evolution.phi(i,j)=min((fluids_parameters.grid.X(i,j)-center).Magnitude()-radius, fluids_parameters.grid.y(j)-(T).21);
//        particle_levelset_evolution.phi(i,j)=grid.y(j)-(T).21;
//        particle_levelset_evolution.phi(i,j)=(T).79-grid.x(i);
//        particle_levelset_evolution.phi(i,j)=BOX_2D<T>(0.3,0.7,0.7,1).Signed_Distance(grid.X(i,j));
}
//#####################################################################
};    
}
#endif
