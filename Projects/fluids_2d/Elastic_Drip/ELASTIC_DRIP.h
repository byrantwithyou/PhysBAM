//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ELASTIC_DRIP
//#####################################################################
#ifndef __ELASTIC_DRIP__
#define __ELASTIC_DRIP__

#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T,class RW=T>
class ELASTIC_DRIP:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::Adjust_Phi_With_Sources;
    using BASE::Get_Source_Reseed_Mask;using BASE::Get_Source_Velocities;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::data_directory;

    BOX_2D<T> source;
    VECTOR<T,2> source_velocity;
    
    ELASTIC_DRIP():SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(1,fluids_parameters.WATER)
    {
        fluids_parameters.grid.Initialize(2*64,2*64,0,1,0,1);
        first_frame=0;last_frame=2000;
        frame_rate=24;
        restart=false;restart_frame=18;
        output_directory="Elastic_Drip/output";
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=false; 
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.reseeding_frame_rate=1;
        //bias_towards_negative_particles=true;
        fluids_parameters.incompressible_iterations=20;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        
        source=BOX_2D<T>(.4,.6,.9,1);
        source_velocity=VECTOR<T,2>(0,-.3);
        
        //fluids_parameters.viscosity=(T)100;
        fluids_parameters.use_strain=true;fluids_parameters.write_strain=true;
        fluids_parameters.elastic_modulus=1000;
        fluids_parameters.plasticity_alpha=1;
        fluids_parameters.cfl/=20;
    }

//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection()    PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)
        phi(i,j)=1;
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<T,VECTOR<int,2> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++){
        VECTOR<T,2> X=grid.X(i,j);
        if(source.Lazy_Inside(X))phi(i,j)=min(phi(i,j),source.Signed_Distance(X));}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& cell_centered_mask,const T time) PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,2> >(grid);

    T padding=3*grid.max_dx_dy;
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)if(!source.Outside(grid.X(i,j),padding))(*cell_centered_mask)(i,j)=true;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time) PHYSBAM_OVERRIDE
{
    for(FACE_ITERATOR iterator(fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(iterator.Location())){
        int axis=iterator.Axis();fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(iterator.Face_Index())=true;
        fluids_parameters.incompressible->projection.face_velocities.Component(axis)(iterator.Face_Index())=source_velocity[axis];}
}
//#####################################################################
};      
}
#endif
