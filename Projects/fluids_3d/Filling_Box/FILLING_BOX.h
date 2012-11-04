//#####################################################################
// Copyright 2002-2004, Sergey Koltakov, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//##################################################################### 
#ifndef __FILLING_BOX__
#define __FILLING_BOX__

#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class FILLING_BOX : public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;
using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;

CYLINDER<T> source;VECTOR<T,3> source_velocity;
    
    FILLING_BOX() 
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.WATER),source(VECTOR<T,3>(0,(T).75,(T).5),VECTOR<T,3>((T).05,(T).75,(T).5),(T).1),source_velocity((T)2,0,0)
    {
        first_frame=0;last_frame=1000;frame_rate=24;
        restart=false;restart_frame=0;
        fluids_parameters.grid.Initialize(TV_INT(76,51,51),RANGE<TV>(TV(0,0,0),TV(1.5,1,1)));
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=false;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;
        output_directory="Filling_Box/output";
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
    }

    ~FILLING_BOX() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    ARRAY<T,VECTOR<int,3> >::copy(fluids_parameters.grid.max_dx_dy_dz,fluids_parameters.particle_levelset_evolution.phi);
}
//#####################################################################
// Function Adjust_Phi_For_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) if(source.Inside(grid.X(i,j,ij))) 
        fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=-grid.dx;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    GRID<TV>& u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
    ARRAY<T,VECTOR<int,3> >& u=fluids_parameters.incompressible.projection.u,&v=fluids_parameters.incompressible.projection.v,&w=fluids_parameters.incompressible.projection.w;
    ARRAY<bool,VECTOR<int,3> >& psi_N_u=fluids_parameters.incompressible.projection.elliptic_solver->psi_N_u,&psi_N_v=fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v,&psi_N_w=fluids_parameters.incompressible.projection.elliptic_solver->psi_N_w;
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) for(int ij=0;ij<u_grid.mn;ij++) if(source.Inside(u_grid.X(i,j,ij))){psi_N_u(i,j,ij)=true;u(i,j,ij)=source_velocity.x;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) for(int ij=0;ij<v_grid.mn;ij++) if(source.Inside(v_grid.X(i,j,ij))){psi_N_v(i,j,ij)=true;v(i,j,ij)=source_velocity.y;}
    for(int i=0;i<w_grid.m;i++) for(int j=0;j<w_grid.n;j++) for(int ij=0;ij<w_grid.mn;ij++) if(source.Inside(w_grid.X(i,j,ij))){psi_N_w(i,j,ij)=true;w(i,j,ij)=source_velocity.z;}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time)
{
    GRID<TV>& p_grid=fluids_parameters.p_grid;CYLINDER<T> source_mask=source;source_mask.radius+=p_grid.dx*3;    
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,3> >(1,p_grid.m,1,p_grid.n,1,p_grid.mn);
    for(int i=0;i<p_grid.m;i++) for(int j=0;j<p_grid.n;j++) for(int ij=0;ij<p_grid.mn;ij++) 
        if(source_mask.Inside(VECTOR<T,3>(p_grid.X(i,j,ij)))) (*cell_centered_mask)(i,j,ij)=true;
}
//#####################################################################
};      
}
#endif
