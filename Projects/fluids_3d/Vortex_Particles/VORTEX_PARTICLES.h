//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTEX_PARTICLES 
//##################################################################### 
#ifndef __VORTEX_PARTICLES__
#define __VORTEX_PARTICLES__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include "RIVER_PHI_BOUNDARY.h"
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NR3.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class VORTEX_PARTICLES:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;

    VORTEX_PARTICLE_EVOLUTION_3D<T,RW> vortex_particle_evolution;
    T vorticity_magnitude;
    T particle_distance_from_surface;
    int vortex_count;
    T inflow;
    RANDOM_NR3 random;

    LEVELSET_IMPLICIT_SURFACE<T>* terrain_implicit;RIGID_BODY<TV> rigid_body;
    TRIANGULATED_SURFACE<T>* terrain_triangulated_surface;
    
    RIVER_PHI_BOUNDARY<T,T> custom_phi_boundary;

    VORTEX_PARTICLES()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.WATER)
    {
        first_frame=0;last_frame=1000;
        frame_rate=24;
        restart=false;restart_frame=18;
        //fluids_parameters.grid.Initialize(226,31,91,0,1.5,0,.2,.2,1);
        fluids_parameters.grid.Initialize(226,31,91,0,1.5,0,.2,.2,1);
        //fluids_parameters.grid.Initialize(61,31,31,0,1.5,0,1,0,1);
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=false;
        output_directory="Vortex_Particles/output-outflow";
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.write_ghost_values=true;
        // Setup vortex particle parameters
        fluids_parameters.use_body_force=true;
        fluids_parameters.confinement_parameter=(T)0;
        vortex_particle_evolution.particle_confinement_parameter=(T)1.00;
        vorticity_magnitude=1000; // 400;
        particle_distance_from_surface=.2;
        BOX_3D<T> box=fluids_parameters.grid.Domain();
        vortex_particle_evolution.Set_Radius(0.1);
        vortex_particle_evolution.force_scaling=(T)3e-4;
        vortex_particle_evolution.renormalize_vorticity_after_stretching_tilting=true;
        vortex_particle_evolution.Initialize(fluids_parameters.grid);
        inflow=1;

        terrain_implicit=new LEVELSET_IMPLICIT_SURFACE<T>(*(new GRID<TV>),*(new ARRAY<T,VECTOR<int,3> >));
        terrain_triangulated_surface=TRIANGULATED_SURFACE<T>::Create();FILE_UTILITIES::Read_From_File<RW>("terrain.tri",*terrain_triangulated_surface);
        FILE_UTILITIES::Read_From_File<RW>("terrain.phi",*terrain_implicit);
        custom_phi_boundary.implicit=terrain_implicit;
        custom_phi_boundary.inflow_height=.075;
        fluids_parameters.phi_boundary=&custom_phi_boundary;
    }
    
    ~VORTEX_PARTICLES() 
    {}

//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution.phi;GRID_3D<T>& p_grid=fluids_parameters.p_grid;
    ARRAY<T,VECTOR<int,3> >& p=fluids_parameters.incompressible.projection.p;
    for(int j=0;j<p_grid.n;j++) for(int ij=0;ij<p_grid.mn;ij++)fluids_parameters.incompressible.projection.elliptic_solver->psi_D(1,j,ij)=0; // unset the dirichlet boundary...
    for(int j=0;j<fluids_parameters.u_grid.n;j++) for(int ij=0;ij<fluids_parameters.u_grid.mn;ij++) if(!terrain_implicit->Lazy_Inside(fluids_parameters.p_grid.X(1,j,ij))) {
        fluids_parameters.incompressible.projection.elliptic_solver->psi_N_u(1,j,ij)=true;fluids_parameters.incompressible.projection.u(1,j,ij)=inflow;}
    // set pressure of each column of fluid on the boundary to rho g h, starting over if one or more air cells are found going down ..
    for(int ij=0;ij<p_grid.mn;ij++){
        T depth=0;T rho=1;
        for(int j=p_grid.n;j>=1;j--){
            fluids_parameters.incompressible.projection.elliptic_solver->psi_D(p_grid.m+1,j,ij)=true;
            if(phi(p_grid.m,j,ij)>0){p(p_grid.m+1,j,ij)=0;depth=0;}else{p(p_grid.m+1,j,ij)=fluids_parameters.gravity*depth*rho;depth+=p_grid.dy;}}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    rigid_body.Initialize_Triangulated_Surface(*terrain_triangulated_surface);rigid_body.Initialize_Implicit_Surface(*terrain_implicit);
    int id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body_And_Geometry(&rigid_body);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;
    Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);
}
//#####################################################################
// Function Update_Force
//#####################################################################
void Get_Body_Force(ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& force,const T time,const T dt)
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > V_ghost(grid,3,false);fluids_parameters.fluid_boundary->Fill_Ghost_Cells(grid,fluids_parameters.incompressible.V,V_ghost,dt,time);
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> >::copy(VECTOR<T,3>(),force);
    // Compute force
    vortex_particle_evolution.Compute_Body_Force(V_ghost,force,dt,time);
    // Advance Particles
    vortex_particle_evolution.Euler_Step(V_ghost,dt,time);
    // delete particles that are outside the water
    VORTICITY_PARTICLES<T,VECTOR<T,3> >& vorticity_particles=vortex_particle_evolution.vorticity_particles;
    for(int p=vorticity_particles.Size();p>=1;p--) if(fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.Phi(vorticity_particles.X(p))>0)
        vorticity_particles.Delete_Particle(p);
    // seed particles
    return;
    VECTOR<T,3> cell_upper=(T).5*VECTOR<T,3>(grid.dx,grid.dy,grid.dz),cell_lower=-cell_upper;
    int i=1;
    for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) if(fluids_parameters.particle_levelset_evolution.phi(i,j,ij)<0)
        if(random.Get_Uniform_Number((T)0,(T)1)<(T).0003){
            int particle_id=vorticity_particles.Add_Element(); 
            vorticity_particles.X(particle_id)=grid.X(i,j,ij)+random.Get_Uniform_Vector(cell_lower,cell_upper);
            int sign=random.Get_Uniform_Integer(0,1);if(sign==0)sign=-1;
            vorticity_particles.vorticity(particle_id)=VECTOR<T,3>(0,(T)sign*vorticity_magnitude,0);
            std::cout<<"ADDING DEFORMABLE_PARTICLES"<<std::endl;}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Velocities()
{
    for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++) for(int ij=0;ij<fluids_parameters.grid.mn;ij++) if(fluids_parameters.particle_levelset_evolution.phi(i,j,ij)<0)
        fluids_parameters.incompressible.V(i,j,ij)=VECTOR<T,3>(inflow,0,0);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    GRID<TV>& grid=fluids_parameters.grid;
//    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++) fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=grid.y(j)-(T).1;
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++) fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=grid.dx;;
    for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) if(!terrain_implicit->Lazy_Inside(grid.X(1,j,ij)) && grid.y(j) < custom_phi_boundary.inflow_height) fluids_parameters.particle_levelset_evolution.phi(1,j,ij)=-grid.dx;

}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    GRID<TV> &u_grid=fluids_parameters.u_grid,&v_grid=fluids_parameters.v_grid,&w_grid=fluids_parameters.w_grid;
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) for(int ij=0;ij<u_grid.mn;ij++) if(terrain_implicit->Lazy_Inside(u_grid.X(i,j,ij))){
        fluids_parameters.incompressible.projection.u(i,j,ij)=0;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_u(i,j,ij)=true;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) for(int ij=0;ij<v_grid.mn;ij++) if(terrain_implicit->Lazy_Inside(v_grid.X(i,j,ij))){
        fluids_parameters.incompressible.projection.v(i,j,ij)=0;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v(i,j,ij)=true;}
    for(int i=0;i<w_grid.m;i++) for(int j=0;j<w_grid.n;j++) for(int ij=0;ij<w_grid.mn;ij++) if(terrain_implicit->Lazy_Inside(w_grid.X(i,j,ij))){
        fluids_parameters.incompressible.projection.w(i,j,ij)=0;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_w(i,j,ij)=true;}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    FILE_UTILITIES::Create_Directory(output_directory);
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Write_Output_Files(frame);
    vortex_particle_evolution.Write_Output_Files(output_directory,frame);
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time)
{
    GRID<TV>& p_grid=fluids_parameters.p_grid;
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,3> >(1,p_grid.m,1,p_grid.n,1,p_grid.mn);
    for(int j=0;j<p_grid.n;j++) for(int ij=0;ij<p_grid.mn;ij++)
        (*cell_centered_mask)(1,j,ij)=(*cell_centered_mask)(2,j,ij)=true;
}
//#####################################################################
};      
}
#endif    


