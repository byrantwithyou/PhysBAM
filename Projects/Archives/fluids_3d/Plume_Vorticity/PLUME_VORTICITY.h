//#####################################################################
// Copyright 2001-2004, Ron Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLUME_VORTICITY
//#####################################################################
#ifndef __PLUME_VORTICITY__
#define __PLUME_VORTICITY__

#include <PhysBAM_Tools/Particles_Interpolation/SCATTERED_INTERPOLATION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>

namespace PhysBAM{

template <class T,class RW=T>
class PLUME_VORTICITY:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;
    T source_vorticity_magnitude;
    T particle_vorticity_minimum_density;
    RANDOM_NR3 random;
    VORTEX_PARTICLE_EVOLUTION_3D<T,RW> vortex_particle_evolution;
    T rho,rho_bottom,rho_top,buoyancy_constant;
    BOX_3D<T> source_domain;
    
    PLUME_VORTICITY()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.SMOKE)
    {
        fluids_parameters.grid.Initialize(100,100,100,0,1,0,1,0,1);
        first_frame=0;
        last_frame=200;
        frame_rate=24;
        fluids_parameters.cfl=3;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.use_vorticity_confinement=false;fluids_parameters.confinement_parameter=(T).3;
        fluids_parameters.kolmogorov=(T)0;
        rho=(T)1;rho_bottom=(T)1;rho_top=(T).65;buoyancy_constant=(T)0;fluids_parameters.gravity=(T)0;
        fluids_parameters.temperature_container.Set_Cooling_Constant((T)1);fluids_parameters.temperature_products=(T)3000;
        write_output_files=true;fluids_parameters.write_debug_data=true;
        source_domain=BOX_3D<T>((T).45,(T).55,(T)0,(T).1,(T).45,(T).55);
        fluids_parameters.Initialize_Domain_Boundary_Conditions();
        output_directory="Plume_Vorticity/output";
        // Setup vortex particle parameters
        fluids_parameters.use_vorticity_confinement=false; // use our particle instead... 
        fluids_parameters.use_body_force=true;
        fluids_parameters.confinement_parameter=(T)0.3;
        vortex_particle_evolution.particle_confinement_parameter=(T)1.00;
        vortex_particle_evolution.renormalize_vorticity_after_stretching_tilting=true;
        particle_vorticity_minimum_density=(T).9;
        source_vorticity_magnitude=(T)100;
        vortex_particle_evolution.Initialize(fluids_parameters.grid);
    }
    
    ~PLUME_VORTICITY()
    {}
    
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time)
{
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++)for(int ij=1;ij<=fluids_parameters.grid.mn;ij++)if(source_domain.Lazy_Inside(fluids_parameters.grid.X(i,j,ij))){
        fluids_parameters.density_container.density_3d(i,j,ij)=rho;fluids_parameters.temperature_container.temperature_3d(i,j,ij)=fluids_parameters.temperature_products;}
    // keep density >= 0 and T >=0
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++)for(int ij=1;ij<=fluids_parameters.grid.mn;ij++){
        fluids_parameters.density_container.density_3d(i,j,ij)=max((T)0,fluids_parameters.density_container.density_3d(i,j,ij));
        fluids_parameters.temperature_container.temperature_3d(i,j,ij)=max((T)fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature_3d(i,j,ij));}

}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    for(int i=1;i<=fluids_parameters.u_grid.m;i++)for(int j=1;j<=fluids_parameters.u_grid.n;j++)for(int ij=1;ij<=fluids_parameters.u_grid.mn;ij++)
        if(source_domain.Lazy_Inside(fluids_parameters.u_grid.X(i,j,ij))) fluids_parameters.incompressible.projection.u(i,j,ij)=0;
    for(int i=1;i<=fluids_parameters.v_grid.m;i++)for(int j=1;j<=fluids_parameters.v_grid.n;j++)for(int ij=1;ij<=fluids_parameters.v_grid.n;ij++)
        if(source_domain.Lazy_Inside(fluids_parameters.v_grid.X(i,j,ij))){fluids_parameters.incompressible.projection.v(i,j,ij)=(T)1;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v(i,j,ij)=true;}
    for(int i=1;i<=fluids_parameters.w_grid.m;i++)for(int j=1;j<=fluids_parameters.w_grid.n;j++)for(int ij=1;ij<=fluids_parameters.w_grid.mn;ij++)
        if(source_domain.Lazy_Inside(fluids_parameters.w_grid.X(i,j,ij))) fluids_parameters.incompressible.projection.w(i,j,ij)=0;
}
//#####################################################################
// Function Update_Forces
//#####################################################################
void Get_Body_Force(ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& force,const T dt,const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > V_ghost(grid,3,false);fluids_parameters.fluid_boundary->Fill_Ghost_Cells(grid,fluids_parameters.incompressible.V,V_ghost,dt,time);
    // Compute force
    vortex_particle_evolution.Compute_Body_Force(V_ghost,force,dt,time);
    // Add particles
    int add_count=0;VORTICITY_PARTICLES<T,VECTOR<T,3> >& vorticity_particles=vortex_particle_evolution.vorticity_particles;
    VECTOR<T,3> cell_upper=(T).5*VECTOR<T,3>(grid.dx,grid.dy,grid.dz),cell_lower=-cell_upper;
    for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++)
        if(source_domain.Lazy_Inside(grid.X(i,j,ij)) && time>(T)1/24 && random.Get_Uniform_Number((T)0,(T)1)<(T).02){
            add_count++;
            int particle_id=vorticity_particles.array_collection->Add_Element(); 
            vorticity_particles.X(particle_id)=grid.X(i,j,ij)+random.Get_Uniform_Vector(cell_lower,cell_upper);
            vorticity_particles.vorticity(particle_id)=(T)source_vorticity_magnitude*VECTOR<T,3>::Cross_Product((vorticity_particles.X(particle_id)-source_domain.Center()).Robust_Normalized(),VECTOR<T,3>(0,1,0)).Normalized();}
            //vorticity_particles.vorticity(particle_id)=(T)source_vorticity_magnitude*vortex_particle_evolution.grid_vorticity(i,j,ij).Robust_Normalized();}
    printf("Added %d particles\n",add_count);
    // Advance Particles
    vortex_particle_evolution.Euler_Step(V_ghost,dt,time);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Write_Output_Files(frame);
    vortex_particle_evolution.Write_Output_Files(output_directory,frame);
}
//#####################################################################
};      
}
#endif


