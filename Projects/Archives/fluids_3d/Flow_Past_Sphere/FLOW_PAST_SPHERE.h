//#####################################################################
// Copyright 2001-2004, Ron Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOW_PAST_SPHERE
//#####################################################################
#ifndef __FLOW_PAST_SPHERE__
#define __FLOW_PAST_SPHERE__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/VORTEX_PARTICLE_EVOLUTION_3D.h>
#include <Level_Sets/EXTRAPOLATION_3D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>

namespace PhysBAM{

template <class T,class RW=T>
class FLOW_PAST_SPHERE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_output_files;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;

    T rho,rho_bottom,rho_top,buoyancy_constant;
    BOX_3D<T> source_domain;
    ARRAY<T,VECTOR<int,3> > phi_object;
    LEVELSET_3D<GRID<TV> > levelset_object;
    
    T source_vorticity_magnitude;
    T particle_vorticity_minimum_density;
    RANDOM_NR3 random;
    VORTEX_PARTICLE_EVOLUTION_3D<T,RW> vortex_particle_evolution;
    bool use_vortex_particles;


    FLOW_PAST_SPHERE()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.SMOKE),levelset_object(fluids_parameters.grid,phi_object)
    {
        //fluids_parameters.grid.Initialize(61,61,61,0,1,0,1,0,1);
        fluids_parameters.grid.Initialize(201,201,201,0,1,0,1,0,1);
        first_frame=0;
        last_frame=400;
        frame_rate=60;
        fluids_parameters.cfl=0.9;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[3][1]=fluids_parameters.domain_walls[3][2]=false;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T).15;
        fluids_parameters.kolmogorov=(T)0;
        rho=(T)1;rho_bottom=(T)1;rho_top=(T).65;buoyancy_constant=(T).01;fluids_parameters.gravity=(T)0;fluids_parameters.use_body_force=true;
        fluids_parameters.temperature_container.Set_Cooling_Constant((T)1);fluids_parameters.temperature_products=(T)3000;
        write_output_files=true;fluids_parameters.write_debug_data=false;
        source_domain=BOX_3D<T>((T).45,(T).55,(T)0,(T).05,(T).45,(T).55);
        fluids_parameters.Initialize_Domain_Boundary_Conditions();
        output_directory="Flow_Past_Sphere/output";

        // setup vortex particles
        fluids_parameters.use_body_force=true;
        use_vortex_particles=true;
        vortex_particle_evolution.particle_confinement_parameter=(T)1.00;
        vortex_particle_evolution.renormalize_vorticity_after_stretching_tilting=true;
        particle_vorticity_minimum_density=(T).9;
        source_vorticity_magnitude=(T)200;
        vortex_particle_evolution.Set_Radius((T).02);
        vortex_particle_evolution.Initialize(fluids_parameters.grid);
        vortex_particle_evolution.force_scaling=1e-4;

    }
    
    ~FLOW_PAST_SPHERE()
    {}
    
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>((data_directory+"/Rigid_Bodies/sphere").c_str(), (T).1, true, true, false);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR<T,3>((T).5,(T).35,(T).5);
    solids_parameters.rigid_body_parameters.list(id)->is_static=true;
    Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    phi_object.Resize(fluids_parameters.grid,3);
    for(int i=-2;i<=fluids_parameters.grid.m+3;i++) for(int j=-2;j<=fluids_parameters.grid.n+3;j++) for(int k=-2;k<=fluids_parameters.grid.mn+3;k++) 
        phi_object(i,j,k)=solids_parameters.rigid_body_parameters.list(1)->Implicit_Surface_Extended_Value(fluids_parameters.grid.X(i,j,k));
    levelset_object.Compute_Cell_Minimum_And_Maximum();
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(solids_parameters.rigid_body_parameters.list.Number_Of_Elements()<1) return;
    for(int i=1;i<=fluids_parameters.u_grid.m;i++)for(int j=1;j<=fluids_parameters.u_grid.n;j++)for(int ij=1;ij<=fluids_parameters.u_grid.mn;ij++)
        if(levelset_object.Lazy_Inside(fluids_parameters.u_grid.X(i,j,ij))){fluids_parameters.incompressible.projection.u(i,j,ij)=0;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_u(i,j,ij)=true;}
    for(int i=1;i<=fluids_parameters.v_grid.m;i++)for(int j=1;j<=fluids_parameters.v_grid.n;j++)for(int ij=1;ij<=fluids_parameters.v_grid.n;ij++)
        if(levelset_object.Lazy_Inside(fluids_parameters.v_grid.X(i,j,ij))){fluids_parameters.incompressible.projection.v(i,j,ij)=0;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v(i,j,ij)=true;}
    for(int i=1;i<=fluids_parameters.w_grid.m;i++)for(int j=1;j<=fluids_parameters.w_grid.n;j++)for(int ij=1;ij<=fluids_parameters.w_grid.mn;ij++)
        if(levelset_object.Lazy_Inside(fluids_parameters.w_grid.X(i,j,ij))){fluids_parameters.incompressible.projection.w(i,j,ij)=0;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_w(i,j,ij)=true;}

}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time){
    // sources
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++)for(int ij=1;ij<=fluids_parameters.grid.mn;ij++)if(source_domain.Lazy_Inside(fluids_parameters.grid.X(i,j,ij))){
        fluids_parameters.density_container.density_3d(i,j,ij)=rho;fluids_parameters.temperature_container.temperature_3d(i,j,ij)=fluids_parameters.temperature_products;}
    // objects
    phi_object*=(T)-1; // because extrapolation will extrapolate into positive phi. changed back below.
    EXTRAPOLATION_3D<T,T> extrapolate(fluids_parameters.grid,phi_object,fluids_parameters.density_container.density_3d);extrapolate.Set_Band_Width(3);extrapolate.Extrapolate();
    phi_object*=(T)-1; // revert to original phi
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
        if(source_domain.Lazy_Inside(fluids_parameters.v_grid.X(i,j,ij))){fluids_parameters.incompressible.projection.v(i,j,ij)=(T).5;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v(i,j,ij)=true;}
    for(int i=1;i<=fluids_parameters.w_grid.m;i++)for(int j=1;j<=fluids_parameters.w_grid.n;j++)for(int ij=1;ij<=fluids_parameters.w_grid.mn;ij++)
    if(source_domain.Lazy_Inside(fluids_parameters.w_grid.X(i,j,ij))) fluids_parameters.incompressible.projection.w(i,j,ij)=0;
}
//#####################################################################
// Function Update_Forces
//#####################################################################
void Get_Body_Force(ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& force,const T dt,const T time)
{
    // control buoyancy force
    for(int i=1;i<=fluids_parameters.grid.m;i++) for(int j=1;j<=fluids_parameters.grid.n;j++) for(int ij=1;ij<=fluids_parameters.grid.mn;ij++){ // y-direction forces only
        T rho_atm=rho_bottom+(rho_top-rho_bottom)*(fluids_parameters.grid.y(j)-fluids_parameters.grid.ymin)/(fluids_parameters.grid.ymax-fluids_parameters.grid.ymin);
        if(fluids_parameters.grid.y(j)<=2)rho_atm=1;
        T difference=fluids_parameters.temperature_container.temperature_3d(i,j,ij)-fluids_parameters.temperature_container.ambient_temperature;
        if(difference>0)force(i,j,ij).y=fluids_parameters.buoyancy_constant*difference;}

    if(use_vortex_particles){
        GRID<TV>& grid=fluids_parameters.grid;
        // Vortex particles
        ARRAY<VECTOR<T,3> ,VECTOR<int,3> > V_ghost(grid,3,false);fluids_parameters.fluid_boundary->Fill_Ghost_Cells(grid,fluids_parameters.incompressible.V,V_ghost,dt,time);
        // Compute force
        vortex_particle_evolution.Compute_Body_Force(V_ghost,force,dt,time);
        // Add particles
//        if(time <= explosion_end_time){
            int add_count=0;VORTICITY_PARTICLES<T,VECTOR<T,3> >& vorticity_particles=vortex_particle_evolution.vorticity_particles;
            VECTOR<T,3> cell_upper=(T).5*VECTOR<T,3>(grid.dx,grid.dy,grid.dz),cell_lower=-cell_upper;
            for(int i=1;i<=grid.m;i++) for(int j=1;j<=grid.n;j++) for(int ij=1;ij<=grid.mn;ij++)
                if(source_domain.Lazy_Inside(grid.X(i,j,ij)) && time>(T)1/24 && random.Get_Uniform_Number((T)0,(T)1)<(T).005){
                    add_count++;
                    int particle_id=vorticity_particles.array_collection->Add_Element(); 
                    vorticity_particles.X(particle_id)=grid.X(i,j,ij)+random.Get_Uniform_Vector(cell_lower,cell_upper);
                    vorticity_particles.vorticity(particle_id)=(T)source_vorticity_magnitude*VECTOR<T,3>::Cross_Product(VECTOR<T,3>(0,1,0),(vorticity_particles.X(particle_id)-source_domain.Center()).Robust_Normalized()).Normalized();
            printf("Added %d particles\n",add_count);}
        // Advance Particles
        vortex_particle_evolution.Euler_Step(V_ghost,dt,time);}

}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
virtual void Read_Output_Files_Fluids(const int frame)
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Read_Output_Files_Fluids(frame);
    if(use_vortex_particles) vortex_particle_evolution.Read_Output_Files(output_directory,frame);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Write_Output_Files(frame);
    if(use_vortex_particles) vortex_particle_evolution.Write_Output_Files(output_directory,frame);
}
//#####################################################################
};      
}
#endif


