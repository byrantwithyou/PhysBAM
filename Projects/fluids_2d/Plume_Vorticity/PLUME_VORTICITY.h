//#####################################################################
// Copyright 2001-2004, Frank Losasso, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLUME
//##################################################################### 
#ifndef __PLUME_VORTICITY__
#define __PLUME_VORTICITY__

#include <PhysBAM_Tools/Particles_Interpolation/SCATTERED_INTERPOLATION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>

namespace PhysBAM{

template <class T,class RW=T>
class PLUME_VORTICITY:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_output_files;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;

    T rho,rho_bottom,rho_top,buoyancy_constant;
    BOX_2D<T> source_domain;
    VORTICITY_PARTICLES<T,VECTOR_2D<T> > vorticity_positive_particles,vorticity_negative_particles;
    T particle_vorticity_confinement_parameter;
    T particle_vorticity_initial_amplification_factor;
    T particle_vorticity_minimum_density;
    RANDOM_NR3 random;
    SCATTERED_INTERPOLATION<T> scattered_interpolation;
    ARRAY<VECTOR_3D<T> ,VECTOR<int,2> > grid_vorticity,grid_vorticity_raw_particles,grid_vorticity_particles;

    PLUME_VORTICITY()
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.SMOKE)
    {
        fluids_parameters.grid.Initialize(100,100,0,1,0,1);

        first_frame=0;
        last_frame=512;
        frame_rate=24;
        fluids_parameters.cfl=3;
        fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.kolmogorov=(T)0;
        rho=1;rho_bottom=1;rho_top=(T).65;buoyancy_constant=0;fluids_parameters.gravity=0;
        fluids_parameters.temperature_container.Set_Cooling_Constant(0);fluids_parameters.temperature_products=3000;
        output_directory="Plume_Vorticity/output";
        write_output_files=true;fluids_parameters.write_debug_data=true;
        source_domain=BOX_2D<T>((T).45,(T).55,(T)0,(T).1);
        fluids_parameters.solve_neumann_regions=false;
        // Setup vortex particle parameters
        fluids_parameters.use_vorticity_confinement=false; // use our particle instead... 
        fluids_parameters.use_body_force=true;
        fluids_parameters.confinement_parameter=(T).05;
        particle_vorticity_confinement_parameter=(T).5;
        particle_vorticity_minimum_density=(T).05;
        particle_vorticity_initial_amplification_factor=(T)3;
        scattered_interpolation.Set_Radius_Of_Influence(3*fluids_parameters.grid.dx);scattered_interpolation.Use_Tent_Weights();
    }
    
    ~PLUME_VORTICITY()
    {}
    
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
void Adjust_Density_And_Temperature_With_Sources(const T time)
{
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++)if(source_domain.Lazy_Inside(fluids_parameters.grid.X(i,j))){
        fluids_parameters.density_container.density_2d(i,j)=rho;fluids_parameters.temperature_container.temperature_2d(i,j)=fluids_parameters.temperature_products;}
    // keep density >= 0 and T >=0
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++){
        fluids_parameters.density_container.density_2d(i,j)=max((T)0,fluids_parameters.density_container.density_2d(i,j));
        fluids_parameters.temperature_container.temperature_2d(i,j)=max((T)fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature_2d(i,j));}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    for(int i=1;i<=fluids_parameters.u_grid.m;i++)for(int j=1;j<=fluids_parameters.u_grid.n;j++)if(source_domain.Lazy_Inside(fluids_parameters.u_grid.X(i,j)))
        fluids_parameters.incompressible.projection.u(i,j)=0;
    for(int i=1;i<=fluids_parameters.v_grid.m;i++)for(int j=1;j<=fluids_parameters.v_grid.n;j++)if(source_domain.Lazy_Inside(fluids_parameters.v_grid.X(i,j))){
        fluids_parameters.incompressible.projection.v(i,j)=(T).2;fluids_parameters.incompressible.projection.elliptic_solver->psi_N_v(i,j)=true;}
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
void Get_Body_Force(ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >& force,const T time)
{
    std::cout<<"body force"<<std::endl;
    ARRAY<VECTOR_2D<T> ,VECTOR<int,2> > V_ghost(fluids_parameters.grid,3,false);fluids_parameters.fluid_boundary->Fill_Ghost_Cells(fluids_parameters.grid,fluids_parameters.incompressible.V,V_ghost,0,time);
    grid_vorticity.Resize(fluids_parameters.grid,2,false,false);grid_vorticity_raw_particles.Resize(fluids_parameters.grid,2,false,false);grid_vorticity_particles.Resize(fluids_parameters.grid,2,false,false);
    ARRAY<T,VECTOR<int,2> > grid_vorticity_magnitude(fluids_parameters.grid,2,false),grid_vorticity_particles_magnitude(fluids_parameters.grid,2,false);
    // Form grid vorticity
    VORTICITY_UNIFORM<TV>::Vorticity(fluids_parameters.grid,V_ghost,grid_vorticity_magnitude);
    for(int i=grid_vorticity.m_start;i<=grid_vorticity.m_end;i++)for(int j=grid_vorticity.n_start;j<=grid_vorticity.n_end;j++){
        grid_vorticity(i,j)=VECTOR_3D<T>(0,0,grid_vorticity_magnitude(i,j));grid_vorticity_magnitude(i,j)=fabs(grid_vorticity_magnitude(i,j));}
    // Transfer particle vorticity to particle grid vorticity
    ARRAY<VECTOR_3D<T> ,VECTOR<int,2> > grid_vorticity_positive_particles(fluids_parameters.grid,2,false),grid_vorticity_negative_particles(fluids_parameters.grid,2,false);
    scattered_interpolation.Transfer_To_Grid(vorticity_positive_particles.X.array,vorticity_positive_particles.vorticity.array,fluids_parameters.grid,grid_vorticity_positive_particles);
    scattered_interpolation.Transfer_To_Grid(vorticity_negative_particles.X.array,vorticity_negative_particles.vorticity.array,fluids_parameters.grid,grid_vorticity_negative_particles);
    for(int i=0;i<=fluids_parameters.grid.m+1;i++) for(int j=0;j<=fluids_parameters.grid.n+1;j++){
        if(grid_vorticity(i,j).z>=0){
            grid_vorticity_particles(i,j)=clamp_min(grid_vorticity_positive_particles(i,j)-grid_vorticity(i,j),VECTOR_3D<T>(0,0,0));
            grid_vorticity_raw_particles(i,j)=grid_vorticity_positive_particles(i,j);
            grid_vorticity_particles_magnitude(i,j)=grid_vorticity_particles(i,j).Magnitude();}
        else{
            grid_vorticity_particles(i,j)=clamp_max(grid_vorticity_negative_particles(i,j)-grid_vorticity(i,j),VECTOR_3D<T>(0,0,0));
            grid_vorticity_raw_particles(i,j)=grid_vorticity_negative_particles(i,j);
            grid_vorticity_particles_magnitude(i,j)=grid_vorticity_particles(i,j).Magnitude();}}
    // Compute confinement force
    T one_over_two_dx=1/(2*fluids_parameters.grid.dx),one_over_two_dy=1/(2*fluids_parameters.grid.dy);
    for(int i=force.m_start;i<=force.m_end;i++) for(int j=force.n_start;j<=force.n_end;j++){
        VECTOR_3D<T> vortex_normal_vector((grid_vorticity_magnitude(i+1,j)-grid_vorticity_magnitude(i-1,j))*one_over_two_dx,(grid_vorticity_magnitude(i,j+1)-grid_vorticity_magnitude(i,j-1))*one_over_two_dy,0);
        VECTOR_3D<T> particle_vortex_normal_vector((grid_vorticity_particles_magnitude(i+1,j)-grid_vorticity_particles_magnitude(i-1,j))*one_over_two_dx,(grid_vorticity_particles_magnitude(i,j+1)-grid_vorticity_particles_magnitude(i,j-1))*one_over_two_dy,0);
        vortex_normal_vector.Robust_Normalize(1e-6,VECTOR_3D<T>());particle_vortex_normal_vector.Robust_Normalize(1e-6,VECTOR_3D<T>());
        VECTOR_3D<T> force_3d=(T)fluids_parameters.confinement_parameter*VECTOR_3D<T>::Cross_Product(vortex_normal_vector,grid_vorticity(i,j))
            +(T)particle_vorticity_confinement_parameter*VECTOR_3D<T>::Cross_Product(particle_vortex_normal_vector,grid_vorticity_particles(i,j));
        force(i,j)=VECTOR_2D<T>(force_3d.x,force_3d.y);}
    // Add particles
    int add_count=0;
    for(int i=1;i<=fluids_parameters.grid.m;i++)for(int j=1;j<=fluids_parameters.grid.n;j++){
        if(time>(T)1/24 && random.Get_Uniform_Number((T)0,(T)1)<(T).01 && grid_vorticity(i,j).Magnitude_Squared()>1e-6 && fluids_parameters.density_container.density_2d(i,j)>(T)particle_vorticity_minimum_density){
            add_count++;
            if(grid_vorticity(i,j).z>=0){
                int particle_id=vorticity_positive_particles.array_collection->Add_Element(); 
                vorticity_positive_particles.X(particle_id)=fluids_parameters.grid.X(i,j)+random.Get_Uniform_Vector(VECTOR_2D<T>(-fluids_parameters.grid.dx,fluids_parameters.grid.dx),VECTOR_2D<T>(-fluids_parameters.grid.dy,fluids_parameters.grid.dy));
                vorticity_positive_particles.vorticity(particle_id)=particle_vorticity_initial_amplification_factor*grid_vorticity(i,j);}
            else{
                int particle_id=vorticity_negative_particles.array_collection->Add_Element(); 
                vorticity_negative_particles.X(particle_id)=fluids_parameters.grid.X(i,j)+random.Get_Uniform_Vector(VECTOR_2D<T>(-fluids_parameters.grid.dx,fluids_parameters.grid.dx),VECTOR_2D<T>(-fluids_parameters.grid.dy,fluids_parameters.grid.dy));
                vorticity_negative_particles.vorticity(particle_id)=particle_vorticity_initial_amplification_factor*grid_vorticity(i,j);}}}
    printf("Added %d particles\n",add_count);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame)
{
    printf("Advecting vorticity particles\n");
    T dt=1/(T)frame_rate;
    vorticity_positive_particles.X.Euler_Step(vorticity_positive_particles,fluids_parameters.grid,fluids_parameters.incompressible.V,dt);
    vorticity_negative_particles.X.Euler_Step(vorticity_negative_particles,fluids_parameters.grid,fluids_parameters.incompressible.V,dt);
    int old_count=vorticity_positive_particles.array_collection->Size()+vorticity_negative_particles.array_collection->Size();
    vorticity_positive_particles.Delete_Particles_Outside_Grid(fluids_parameters.grid,vorticity_positive_particles.X.array);
    vorticity_negative_particles.Delete_Particles_Outside_Grid(fluids_parameters.grid,vorticity_negative_particles.X.array);
    int deleted_count=old_count-vorticity_positive_particles.array_collection->Size()-vorticity_negative_particles.array_collection->Size();
    printf("Deleted %d particles\n",deleted_count);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const T time,const int frame)
{
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Write_Output_Files(time,frame);
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/vorticity_particles.%d",output_directory.c_str(),frame),vorticity_positive_particles,vorticity_negative_particles);
    // Output non-zero portion of vorticity for easy visualization.
    ARRAY<T,VECTOR<int,2> > vorticity_z(fluids_parameters.grid,0,false);
    for(int i=1;i<=fluids_parameters.grid.m;i++) for(int j=1;j<=fluids_parameters.grid.n;j++) vorticity_z(i,j)=grid_vorticity(i,j).z;
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/grid_vorticity.%d",output_directory.c_str(),frame),vorticity_z);
    for(int i=1;i<=fluids_parameters.grid.m;i++) for(int j=1;j<=fluids_parameters.grid.n;j++) vorticity_z(i,j)=grid_vorticity_raw_particles(i,j).z;
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/grid_vorticity_raw_particles.%d",output_directory.c_str(),frame),vorticity_z);
    for(int i=1;i<=fluids_parameters.grid.m;i++) for(int j=1;j<=fluids_parameters.grid.n;j++) vorticity_z(i,j)=grid_vorticity_particles(i,j).z;
    FILE_UTILITIES::Write_To_File<RW>(STRING_UTILITIES::string_sprintf("%s/grid_vorticity_particles.%d",output_directory.c_str(),frame),vorticity_z);
}
//#####################################################################
};      
}
#endif


