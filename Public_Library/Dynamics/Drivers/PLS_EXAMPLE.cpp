//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Incompressible/Forces/FLUID_GRAVITY.h>
#include <Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <Dynamics/Drivers/PLS_EXAMPLE.h>
using namespace PhysBAM;
//#####################################################################
// PLS_EXAMPLE
//#####################################################################
template<class TV_input> PLS_EXAMPLE<TV_input>::
PLS_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    write_substeps_level(-1),write_output_files(true),output_directory("output"),restart(0),
    number_of_ghost_cells(3),cfl(.9),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac_grid),
    mpi_grid(0),projection(mac_grid),particle_levelset_evolution(mac_grid,collision_bodies_affecting_fluid,number_of_ghost_cells,false),incompressible(mac_grid,projection),boundary(0),
    phi_boundary(0),collision_bodies_affecting_fluid(mac_grid)
{
    incompressible.Set_Custom_Advection(advection_scalar);
    for(int i=0;i<TV::m;i++){domain_boundary(i)(0)=true;domain_boundary(i)(1)=true;}
    domain_boundary(1)(1)=false;
}
//#####################################################################
// ~PLS_EXAMPLE
//#####################################################################
template<class TV_input> PLS_EXAMPLE<TV_input>::
~PLS_EXAMPLE()
{
    if(mpi_grid){
        delete boundary;
        delete phi_boundary;}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV_input> void PLS_EXAMPLE<TV_input>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=LOG::sprintf("%d",frame);
    Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    if(mpi_grid) Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    Write_To_File(stream_type,output_directory+"/"+f+"/pressure",incompressible.projection.p);
    Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
    Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
    PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=particle_levelset_evolution.Particle_Levelset(0);
    Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    Write_To_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    Write_To_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    Write_To_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    Write_To_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV_input> void PLS_EXAMPLE<TV_input>::
Read_Output_Files(const int frame)
{
    std::string f=LOG::sprintf("%d",frame);
    PARTICLE_LEVELSET_UNIFORM<TV>& particle_levelset=particle_levelset_evolution.Particle_Levelset(0);
    Read_From_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    Read_From_File(stream_type,LOG::sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    Read_From_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
    std::string filename;
    filename=output_directory+"/"+f+"/pressure";
    if(File_Exists(filename)){LOG::cout<<"Reading pressure "<<filename<<std::endl;Read_From_File(stream_type,filename,incompressible.projection.p);}
    filename=output_directory+"/"+f+"/mac_velocities";
    if(File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;Read_From_File(stream_type,filename,face_velocities);}
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV_input> void PLS_EXAMPLE<TV_input>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=particle_levelset_evolution.Particle_Levelset(0).Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T min_collision_distance=particle_levelset_evolution.Particle_Levelset(0).min_collision_distance_factor*max_collision_distance;
    TV min_corner=mac_grid.domain.Minimum_Corner(),max_corner=mac_grid.domain.Maximum_Corner();
    for(int axis=0;axis<TV::m;axis++){
        if(domain_boundary[axis][0] && X_new[axis]<min_corner[axis]+max_collision_distance){
            T collision_distance=X[axis]-min_corner[axis];
            if(collision_distance>max_collision_distance)collision_distance=X_new[axis]-min_corner[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
            V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
        if(domain_boundary[axis][1] && X_new[axis]>max_corner[axis]-max_collision_distance){
            T collision_distance=max_corner[axis]-X[axis];
            if(collision_distance>max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
            V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
}
//#####################################################################
namespace PhysBAM{
template class PLS_EXAMPLE<VECTOR<float,2> >;
template class PLS_EXAMPLE<VECTOR<float,3> >;
template class PLS_EXAMPLE<VECTOR<double,2> >;
template class PLS_EXAMPLE<VECTOR<double,3> >;
}
