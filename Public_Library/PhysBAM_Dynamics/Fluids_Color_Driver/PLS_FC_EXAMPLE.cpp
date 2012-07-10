//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Fluids_Color_Driver/PLS_FC_EXAMPLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PLS_FC_EXAMPLE<TV>::
PLS_FC_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),last_frame(100),
    write_substeps_level(-1),write_output_files(true),output_directory("output"),restart(0),
    number_of_ghost_cells(3),dt(0),time_steps_per_frame(1),use_preconditioner(true),max_iter(1000),
    dump_matrix(false),wrap(true),grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    particle_levelset_evolution(grid,number_of_ghost_cells),boundary(0),
    levelset_color(grid,*new ARRAY<T,TV_INT>,*new ARRAY<int,TV_INT>),collision_bodies_affecting_fluid(grid)
{
    for(int i=0;i<TV::dimension;i++){domain_boundary(i)(0)=true;domain_boundary(i)(1)=true;}
    domain_boundary(1)(1)=false;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PLS_FC_EXAMPLE<TV>::
~PLS_FC_EXAMPLE()
{
    delete &levelset_color.color;
    delete &levelset_color.phi;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void PLS_FC_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
//    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
//    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class TV> void PLS_FC_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Read_From_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
    std::string filename;
//    filename=output_directory+"/"+f+"/pressure";
//    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading pressure "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);}
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){LOG::cout<<"Reading mac_velocities "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV> void PLS_FC_EXAMPLE<TV>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T min_collision_distance=particle_levelset_evolution.particle_levelset.min_collision_distance_factor*max_collision_distance;
    TV min_corner=grid.domain.Minimum_Corner(),max_corner=grid.domain.Maximum_Corner();
    for(int axis=0;axis<GRID<TV>::dimension;axis++){
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
template class PLS_FC_EXAMPLE<VECTOR<float,2> >;
template class PLS_FC_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_FC_EXAMPLE<VECTOR<double,2> >;
template class PLS_FC_EXAMPLE<VECTOR<double,3> >;
#endif
