//#####################################################################
// Copyright 2005, Ron Fedkiw, Eran Guendelman, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_MULTIPLE_UNIFORM.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
namespace PhysBAM{

//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
PARTICLE_LEVELSET_MULTIPLE_UNIFORM(GRID<TV>& grid_input,ARRAY<ARRAY<T,TV_INT>>& phis_input,const int number_of_ghost_cells_input)
    :levelset_multiple(grid_input,phis_input,true),number_of_ghost_cells(number_of_ghost_cells_input)
{
    Set_Collision_Distance_Factors(); // TODO: use this from a normal particle levelset
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
~PARTICLE_LEVELSET_MULTIPLE_UNIFORM()
{
    particle_levelsets.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize_Particle_Levelsets_And_Grid_Values
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Initialize_Particle_Levelsets_And_Grid_Values(GRID<TV>& grid,ARRAY<ARRAY<T,TV_INT>>& phis,
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_body_list_input,
    const int number_of_regions,const bool only_use_negative_particles)
{
    if(particle_levelsets.m!=number_of_regions){
        for(int i=0;i<particle_levelsets.m;i++)delete particle_levelsets(i);
        particle_levelsets.Resize(number_of_regions);levelset_multiple.levelsets.Resize(particle_levelsets.m);
        for(int i=0;i<particle_levelsets.m;i++){
            particle_levelsets(i)=new PARTICLE_LEVELSET_UNIFORM<TV>(grid,phis(i),collision_body_list_input,number_of_ghost_cells);
            particle_levelsets(i)->only_use_negative_particles=only_use_negative_particles;
            particle_levelsets(i)->Initialize_Particle_Levelset_Grid_Values();
            levelset_multiple.levelsets(i)=&particle_levelsets(i)->levelset;}}
    else for(int i=0;i<particle_levelsets.m;i++) particle_levelsets(i)->Initialize_Particle_Levelset_Grid_Values();
}
//#####################################################################
// Function Store_Unique_Particle_Id
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Store_Unique_Particle_Id()
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Store_Unique_Particle_Id();
}
//#####################################################################
// Function Set_Band_Width
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Set_Band_Width(const T number_of_cells) 
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Set_Band_Width(number_of_cells);
}
//#####################################################################
// Function Use_Removed_Negative_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Use_Removed_Negative_Particles(const bool use_removed_negative_particles_input)
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Use_Removed_Negative_Particles(use_removed_negative_particles_input);
}
//#####################################################################
// Function Use_Removed_Positive_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Use_Removed_Positive_Particles(const bool use_removed_positive_particles_input)
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Use_Removed_Positive_Particles(use_removed_positive_particles_input);
}
//#####################################################################
// Function Set_Collision_Distance_Factors
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Set_Collision_Distance_Factors(const T min_collision_distance_factor_input,const T max_collision_distance_factor_input)
{
    min_collision_distance_factor=min_collision_distance_factor_input;
    max_collision_distance_factor=max_collision_distance_factor_input;
    max_minus_min_collision_distance_factor_over_max_short=(max_collision_distance_factor-min_collision_distance_factor)/USHRT_MAX;
}
//#####################################################################
// Function Particle_Collision_Distance
//#####################################################################
template<class TV> auto PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Particle_Collision_Distance(const unsigned short quantized_collision_distance) -> T
{
    return (min_collision_distance_factor+(T)quantized_collision_distance*max_minus_min_collision_distance_factor_over_max_short)*levelset_multiple.grid.dX.Min();
}
//#####################################################################
// Function Seed_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Seed_Particles(const T time,const bool verbose)
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Seed_Particles(time,verbose);
}
//#####################################################################
// Function Adjust_Particle_Radii
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Adjust_Particle_Radii()
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Adjust_Particle_Radii();
}
//#####################################################################
// Function Modify_Levelset_Using_Escaped_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Modify_Levelset_Using_Escaped_Particles(ARRAY<T,FACE_INDEX<TV::m> >* face_velocities)
{
    for(int i=0;i<particle_levelsets.m;i++){
        ARRAY<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES*> other_positive_particles(particle_levelsets.m-1);
        int index=0;
        for(int j=0;j<particle_levelsets.m;j++)
            if(i!=j)
                other_positive_particles(index++)=&particle_levelsets(j)->negative_particles;
        particle_levelsets(i)->Modify_Levelset_Using_Escaped_Particles(face_velocities,&other_positive_particles);}
}
//#####################################################################
// Function Euler_Step_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Euler_Step_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const T dt,const T time,
    const bool use_second_order_for_nonremoved_particles,const bool update_particle_cells_after_euler_step,
    const bool verbose)
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Euler_Step_Particles(V,dt,time,use_second_order_for_nonremoved_particles,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Euler_Step_Removed_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step,const bool verbose)
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Euler_Step_Removed_Particles(dt,time,update_particle_cells_after_euler_step,verbose);
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class TV> int PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Reseed_Particles(const T time,ARRAY<bool,TV_INT>* cell_centered_mask)
{
    int new_particles=0;
    for(int i=0;i<particle_levelsets.m;i++)
        new_particles+=particle_levelsets(i)->Reseed_Particles(time,cell_centered_mask);
    return new_particles;
}
//#####################################################################
// Function Delete_Particles_Outside_Grid
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Delete_Particles_Outside_Grid()
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Delete_Particles_Outside_Grid();
}
//#####################################################################
// Function Identify_And_Remove_Escaped_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Identify_And_Remove_Escaped_Particles(const ARRAY<T,FACE_INDEX<TV::m> >& V,const T radius_fraction,const bool verbose)
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Identify_And_Remove_Escaped_Particles(V,radius_fraction,verbose);
}
//#####################################################################
// Function Reincorporate_Removed_Particles
//#####################################################################
template<class TV> void PARTICLE_LEVELSET_MULTIPLE_UNIFORM<TV>::
Reincorporate_Removed_Particles(const T radius_fraction,const T mass_scaling,ARRAY<T,FACE_INDEX<TV::m> >* V,const bool conserve_momentum_for_removed_negative_particles)
{
    for(int i=0;i<particle_levelsets.m;i++)
        particle_levelsets(i)->Reincorporate_Removed_Particles(radius_fraction,mass_scaling,V,
            conserve_momentum_for_removed_negative_particles);
}
template class PARTICLE_LEVELSET_MULTIPLE_UNIFORM<VECTOR<double,1> >;
template class PARTICLE_LEVELSET_MULTIPLE_UNIFORM<VECTOR<double,2> >;
template class PARTICLE_LEVELSET_MULTIPLE_UNIFORM<VECTOR<double,3> >;
template class PARTICLE_LEVELSET_MULTIPLE_UNIFORM<VECTOR<float,1> >;
template class PARTICLE_LEVELSET_MULTIPLE_UNIFORM<VECTOR<float,2> >;
template class PARTICLE_LEVELSET_MULTIPLE_UNIFORM<VECTOR<float,3> >;
}
