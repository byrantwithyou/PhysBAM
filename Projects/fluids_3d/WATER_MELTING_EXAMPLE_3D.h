//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_MELTING_EXAMPLE_3D
//#####################################################################
#ifndef __WATER_MELTING_EXAMPLE_3D__
#define __WATER_MELTING_EXAMPLE_3D__

#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Dynamics/Particles/SPH_PARTICLES.h>
#include <Level_Sets/LEVELSET_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Level_Sets/LEVELSET_RED_GREEN.h>
#include <Solids_And_Fluids/MELTING_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW>
class WATER_MELTING_EXAMPLE_3D:public MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_collision_body_list_in_particle_levelset;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_collision_aware_velocity_extrapolation;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_phi;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::thin_shells_semi_lagrangian_velocity;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_collision_aware_signed_distance;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::clamp_phi_with_collision_bodies;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::use_melting;
    using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::melting_parameters;using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::Update_Solids_Topology_For_Melting;

    ARRAY<ARRAY<VECTOR_3D<T> > > objects_particles;
    ARRAY<ARRAY<unsigned short> > objects_particles_collision_distance;
    T current_phi;

    ARRAY<T,VECTOR<int,3> > phi_objects; // level set for the objects
    ARRAY<VECTOR_3D<T> ,VECTOR<int,3> > V_objects; // level set for the objects
    bool initialized;

    WATER_MELTING_EXAMPLE_3D()
        :MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >(fluids_parameters.WATER,false),initialized(false)
    {
        // make sure we turn off particle thin shell handling in PLS (too slow with it)
        use_collision_body_list_in_particle_levelset=false;
        fluids_parameters.particle_levelset_evolution.particle_levelset.reincorporate_removed_particles_everywhere=true;

    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::Initialize_Bodies();       
    if(!initialized){initialized=true;
        objects_particles.Resize(melting_parameters.levelsets.m);
        objects_particles_collision_distance.Resize(melting_parameters.levelsets.m);
        if(!restart)for(int object=1;object<=melting_parameters.levelsets.m;object++){
            Seed_Particles_Inside_Object(objects_particles(object),objects_particles_collision_distance(object),object);}}
}
//#####################################################################
// Function Seed_Particles_Inside_Object
//#####################################################################
void Seed_Particles_Inside_Object(ARRAY<VECTOR_3D<T> >& object_particles,ARRAY<unsigned short>& object_particles_collision_distance,int object)
{
}
//#####################################################################
// Function Find_Escaped_Particles
//#####################################################################
void Find_Escaped_Particles(ARRAY<VECTOR_3D<T> >& object_particles,ARRAY<unsigned short>& object_particles_collision_distance,int object,ARRAY<VECTOR_3D<T> >& escaped_particles,ARRAY<unsigned short>& escaped_particles_collision_distance)
{
    melting_parameters.levelsets(object)->levelset.Lazy_Update_Overlay_Levelset();
    LEVELSET_RED_GREEN_3D<T>& levelset=melting_parameters.levelsets(object)->levelset;
    for(int i=object_particles.m;i>=1;i--)if(levelset.Phi(object_particles(i))>((T)object_particles_collision_distance(i)/USHRT_MAX)*fluids_parameters.grid.max_dx_dy_dz){
        escaped_particles.Append(object_particles(i));escaped_particles_collision_distance.Append(object_particles_collision_distance(i));object_particles.Remove_Index_Lazy(i);}
    std::cout << "Identified " << escaped_particles.m << " as escaped particles" << std::endl;
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
virtual void Extrapolate_Phi_Into_Objects(const T time)
{
    fluids_parameters.Extrapolate_Phi_Into_Object(phi_objects);
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    fluids_parameters.Adjust_Phi_With_Object(phi_objects,V_objects,time);
}
//#####################################################################
// Function Melting_Substep
//#####################################################################
virtual void Melting_Substep(const T dt,const T time)
{
    Update_Solids_Topology_For_Melting(dt,time);
}
//#####################################################################
// Function Modify_Fluid_For_Melting
//#####################################################################
void Modify_Fluid_For_Melting(const T dt,const T time)
{
    GRID<TV>& grid=fluids_parameters.grid;
/*
    ARRAY<bool,VECTOR<int,3> > nodes_inside_objects_new(grid,0);
    Mark_Nodes_Inside_Objects(nodes_inside_objects_new);
    for(int i=1;i<=grid.m;i++)for(int j=1;j<=grid.n;j++)for(int ij=1;ij<=grid.mn;ij++) if((*nodes_inside_objects)(i,j,ij)&&!nodes_inside_objects_new(i,j,ij)) 
        fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=-grid.max_dx_dy_dz;
    ARRAY<bool,VECTOR<int,3> >::Exchange_Arrays(nodes_inside_objects_new,*nodes_inside_objects);
*/
    // remove all the air particles that are in/next to the object
    PARTICLE_LEVELSET_3D<GRID<TV> >& particle_levelset=fluids_parameters.particle_levelset_evolution.particle_levelset;
    for(int i=1;i<=grid.number_of_cells_x;i++)for(int j=1;j<=grid.number_of_cells_y;j++)for(int ij=1;ij<=grid.number_of_cells_z;ij++){
        if(particle_levelset.positive_particles(i,j,ij))
            for(int ii=i;ii<=i+1;ii++)for(int jj=j;jj<=j+1;jj++)for(int iijj=ij;iijj<=ij+1;iijj++)
                if(phi_objects(ii,jj,iijj)<0){
                    delete particle_levelset.positive_particles(i,j,ij);particle_levelset.positive_particles(i,j,ij)=0;
                    goto HERE;}
        HERE:;
    }
    // add new water particles from the melting of the object
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        int index=melting_parameters.body_index(object);
        if(!index)continue;
        if(melting_parameters.body_type(object)==melting_parameters.DEFORMABLE){
            LEVELSET_TETRAHEDRALIZED_VOLUME<T>& levelset=*melting_parameters.levelsets(object);
            RED_GREEN_GRID_3D<T>& grid=levelset.grid;
            EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume=levelset.embedded_tetrahedralized_volume;
            PARTICLES<T,VECTOR_3D<T> >& particles=embedded_tetrahedralized_volume.particles;
            ARRAY<VECTOR_3D<T> > cell_based_X(grid.number_of_nodes);
            ARRAY<VECTOR_3D<T> > cell_based_V(grid.number_of_nodes);
            for(int i=1;i<=grid.number_of_nodes;i++)if(levelset.node_to_particle_mapping(i)){
                cell_based_X(i)=particles.X.array(levelset.node_to_particle_mapping(i));
                cell_based_V(i)=particles.V.array(levelset.node_to_particle_mapping(i));}
            ARRAY<VECTOR_3D<T> > escaped_particles;ARRAY<unsigned short> escaped_particles_collision_distance;
            Find_Escaped_Particles(objects_particles(object),objects_particles_collision_distance(object),object,escaped_particles,escaped_particles_collision_distance);
            for(int i=1;i<=escaped_particles.m;i++)
                fluids_parameters.particle_levelset_evolution.particle_levelset.Add_Negative_Particle(grid.Interpolate_Nodes(cell_based_X,escaped_particles(i)),grid.Interpolate_Nodes(cell_based_V,escaped_particles(i)),escaped_particles_collision_distance(i));}
        else if(melting_parameters.body_type(object)==melting_parameters.RIGID){
            RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(index);
            ARRAY<VECTOR_3D<T> > escaped_particles;ARRAY<unsigned short> escaped_particles_collision_distance;
            Find_Escaped_Particles(objects_particles(object),objects_particles_collision_distance(object),object,escaped_particles,escaped_particles_collision_distance);
            FRAME<T> frame=rigid_body.Frame()*melting_parameters.rigid_body_grid_frames(object).Inverse();
            for(int i=1;i<=escaped_particles.m;i++)
                fluids_parameters.particle_levelset_evolution.particle_levelset.Add_Negative_Particle(frame*escaped_particles(i),VECTOR_3D<T>(0,0,0),escaped_particles_collision_distance(i));}}

    Update_Solids_Topology_For_Melting_Part2(dt,time);
    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        ARRAY<VECTOR_3D<T> > escaped_particles;ARRAY<unsigned short> escaped_particles_collision_distance;
        Find_Escaped_Particles(objects_particles(object),objects_particles_collision_distance(object),object,escaped_particles,escaped_particles_collision_distance);}

    int number_of_particles_deleted=fluids_parameters.particle_levelset_evolution.particle_levelset.Reseed_Delete_Particles(fluids_parameters.particle_levelset_evolution.particle_levelset.negative_particles,-1);
    std::cout << "deleted " << number_of_particles_deleted << " negative_particles by calling reseed_delete_particles" << std::endl;
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Construct_Levelsets_For_Objects(const T time,ARRAY<T,VECTOR<int,3> >* static_collision_bodies_phi)
{
    Construct_Levelset_For_Objects(fluids_parameters.grid,phi_objects,V_objects,static_collision_bodies_phi);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
virtual void Get_Object_Velocities(const T dt,const T time)
{
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_objects,V_objects,3,true,time);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
virtual bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_3D<T> >& particles,const int index,VECTOR_3D<T>& V,const typename PARTICLE_LEVELSET<T,VECTOR_3D<T> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{
    LEVELSET_3D<GRID<TV> > levelset(fluids_parameters.grid,phi_objects);
    fluids_parameters.Adjust_Particle_For_Object(levelset,V_objects,particles,index,V,particle_type,dt,time);
    return true;
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
virtual void Read_Output_Files_Solids(const int frame)
{
    MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::Read_Output_Files_Solids(frame);
//    std::string prefix=output_directory+"/";
//    FILE_UTILITIES::Read_From_File<RW>(prefix+STRING_UTILITIES::string_sprintf("object_particles.%d",frame),objects_particles);

    for(int object=1;object<=melting_parameters.levelsets.m;object++){
        Seed_Particles_Inside_Object(objects_particles(object),objects_particles_collision_distance(object),object);}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
virtual void Write_Output_Files(const int frame) const
{
    MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::Write_Output_Files(frame);
//    std::string prefix=output_directory+"/";
//    FILE_UTILITIES::Write_To_File<RW>(prefix+STRING_UTILITIES::string_sprintf("object_particles.%d",frame),objects_particles);
}
//#####################################################################
};
}
#endif
