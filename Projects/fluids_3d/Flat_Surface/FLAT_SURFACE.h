//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLAT_SURFACE 
//##################################################################### 
#ifndef __FLAT_SURFACE__
#define __FLAT_SURFACE__

#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW=T>
class FLAT_SURFACE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_substeps;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::write_frame_title;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;

    ARRAY<T,VECTOR<int,3> > phi_object;
    LEVELSET_3D<GRID<TV> > levelset_object;
    ARRAY<VECTOR<T,3> ,VECTOR<int,3> > V_object;
    T water_level;
    bool use_objects;
    bool use_cylinder;
    bool moving_objects;
    int example;

    FLAT_SURFACE()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.WATER),levelset_object(fluids_parameters.grid,phi_object),water_level(.51234),
        use_objects(true),use_cylinder(true),moving_objects(false),example(0)
    {
        first_frame=0;last_frame=300;
        frame_rate=24;
        fluids_parameters.domain_walls[0][0]=true;fluids_parameters.domain_walls[0][1]=true;fluids_parameters.domain_walls[1][0]=true;
        fluids_parameters.domain_walls[1][1]=false;fluids_parameters.domain_walls[2][0]=true;fluids_parameters.domain_walls[2][1]=true;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=true;
        fluids_parameters.write_debug_data=true;
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.incompressible_iterations=100;

        restart=false;restart_frame=0;
        //fluids_parameters.use_strain=true;fluids_parameters.write_strain=true;
        //fluids_parameters.levelset_substeps=1;
        fluids_parameters.object_friction=1;
        fluids_parameters.adhesion_coefficient=0;
        fluids_parameters.use_old_velocities_for_boundary_conditions=true;
        fluids_parameters.second_order_pressure=true;

        fluids_parameters.bias_towards_negative_particles=true;fluids_parameters.number_particles_per_cell=16;
        fluids_parameters.use_removed_positive_particles=false;fluids_parameters.use_removed_negative_particles=true;
        
        //write_substeps=true;write_frame_title=true;

        example=3;
        output_directory=STRING_UTILITIES::string_sprintf("Flat_Surface/output%d",example);

        if(example==1){
            int r=7;
            fluids_parameters.grid.Initialize(5*r+1,5*r+1,5*r+1,0,1,0,1,0,1);
            fluids_parameters.object_friction=0;
            water_level=.5;}
        else if(example==2){
            //restart=true;restart_frame=19;
            int r=20;
            fluids_parameters.grid.Initialize(5*r+1,5*r+1,5*r+1,0,1,0,1,0,1);
            fluids_parameters.object_friction=0;
            solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",(T).15,true,true,false);}
        else if(example==3){
            //restart=true;restart_frame=19;
            int r=10;
            fluids_parameters.grid.Initialize(5*r+1,5*r+1,5*r+1,0,1,0,1,0,1);
            fluids_parameters.object_friction=1;
            //solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",(T).15,true,true,false);
            fluids_parameters.gravity_direction=VECTOR<T,3>(0,0,-1);
            fluids_parameters.domain_walls[1][1]=true;
            use_cylinder=false;}
        else{std::cerr<<"Dying inexplicably.\n";exit(1);}
    }

    ~FLAT_SURFACE()
    {}

///####################################################################
// Function Construct_Levelset_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    CYLINDER<T> cylinder(VECTOR<T,3>(.5,.1,.5),VECTOR<T,3>(.5,2,.5),.4);
    RIGID_BODY_LIST_3D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;

    if(!use_cylinder) cylinder.radius=30;

    if(example==2){
        moving_objects=true;
        rigid_body_particles.Rigid_Body(1).position=VECTOR<T,3>(.5,.3,.5);
        rigid_body_particles.Rigid_Body(1).velocity=VECTOR<T,3>(0,1,0);
        for(int r=0;r<rigid_body_list.rigid_bodies.m;r++)rigid_body_particles.Rigid_Body(r).position+=time*rigid_body_particles.Rigid_Body(r).velocity;}

    static bool initialized=false;if(!moving_objects && initialized) return;

    LOG::Time("Constructing levelsets");
    GRID<TV>& grid=fluids_parameters.grid;
    phi_object.Resize(grid,3,false,false);
    V_object.Resize(grid,3,false,false);
    for(int i=-2;i<=grid.m+3;i++)for(int j=-2;j<=grid.n+3;j++)for(int ij=-2;ij<=grid.mn+3;ij++){
        VECTOR<T,3> X=grid.X(i,j,ij); 
        T min_phi=-cylinder.Signed_Distance(grid.X(i,j,ij));int min_r=0;
        for(int r=0;r<rigid_body_list.rigid_bodies.m;r++){
            T phi=rigid_body_particles.Rigid_Body(r).Implicit_Surface_Value(X);
            if(min_phi>phi){min_phi=phi;min_r=r;}}
        phi_object(i,j,ij)=min_phi;
        V_object(i,j,ij)=min_r?rigid_body_particles.Rigid_Body(min_r).Pointwise_Object_Velocity(X):VECTOR<T,3>();}
    levelset_object.Compute_Cell_Minimum_And_Maximum();
    LOG::Stop_Time();
    initialized=true;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
        fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=grid.y(j)-water_level;
    if(example==3)
        for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)for(int ij=0;ij<grid.mn;ij++)
            fluids_parameters.particle_levelset_evolution.phi(i,j,ij)=grid.z(ij)-water_level;
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects(const T time)
{
    if(!use_objects) return;
    fluids_parameters.Extrapolate_Phi_Into_Object(phi_object);
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    if(!use_objects) return;
    fluids_parameters.Adjust_Phi_With_Object(phi_object,V_object,time);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const int index,VECTOR<T,3>& V,
    const typename PARTICLE_LEVELSET<T,VECTOR<T,3> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(!use_objects) return true;
    fluids_parameters.Adjust_Particle_For_Object(levelset_object,V_object,particles,index,V,particle_type,dt,time);
    return true;
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const typename PARTICLE_LEVELSET<T,VECTOR<T,3> >::PARTICLE_TYPE particle_type,const T time)
{
    if(!use_objects) return;
    T contour_value=particle_type == PARTICLE_LEVELSET<T,VECTOR<T,3> >::NEGATIVE || particle_type == PARTICLE_LEVELSET<T,VECTOR<T,3> >::REMOVED_NEGATIVE ? -fluids_parameters.grid.dx : 0;
    for(int k=particles.array_collection->Size();k>=1;k--)if(levelset_object.Lazy_Inside(particles.X(k),contour_value)) particles.Delete_Particle(k);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(!use_objects) return;
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_object,V_object,time);
}
//#####################################################################
// Function Adjust_Strain
//#####################################################################
void Adjust_Strain(ARRAY<SYMMETRIC_MATRIX<T,3> ,VECTOR<int,3> >& e_ghost,const T time)
{
    if(use_objects && fluids_parameters.use_strain)
        fluids_parameters.Adjust_Strain_For_Object(phi_object,e_ghost,time);
}
//#####################################################################
};      
}
#endif    


