//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPLASH 
//##################################################################### 
#ifndef __SPLASH__
#define __SPLASH__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>

namespace PhysBAM{

template<class T,class RW=T>
class SPLASH:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_output_files;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;

    ARRAY<T,VECTOR<int,2> > phi_object;

    SPLASH()
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER)
    {
        first_frame=0;last_frame=300;
        frame_rate=24*3;
        restart=false;restart_frame=18;
        //fluids_parameters.grid.Initialize(281,221,0,(T)1.4,0,(T)1.1);
        fluids_parameters.grid.Initialize(141,111,0,(T)1.4,0,(T)1.1);
        //fluids_parameters.grid.Initialize(211,166,136,0,(T)1.4,0,(T)1.1,0,(T).9);
        //fluids_parameters.grid.Initialize(71,56,46,0,(T)1.4,0,(T)1.1,0,(T).9);
        //fluids_parameters.grid.Initialize(36,27,22,0,(T)1.4,0,(T)1.1,0,(T).9);
        fluids_parameters.domain_walls[1][1]=true;fluids_parameters.domain_walls[1][2]=true;fluids_parameters.domain_walls[2][1]=true;fluids_parameters.domain_walls[2][2]=false;
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.viscosity=(T)0;fluids_parameters.implicit_viscosity=false;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;
        fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=false;
        output_directory="Splash/output";
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;

        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
    }
    
    ~SPLASH() 
    {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<float>((data_directory+"/Rigid_Bodies_2D/circle").c_str(), (T).1, true, true, false);
    solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>((T)1.25,(T).55);
    solids_parameters.rigid_body_parameters.list(id)->is_kinematic=true;
    Add_To_Fluid_Simulation(*solids_parameters.rigid_body_parameters.list(id),true,false,100);
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    Update_Rigid_Bodies(time);
    phi_object.Resize(fluids_parameters.grid,3);
    for(int i=-2;i<=fluids_parameters.grid.m+3;i++) for(int j=-2;j<=fluids_parameters.grid.n+3;j++)
        phi_object(i,j)=solids_parameters.rigid_body_parameters.list(1)->Implicit_Curve_Extended_Value(fluids_parameters.grid.X(i,j));
}
//#####################################################################
// Function Update_Rigid_Bodies
//#####################################################################
void Update_Rigid_Bodies(const T time)
{
    INTERPOLATION_CURVE<T,VECTOR_2D<T> > motion_curve;
    motion_curve.Add_Control_Point(0,VECTOR_2D<T>((T)1.25,(T).55));
    motion_curve.Add_Control_Point(.2,VECTOR_2D<T>((T).8,(T).1)); // .03 was old
    motion_curve.Add_Control_Point(3,VECTOR_2D<T>((T).8,(T).1));
    for(int i=1;i<=solids_parameters.rigid_body_parameters.list.Number_Of_Elements();i++){
        solids_parameters.rigid_body_parameters.list(i)->position=motion_curve.Value(time);
        solids_parameters.rigid_body_parameters.list(i)->velocity=motion_curve.Derivative(time);}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++)
        fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.y(j)-(T).4;
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects()
{
    fluids_parameters.Extrapolate_Phi_Into_Object(fluids_parameters.particle_levelset_evolution.particle_levelset.levelset.phi,phi_object);
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    if(solids_parameters.rigid_body_parameters.list.Number_Of_Elements()<1) return;
    fluids_parameters.Adjust_Phi_With_Object(phi_object,*solids_parameters.rigid_body_parameters.list(1),time);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
bool Adjust_Particle_For_Objects(VECTOR_2D<T>& X,VECTOR_2D<T>& V,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T dt,const T time)
{          
    if(solids_parameters.rigid_body_parameters.list.Number_Of_Elements()<1) return true;
    fluids_parameters.Adjust_Particle_For_Object(*solids_parameters.rigid_body_parameters.list(1),X,V,particle_type,dt,time);
    return true;
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR_2D<T> >& particles,const typename PARTICLE_LEVELSET<T,VECTOR_2D<T> >::PARTICLE_TYPE particle_type,const T time)
{
    if(solids_parameters.rigid_body_parameters.list.Number_Of_Elements()<1) return;
    if(particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::NEGATIVE || particle_type == PARTICLE_LEVELSET<T,VECTOR_2D<T> >::REMOVED_NEGATIVE){
        for(int k=particles.array_collection->Size();k>=1;k--) if(solids_parameters.rigid_body_parameters.list(1)->Implicit_Curve_Lazy_Inside_Extended_Levelset(particles.X(k),-fluids_parameters.grid.dx)) particles.Delete_Particle(k);}
    else for(int k=particles.array_collection->Size();k>=1;k--) if(solids_parameters.rigid_body_parameters.list(1)->Implicit_Curve_Lazy_Inside_Extended_Levelset(particles.X(k))) particles.Delete_Particle(k);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(solids_parameters.rigid_body_parameters.list.Number_Of_Elements()<1) return;
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_object,*solids_parameters.rigid_body_parameters.list(1),3,true,time);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Limit_Dt(T& dt,const T time)
{
    VECTOR_2D<T> velocity=solids_parameters.rigid_body_parameters.list(1)->velocity;
    dt=min(dt,1/((T)fabs(velocity.x)/fluids_parameters.grid.dx+(T)fabs(velocity.y)/fluids_parameters.grid.dy));
}
//#####################################################################
};      
}
#endif    


