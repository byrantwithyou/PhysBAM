//#####################################################################
// Copyright 2003-2004, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOW_OVER_BLOB 
//##################################################################### 
#ifndef __FLOW_OVER_BLOB__
#define __FLOW_OVER_BLOB__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_CYLINDER.h>
#include "../WATER_FREE_SURFACE_3D_EXAMPLE.h"
namespace PhysBAM{

template<class T,class RW=T>
class FLOW_OVER_BLOB:public WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>
{
public:
    RIGID_BODY_LIST_3D<T> rigid_body_list;
    ARRAY<RIGID_BODY<TV> *>& rigid_bodies;
    ARRAY<T,VECTOR<int,3> > phi_object;
    RANDOM_NUMBERS random;
    bool use_cylinder_source;
    RENDERING_CYLINDER<T> cylinder;
    T source_speed;

    FLOW_OVER_BLOB()
        :rigid_bodies(rigid_body_list.rigid_bodies)
    {
        T width=(T)61.5,height=(T)36,depth=(T)43.5,xoff=(T)0,yoff=(T)7.42238,zoff=(T).915671;
        grid.Initialize(124,73,88,(T)(-.5*width+xoff),(T)(.5*width+xoff),(T)(-height*.5+yoff),(T)(height*.5+yoff),(T)(-depth*.5+zoff),(T)(depth*.5+zoff));
        first_frame=0;last_frame=5*24*4;
        frame_rate=24;
        restart=true;restart_frame=19;
        domain_walls[0][0]=true;domain_walls[0][1]=true;domain_walls[1][0]=true;domain_walls[1][1]=false;domain_walls[2][0]=true;domain_walls[2][1]=true;
        reseeding_frame_rate=1;
        bias_towards_negative_particles=true;number_particles_per_cell=16;use_removed_positive_particles=false;use_removed_negative_particles=true;
        viscosity=0;//(T)(.001137*2e5);
        implicit_viscosity=false;variable_viscosity=false;second_order_pressure=false;
        write_levelset=true;write_velocity=true;write_particles=true;write_removed_positive_particles=false;write_removed_negative_particles=true;
        write_debug_data=true;
        store_particle_ids=true;
        output_directory="Flow_Over_Blob/output";
        delete_fluid_inside_objects=true;
        enforce_divergence_free_extrapolation=false;

        // initialize rigid bodies
        rigid_body_list.template Add_Rigid_Body<float>("Flow_Over_Blob/blob",1,true,true,false);
        rigid_bodies(1)->implicit_surface->Compute_Cell_Minimum_And_Maximum(); // !!!!!!!!!!!!!!!!!!!!!!

        use_cylinder_source=true;
        cylinder.cylinder.radius=2;
        cylinder.cylinder.Set_Height(2);
        cylinder.Update_Transform(MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(0,22,2.5)));
        source_speed=-5;
    }

    ~FLOW_OVER_BLOB() 
    {}

//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    if(rigid_bodies.m<1) return;
    static bool initialized=false;
    if(!initialized){ // only need to this once since the object never changes
    phi_object.Resize(grid,3);
        for(int i=-2;i<=grid.m+3;i++) for(int j=-2;j<=grid.n+3;j++) for(int k=-2;k<=grid.mn+3;k++) phi_object(i,j,k)=rigid_bodies(1)->Implicit_Surface_Value(grid.X(i,j,k));
        initialized=true;}
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) phi(i,j,ij)=grid.dx;
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects()
{
    if(rigid_bodies.m<1) return;
    fluids_parameters.Extrapolate_Phi_Into_Object(phi,phi_object);
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    if(rigid_bodies.m<1) return;
    fluids_parameters.Adjust_Phi_With_Object(phi_object,*rigid_bodies(1),phi,time);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    if(!use_cylinder_source) return;
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++){
        VECTOR<T,3> X(grid.X(i,j,ij));
        if(cylinder.Lazy_Inside(X)) phi(i,j,ij)=min(phi(i,j,ij),cylinder.Signed_Distance(X));}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time)
{
    if(!use_cylinder_source) return;
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,3> >(grid);

    T padding=3*grid.dx;
    RENDERING_CYLINDER<T> cylinder_mask=cylinder;cylinder_mask.cylinder.radius+=padding;cylinder_mask.cylinder.Set_Height(cylinder.cylinder.height+2*padding);
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) if(cylinder_mask.Lazy_Inside(grid.X(i,j,ij))) (*cell_centered_mask)(i,j,ij)=true;
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
void Adjust_Particle_For_Objects(VECTOR<T,3>& X,VECTOR<T,3>& V,const typename PARTICLE_LEVELSET<T>::PARTICLE_TYPE particle_type,const T dt,const T time)
{          
    if(rigid_bodies.m<1) return;
    fluids_parameters.Adjust_Particle_For_Object(*rigid_bodies(1),X,V,particle_type,dt,time);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<T,VECTOR<T,3> >& particles,const typename PARTICLE_LEVELSET<T>::PARTICLE_TYPE particle_type,const T time)
{
    if(rigid_bodies.m<1) return;
    if(particle_type == NEGATIVE || particle_type == REMOVED_NEGATIVE){
        for(int k=particles.Size();k>=1;k--) if(rigid_bodies(1)->Implicit_Surface_Lazy_Inside(particles.X(k),-grid.dx)) particles.Delete_Particle(k);}
   else for(int k=particles.Size();k>=1;k--) if(rigid_bodies(1)->Implicit_Surface_Lazy_Inside(particles.X(k))) particles.Delete_Particle(k);
}
//#####################################################################
// Function Get_Variable_Viscosity
//#####################################################################
void Get_Variable_Viscosity(ARRAY<T,VECTOR<int,3> >& variable_viscosity,const T time)
{
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) 
        if(grid.x(i) > .5) variable_viscosity(i,j,ij)=1e6*.001137;
        else variable_viscosity(i,j,ij)=(T)(1e6*.001137*1e-6);
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    if(!use_cylinder_source) return;
    VECTOR<T,3> V_source(0,source_speed,0);
    for(int i=0;i<u_grid.m;i++) for(int j=0;j<u_grid.n;j++) for(int ij=0;ij<u_grid.mn;ij++) if(cylinder.Lazy_Inside(u_grid.X(i,j,ij))){psi_N_u(i,j,ij)=true;u_fixed(i,j,ij)=V_source.x;}
    for(int i=0;i<v_grid.m;i++) for(int j=0;j<v_grid.n;j++) for(int ij=0;ij<v_grid.mn;ij++) if(cylinder.Lazy_Inside(v_grid.X(i,j,ij))){psi_N_v(i,j,ij)=true;v_fixed(i,j,ij)=V_source.y;}
    for(int i=0;i<w_grid.m;i++) for(int j=0;j<w_grid.n;j++) for(int ij=0;ij<w_grid.mn;ij++) if(cylinder.Lazy_Inside(w_grid.X(i,j,ij))){psi_N_w(i,j,ij)=true;w_fixed(i,j,ij)=V_source.z;}
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time)
{
    if(rigid_bodies.m<1) return;
    fluids_parameters.Extrapolate_Velocity_Into_Object(phi_object,*rigid_bodies(1),psi_N_u,psi_N_v,psi_N_w,u_fixed,v_fixed,w_fixed,phi,time);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
void Read_Output_Files_Fluids(const int frame)
{      
    WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>::Read_Output_Files(frame);
    rigid_body_list.template Read_Dynamic_Variables<RW>(output_directory,frame);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>::Write_Output_Files(frame);
    if(frame==first_frame){rigid_body_list.template Write_Static_Variables<RW>(output_directory);}
    rigid_body_list.template Write_Dynamic_Variables<RW>(output_directory,frame);
}  
//#####################################################################
};    
}
#endif  
