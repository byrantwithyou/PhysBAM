//#####################################################################
// Copyright 2003, Doug Enright, Ronald Fedkiw, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLOW_OVER_SPHERE 
//##################################################################### 
#ifndef __FLOW_OVER_SPHERE__
#define __FLOW_OVER_SPHERE__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_CYLINDER.h>
#include "../WATER_FREE_SURFACE_3D_EXAMPLE.h"
#include <Geometry/IMPLICIT_SURFACE_LIST.h>
#include <Geometry/TRIANGULATED_SURFACE_LIST.h>
namespace PhysBAM{

template<class T,class RW=T>
class FLOW_OVER_SPHERE:public WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>
{
public:
    RIGID_BODY_LIST_3D<T> rigid_body_list;
    ARRAY<RIGID_BODY<TV> *>& rigid_bodies;
    TRIANGULATED_SURFACE_LIST<T> triangulated_surface_list;
    IMPLICIT_SURFACE_LIST<T> implicit_surface_list;
    T cylinder_radius,cylinder_height;
    ARRAY<T,VECTOR<int,3> > phi_object;
    RANDOM_NUMBERS random;

    FLOW_OVER_SPHERE():rigid_bodies(rigid_body_list.rigid_bodies),cylinder_radius((T).13),cylinder_height((T).15)
    {
        first_frame=0;last_frame=200;
        frame_rate=200;
        restart=false;restart_frame=21;
        m=80;n=80;mn=80;
        lagrangian_steps=(T)4.9;
        xmin=(T)-.5;xmax=(T).5;ymin=(T)-.4;ymax=(T).6;zmin=(T)-.5;zmax=(T).5;
        domain_walls[1][1]=true;domain_walls[1][2]=true;domain_walls[2][1]=true;domain_walls[2][2]=false;domain_walls[3][1]=true;domain_walls[3][2]=true;
        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
        bias_towards_negative_particles=true;number_particles_per_cell=16;particle_restitution=0;use_removed_positive_particles=false;use_removed_negative_particles=false;
        viscosity=0;//(T)(.001137*2e5);
        implicit_viscosity=true;variable_viscosity=true;use_central_differencing=false;second_order_pressure=false;
        write_matlab_file=false;write_levelset=true;write_velocity=true;write_particles=true;write_removed_positive_particles=false;write_removed_negative_particles=true;write_velocity=true;
        matlab_directory="Flow_Over_Sphere/matlab";output_directory="Flow_Over_Sphere/output1";
        delete_fluid_inside_objects=true;
        enforce_divergence_free_extrapolation=false;

        // initialize rigid bodies
        rigid_body_list.read_rgd_file=false;    // Our blob doesn't have an rgd file
        rigid_body_list.template Add_Rigid_Body<float>("Flow_Over_Sphere/sphere", (T).175);
        rigid_bodies(1)->implicit_surface->Compute_Cell_Minimum_And_Maximum(); // !!!!!!!!!!!!!!!!!!!!!!
    }

    ~FLOW_OVER_SPHERE() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) phi(i,j,ij)=grid.dx;
}
//#####################################################################
// Function Specify_Variable_Viscosity_Values
//#####################################################################
void Specify_Variable_Viscosity_Values(ARRAY<T,VECTOR<int,3> >& viscosity,const T time)
{
    for(int i=0;i<=m+1;i++) for(int j=0;j<=n+1;j++) for(int ij=0;ij<=mn+1;ij++)
        if(grid.x(i) > 0) viscosity(i,j,ij)=(T)1e6*.001137;
        else viscosity(i,j,ij)=(T)(1e6*.001137*1e-6);
}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time)
{
    RENDERING_CYLINDER<T> cylinder; 
    cylinder.cylinder.radius=cylinder_radius;cylinder.cylinder.Set_Height(cylinder_height);
    cylinder.Update_Transform(MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>((grid.xmin+grid.xmax)/2,grid.ymax,(grid.zmin+grid.zmax)/2)));
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++) if(cylinder.Lazy_Inside(VECTOR<T,3>(grid.x(i),grid.y(j),grid.z(ij)))) phi(i,j,ij)=-grid.dx;
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time)
{
    RENDERING_CYLINDER<T> cylinder; 
    cylinder.cylinder.radius=cylinder_radius;cylinder.cylinder.Set_Height(cylinder_height);
    cylinder.Update_Transform(MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>((grid.xmin+grid.xmax)/2,grid.ymax,(grid.zmin+grid.zmax)/2)));
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++) if(cylinder.Lazy_Inside(VECTOR<T,3>(grid.x(i),grid.y(j),grid.z(ij)))){
        psi_N(i,j,ij)=true;V_source(i,j,ij)=VECTOR<T,3>(0,-4,0);}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,3> >*& cell_centered_mask,const T time)
{
    if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool,VECTOR<int,3> >(1,grid.m,1,grid.n,1,grid.mn);

    T padding=3*grid.dx;
    RENDERING_CYLINDER<T> cylinder_mask;cylinder_mask.cylinder.radius=(T)(cylinder_radius+padding);
    cylinder_mask.cylinder.Set_Height((T)(cylinder_height+2*padding));
    cylinder_mask.Update_Transform(MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>((grid.xmin+grid.xmax)/2,grid.ymax+padding,(grid.zmin+grid.zmax)/2)));
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) if(cylinder_mask.Lazy_Inside(VECTOR<T,3>(grid.x(i),grid.y(j),grid.z(ij)))) 
        (*cell_centered_mask)(i,j,ij)=true;
}
//#####################################################################
// Function Construct_Levelsets_For_Objects
//#####################################################################
void Construct_Levelsets_For_Objects(const T time)
{
    static bool initialized=false;
    if(!initialized){ // only need to this once since the object never changes
        phi_object.Resize(-2,m+3,-2,n+3,-2,mn+3);
        for(int i=-2;i<=m+3;i++) for(int j=-2;j<=n+3;j++) for(int k=-2;k<=mn+3;k++) phi_object(i,j,k)=-rigid_bodies(1)->Implicit_Surface_Value(grid.X(i,j,k));
        initialized=true;}
}
//#####################################################################
// Function Adjust_Phi_with_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    T one_over_two_dx=1/(2*grid.dx),one_over_two_dy=1/(2*grid.dy),one_over_two_dz=1/(2*grid.dz);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int k=0;k<mn;k++) if(phi(i,j,k) < 0 && phi_object(i,j,k) > 0){
        VECTOR<T,3> V_fluid=VECTOR<T,3>((T).5*(u(i,j,k)+u(i+1,j,k)),(T).5*(v(i,j,k)+v(i,j+1,k)),(T).5*(w(i,j,k)+w(i,j,k+1)));
        //VECTOR<T,3> V_object=rigid_bodies(1)->Pointwise_Object_Velocity(grid.X(i,j,k)); // velocity object should be spatially varying
        VECTOR<T,3> V_relative=V_fluid;// should be V_fluid-V_object, but ibject is sitting still
        VECTOR<T,3> normal=-VECTOR<T,3>((phi_object(i+1,j,k)-phi_object(i-1,j,k))*one_over_two_dx,(phi_object(i,j+1,k)-phi_object(i,j-1,k))*one_over_two_dy,
                                                                        (phi_object(i,j,k+1)-phi_object(i,j,k-1))*one_over_two_dz);
        T denominator=normal.Magnitude();if(denominator > 1e-8) normal/=denominator;else normal=VECTOR<T,3>(1,0,0);
        T VN=VECTOR<T,3>::Dot_Product(V_relative,normal),magnitude=V_relative.Magnitude();
        if(VN > .1*magnitude) phi(i,j,k)=phi_object(i,j,k);}
}
//#####################################################################
// Function Extrapolate_Levelset_Into_Objects
//#####################################################################
void Extrapolate_Levelset_Into_Objects(const T time)
{
    WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>::Extrapolate_Levelset_Into_Object(grid,phi_object,phi);
}
//#####################################################################
// Function Get_Velocities_For_Objects
//#####################################################################
void Get_Velocities_For_Objects(const T dt,const T time)
{
    fluids_parameters.Extrapolate_Velocity_Into_Object(grid,phi_object,u,v,w,*rigid_bodies(1),psi_N,V_object,time);
}
//#####################################################################
// Function Adjust_Particle_For_Objects
//#####################################################################
// doesn't take into account the object velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void Adjust_Particle_For_Objects(VECTOR<T,3>& X,VECTOR<T,3>& V,T& radius,PARTICLE_LEVELSET_3D<GRID<TV> >& particle_levelset,
    const typename WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>::PARTICLE_TYPE particle_type,const T dt,const T time)
{           
    if(particle_type == POSITIVE || particle_type == REMOVED_POSITIVE) return;
    
    VECTOR<T,3> X_new=X,V_new=V;T phi;
    T collision_distance=grid.dx*(radius-particle_levelset.minimum_particle_radius)/(particle_levelset.maximum_particle_radius-particle_levelset.minimum_particle_radius); // .05 - .25 -> 0-1
    if(rigid_bodies(1)->Implicit_Surface_Lazy_Inside_And_Value(X_new,phi,collision_distance)){
        if(radius == maximum_particle_radius){
            radius=random.Get_Uniform_Number(particle_levelset.minimum_particle_radius,particle_levelset.maximum_particle_radius);
            collision_distance=grid.dx*(radius-particle_levelset.minimum_particle_radius)/(particle_levelset.maximum_particle_radius-particle_levelset.minimum_particle_radius);} 
         VECTOR<T,3> N=rigid_bodies(1)->Implicit_Surface_Normal(X_new);
         X+=(-phi+collision_distance)*N; // could work to make this more accurate 
         T VN=VECTOR<T,3>::Dot_Product(V_new,N);
        if(VN < 0){
            V_new-=VN*N;
            if(particle_restitution && particle_type != REMOVED_NEGATIVE){
                T magnitude=V_new.Magnitude();if(magnitude != 0) V_new*=(magnitude+particle_restitution*(V.Magnitude()-magnitude))/magnitude;}
            V=V_new;}}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(HEAVY_PARTICLES<T,VECTOR<T,3> >& particles,const typename WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>::PARTICLE_TYPE particle_type,const T time)
{
    if(particle_type == NEGATIVE || particle_type == REMOVED_NEGATIVE){
        for(int k=0;k<particles.array_collection->Size();k++) if(rigid_bodies(1)->Implicit_Surface_Lazy_Inside_Extended_Levelset(particles.X(k),-grid_3d.dx)) particles.array_collection->Add_To_Deletion_List(k);}
   else for(int k=0;k<particles.array_collection->Size();k++) if(rigid_bodies(1)->Implicit_Surface_Lazy_Inside_Extended_Levelset(particles.X(k))) particles.array_collection->Add_To_Deletion_List(k);
   particles.array_collection->Delete_Elements_On_Deletion_List(false,true); // already sorted
}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
void Read_Output_Files_Fluids(const int frame)
{      
    WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>::Read_Output_Files_Fluids(frame);
    rigid_body_list.template Read_Frame<RW>(output_directory,frame);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const
{
    WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>::Write_Output_Files(frame);
    if(frame==first_frame){rigid_body_list.template Write_Initial_Data<RW>(output_directory);}
    rigid_body_list.template Write_Frame<RW>(output_directory,frame);
}  
//#####################################################################
};      
}
#endif    


