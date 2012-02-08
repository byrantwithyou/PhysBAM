//#####################################################################
// Copyright 2003, Ronald Fedkiw, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CYLINDER_IN_BOX
//##################################################################### 
//
//#####################################################################
// Koltakov - August 11, 2003
//#####################################################################
#ifndef __CYLINDER_IN_BOX__
#define __CYLINDER_IN_BOX__

#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_CYLINDER.h>
#include "../WATER_FREE_SURFACE_3D_EXAMPLE.h"
#include <Geometry/IMPLICIT_SURFACE_LIST.h>
#include <Geometry/TRIANGULATED_SURFACE_LIST.h>
namespace PhysBAM{

template<class T,class RW=T>
class CYLINDER_IN_BOX:public WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>
{
public:
    RENDERING_CYLINDER<T> cylinder; 

    CYLINDER_IN_BOX()
    {
        first_frame=0;
        last_frame=5*(3*24);
        restart=false;
        restart_frame=20;
        m=100;n=20;mn=100;
        xmin=(T)-1.25;xmax=(T)1.25;ymin=(T)-.25;ymax=(T).25;zmin=(T)-1.25;zmax=(T)1.25;
        domain_walls[0][0]=true;domain_walls[0][1]=true;domain_walls[1][0]=true;domain_walls[1][1]=false;domain_walls[2][0]=true;domain_walls[2][1]=true;
        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
        number_particles_per_cell=16;
        use_removed_positive_particles=false;use_removed_negative_particles=false;
        implicit_viscosity=true;use_central_differencing=false;second_order_pressure=false;variable_viscosity=true;
        write_levelset=true;write_velocity=false;write_removed_positive_particles=false;write_removed_negative_particles=false;write_particles=false;
        matlab_directory="Cylinder_in_box/matlab1";output_directory="Cylinder_in_box/output1";
        delete_fluid_inside_objects=false;

        // initialize the cylinder
        cylinder.cylinder.radius=(T).3;cylinder.cylinder.Set_Height((T)4.);
        cylinder.Update_Transform(MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(xmin+cylinder.cylinder.radius,(ymin+ymax)/2,(zmin+zmax)/2)));
    }

    ~CYLINDER_IN_BOX() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) for(int ij=0;ij<grid.mn;ij++) {
        if (j<(int)(grid.n*3.5/5)) phi(i,j,ij)=-grid.dx; else phi(i,j,ij)=grid.dx;
        if(cylinder.Inside(VECTOR<T,3>(grid.x(i),grid.y(j),grid.z(ij)))) phi(i,j,ij)=-grid.dx;}
}
//#####################################################################
// Function Adjust_Phi_And_Get_Velocities_For_Objects
//#####################################################################
// dt, u, v and w are not used
// doesn't use the obejct velocities!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void Adjust_Phi_And_Get_Velocities_For_Objects(const T dt,const T time,const ARRAY<T,VECTOR<int,3> >& u,const ARRAY<T,VECTOR<int,3> >& v,const ARRAY<T,VECTOR<int,3> >& w)
{
    T velocity = (T)2.0;
    cylinder.Set_Transform(MATRIX<T,4>::Identity_Matrix());
    cylinder.Update_Transform(MATRIX<T,4>::Translation_Matrix(VECTOR<T,3>(xmin+cylinder.cylinder.radius+velocity*time,(ymin+ymax)/2,(zmin+zmax)/2)));
    T epsilon=(T)1e-4*max(grid.dx,grid.dy,grid.dz);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++){
        VECTOR<T,3> location(grid.x(i),grid.y(j),grid.z(ij));
        if(cylinder.Inside(location)){
            phi(i,j,ij)=-grid.dx;
            VECTOR<T,3> surface=cylinder.Surface(location);
            VECTOR<T,3> normal=cylinder.Normal(surface);
            surface+=epsilon*normal;
            //phi(i,j,ij)=grid.interpolate(phi,surface); // should be using this to extrapolate phi inside the object ?????
            if(delete_fluid_inside_objects)    phi(i,j,ij)=(T)(-(cylinder.Surface(location)-location).Magnitude());
            psi_N(i,j,ij)=true;V_object(i,j,ij)=VECTOR<T,3>(velocity,0,0);}}
}
//#####################################################################
// Function Adjust_Particle_Velocity_For_Objects
//#####################################################################
// doesn't take into account the object velocity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void Adjust_Particle_Velocity_For_Objects(const GRID<TV>& grid,HEAVY_PARTICLES<T,VECTOR<T,3> >& particles,const T dt,const T time)
{
    for(int k=0;k<particles.array_size;k++) if(particles.active(k)){
        VECTOR<T,3> X=particles.X(k),V=particles.V(k),Xnew(X+dt*V);
        if(cylinder.Inside(Xnew)){
            VECTOR<T,3> N=cylinder.Normal(X);V-=VECTOR<T,3>::Dot_Product(V,N)*N;
            if(particle_restitution){T magnitude=V.Magnitude();if(magnitude != 0) V*=(magnitude+particle_restitution*(particles.V(k).Magnitude()-magnitude))/magnitude;}
            particles.V(k)=V;}}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(const GRID<TV>& grid,HEAVY_PARTICLES<T,VECTOR<T,3> >& particles,const T time)
{
    for(int k=0;k<particles.array_size;k++) if(particles.active(k)) 
        if(cylinder.Inside(particles.X(k))) particles.Delete_Particle(k);
}
//#####################################################################
};      
}
#endif    
