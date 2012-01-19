//#####################################################################
// Copyright 2002, Sergey Koltakov, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_OVER_CIRCLE
//##################################################################### 
#ifndef __WATER_OVER_CIRCLE__
#define __WATER_OVER_CIRCLE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>

#include "../WATER_FREE_SURFACE_2D_EXAMPLE.h"

namespace PhysBAM{
template <class T,class RW=T>
class WATER_OVER_CIRCLE : public WATER_FREE_SURFACE_2D_EXAMPLE<T,RW>
{
public:
    T radius;
    VECTOR_2D<T> velocity_object, center;

    WATER_OVER_CIRCLE() {

        // grid
        m=100;n=100;
        xmin=0;xmax=1;ymin=0;ymax=1;
        // time
        initial_time=0;
        first_frame=0;last_frame=720;
        frame_rate=3*24;
        //circle specs.
        center=VECTOR_2D<T>(xmin+(T).5*(xmax-xmin), ymin+(T).5*(ymax-ymin));radius=(T).1*(xmax-xmin);
        velocity_object = VECTOR_2D<T>((T)0,(T)0); slip_coefficient = (T)1.;
        //output
        write_levelset=true;write_velocity=true;write_particles=true;
        matlab_directory="Water_Over_Circle/matlab";data_directory="Water_Over_Circle/output";
        write_matlab_file=false;
    }

    virtual ~WATER_OVER_CIRCLE() {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(){
    int m=grid.m,n=grid.n;
    T dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        T x=grid.x(i),y=grid.y(j);
        if((fabs(x-(xmin+xmax)/2) <= .1*(xmax-xmin)) && y >= grid.y(n)-3.*dy) phi(i,j)=-dx;
        else phi(i,j)=dx;}}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
void Get_Source_Velocities(const T time){
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        T x=grid.x(i),y=grid.y(j),dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
        if((fabs(x-(xmin+xmax)/2) <= .1*(xmax-xmin)) && y >= grid.y(n)-3.*dy){
            psi_N(i,j)=true;V_fixed(i,j)=VECTOR_2D<T>(0,-4);}}}
//#####################################################################
// Function Adjust_Phi_With_Sources
//#####################################################################
void Adjust_Phi_With_Sources(const T time){
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        T x=grid.x(i),y=grid.y(j),dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
        if((fabs(x-(xmin+xmax)/2) <= .1*(xmax-xmin)) && y >= grid.y(n)-3.*dy) phi(i,j)=((time<1.0/frame_rate*4)?-dx:dx);}}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void Get_Source_Reseed_Mask(ARRAY<bool,VECTOR<int,2> >*& reseed, const T time){
    if(reseed) delete reseed;reseed=new ARRAY<bool,VECTOR<int,2> >(1,grid.m,1,grid.n);

    T dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    T padding=3*dx;
    ARRAY<bool,VECTOR<int,2> >::copy(false,*reseed);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        T x=grid.x(i),y=grid.y(j);
        if((fabs(x-(xmin+xmax)/2) <= .1*(xmax-xmin)+padding) && y >= grid.y(n)-3.*dy-2*padding){
              if(time < 1.0/frame_rate*4) (*reseed)(i,j)=true;}}}
//#####################################################################
// Function Adjust_Phi_And_Get_Velocities_For_Objects
//#####################################################################
virtual void Adjust_Phi_And_Get_Velocities_For_Objects(ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >& V_object,const T time) 
{
    // creating a levelset for the circle
    ARRAY<T,VECTOR<int,2> > phi_object(-2,m+3,-2,n+3);
    VECTOR_2D<T> updated_center=center+velocity_object*time;
    for(int i=-2;i<=m+3;i++) for(int j=-2;j<=n+3;j++) phi_object(i,j)=(grid.X(i,j)-updated_center).Magnitude()-radius;
    
    Extrapolate_Velocity_Into_Objects(grid,phi_object,u,v,velocity_object,V_object,psi_N,time);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
void Delete_Particles_Inside_Objects(HEAVY_PARTICLES<T,VECTOR_2D<T> >& particles,const T time)
{
    VECTOR_2D<T> updated_center = center+velocity_object*time;
    for(int k=0;k<particles.array_size;k++) if(particles.active(k)) 
        if(sqrt(sqr(particles.X(k).x-updated_center.x)+sqr(particles.X(k).y-updated_center.y)) < radius) particles.Delete_Particle(k);
}
//#####################################################################
};      
}
#endif
