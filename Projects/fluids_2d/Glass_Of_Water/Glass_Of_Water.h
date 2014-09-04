//#####################################################################
// Copyright 2002, Doug Enright, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GLASS_OF_WATER
//##################################################################### 
//
//#####################################################################
// Enright - June 19, 2002
// Nguyen - June 27, 2002
//#####################################################################
#ifndef __GLASS_OF_WATER__
#define __GLASS_OF_WATER__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Particles/HEAVY_PARTICLES.h>

#include <Tools/Ordinary_Differential_Equations/EXAMPLE.h>

namespace PhysBAM{

class GLASS_OF_WATER:public EXAMPLE
{
public:
    GLASS_OF_WATER() {

        // grid
        m=100;n=100;
        xmin=0;xmax=1;ymin=0;ymax=1;
        // time
        start_time=0;final_time=10;
        frame_rate=3*24;
        // particle
        use_deep_particles=0;
        //output
        matlab_directory_name="Glass_of_Water/matlab";}

    ~GLASS_OF_WATER() {}

//#####################################################################
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void GLASS_OF_WATER::Initialize_Phi(){
    int i,j;
    int m=grid.m,n=grid.n;
    double dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    for(i=0;i<m;i++) for(j=0;j<n;j++){
        double x=grid.x(i),y=grid.y(j);
        if((fabs(x-(xmin+xmax)/2+.05*(xmax-xmin)) <= .05*(xmax-xmin)) && y >= grid.y(n)-3.*dy) (phi)(i,j)=-dx;
        else (phi)(i,j)=dx;}}
//#####################################################################
// Function Get_Sources
//#####################################################################
void GLASS_OF_WATER::Get_Sources(double time){
    double dx=grid.dx,dy=grid.dy;
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        double x=grid.x(i),y=grid.y(j),xmin=grid.xmin,xmax=grid.xmax;
        if((fabs(x-(xmin+xmax)/2+.05*(xmax-xmin)) <= .05*(xmax-xmin)) && y >= grid.y(n)-3.*dy){
            if(time <= 4){(phi)(i,j)=-dx;(psi_N)(i,j)=1;(u_fixed)(i,j)=-1;(v_fixed)(i,j)=-1;}}}}
//#####################################################################
// Function Get_Objects
//#####################################################################

void GLASS_OF_WATER::Get_Objects(double time){
    // stationary circle
    double epsilon=.0001*max(grid.dx,grid.dy);
    double xmin=grid.xmin,xmax=grid.xmax;
    // walls
    double x_center=(xmin+xmax)/2,wall_radius=(xmax-xmin)/4;
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
        double radius=fabs(grid.x(i)-x_center);
        if(radius > wall_radius){
        // extrapolate phi inward from outside the object
            double nx=(x_center-grid.x(i))/fabs(grid.x(i)-x_center);
            double magnitude=radius-wall_radius+epsilon;
            double closest_point_x=grid.x(i)+magnitude*nx,closest_point_y=grid.y(j);
            (phi)(i,j)=interpolate(grid,phi,closest_point_x,closest_point_y); 
            // set the velocity
            (psi_N)(i,j)=1;(u_fixed)(i,j)=0;(v_fixed)(i,j)=0;}}}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void GLASS_OF_WATER::Get_Source_Reseed_Mask(ARRAY<int,VECTOR<int,2> >& reseed_mask,double time){
    double dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    double padding=3*dx;
    copy(0,reseed_mask);
    for(int i=0;i<m;i++) for(int j=0;j<n;j++){
      double x=grid.x(i),y=grid.y(j);
      if((fabs(x-(xmin+xmax)/2+.05*(xmax-xmin)) <= .05*(xmax-xmin)+padding) && y >= grid.y(n)-3.*dy-2*padding){
          (reseed_mask)(i,j)=1;}}}
//#####################################################################
// Function Delete_Particles_Inside_Geometry
//#####################################################################
void GLASS_OF_WATER::Delete_Particles_Inside_Geometry(){
    // stationary circle
    int k;
    double dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    double epsilon=.0001*max(dx,dy);
    // walls
    double x_center=(xmin+xmax)/2,wall_radius=(xmax-xmin)/4;
    // delete negative particles inside object
    for(k=0;k<particle_levelset.negative_particles.array_size;k++) if(particle_levelset.negative_particles.active(k)){
        double x=particle_levelset.negative_particles.x(k),y=particle_levelset.negative_particles.y(k);
        double radius=fabs(x-x_center);
        if(radius > wall_radius+epsilon) particle_levelset.negative_particles.Delete_Particle(k);} // include buffer zone
    // delete positive particles inside object
    for(k=0;k<particle_levelset.positive_particles.array_size;k++) if(particle_levelset.positive_particles.active(k)){
        double x=particle_levelset.positive_particles.x(k),y=particle_levelset.positive_particles.y(k);
        double radius=fabs(x-x_center);
        if(radius > wall_radius+epsilon) particle_levelset.positive_particles.Delete_Particle(k);}}// include buffer zone
//#####################################################################
// Function Write_Matlab_Data_File
//#####################################################################
void Write_Matlab_Data_File(const char *directory, const int stepnumber){
    MATLAB_OUTPUT matlab_output;
    ARRAY_2D output(1,m,1,n);
    char file_name[256];
    
    // main data
    sprintf(file_name,"%s/ooo_header",directory);
    matlab_output.Write_Header_File(file_name,grid,stepnumber);
    sprintf(file_name,"%s/ooo_levelset",directory);
    matlab_output.Write_Output_File(file_name,phi,stepnumber);
}
//#####################################################################
};      
}
#endif
