//#####################################################################
// Copyright 2002, Doug Enright, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WATER_INTO_BOX
//##################################################################### 
//
//#####################################################################
// Enright - June 19, 2002
// Nguyen - June 21, 2002
//#####################################################################
#ifndef __WATER_INTO_BOX__
#define __WATER_INTO_BOX__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Particles/HEAVY_PARTICLES.h>

#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>

namespace PhysBAM{

class WATER_INTO_BOX : public EXAMPLE
{
public:
    WATER_INTO_BOX() {

        // grid
        m=100;n=100;
        xmin=0;xmax=1;ymin=0;ymax=1;
        // time
        start_time=0;final_time=10;
        frame_rate=3*24;
        // particle
        use_deep_particles=0;
        //output
        matlab_directory_name="Water_Into_Box/matlab";}

    ~WATER_INTO_BOX() {}

//#####################################################################
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void WATER_INTO_BOX::Initialize_Phi(){
    int i,j;
    int m=grid.m,n=grid.n;
    double dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    for(i=1;i<=m;i++) for(j=1;j<=n;j++){
        double x=grid.x(i),y=grid.y(j);
        if((fabs(x-(xmin+xmax)/2) <= .1*(xmax-xmin)) && y >= grid.y(n)-3.*dy) phi(i,j)=-dx;
        else phi(i,j)=dx;}}
//#####################################################################
// Function Get_Sources
//#####################################################################
void WATER_INTO_BOX::Get_Sources(double time){
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        double x=grid.x(i),y=grid.y(j),dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
        if((fabs(x-(xmin+xmax)/2) <= .1*(xmax-xmin)) && y >= grid.y(n)-3.*dy){
            phi(i,j)=-dx;psi_N(i,j)=1;u_fixed(i,j)=0;v_fixed(i,j)=-2;}}}
//#####################################################################
// Function Get_Objects
//#####################################################################
void WATER_INTO_BOX::Get_Objects(double time){
    // stationary circle
    double epsilon=.0001*max(grid.dx,grid.dy);
    double xmin=grid.xmin,xmax=grid.xmax;
    double x_center=xmin+.4*(xmax-xmin),y_center=ymin+.1*(ymax-ymin),radius=.1*(xmax-xmin);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++)
        if(sqrt(sqr(grid.x(i)-x_center)+sqr(grid.y(j)-y_center)) < radius){
            // extrapolate phi inward
            double nx=grid.x(i)-x_center,ny=grid.y(j)-y_center,magnitude=sqrt(sqr(nx)+sqr(ny));nx/=magnitude;ny/=magnitude;
            double distance=radius-magnitude+epsilon;
            double closest_point_x=grid.x(i)+distance*nx,closest_point_y=grid.y(j)+distance*ny;
            phi(i,j)=interpolate(grid,phi,closest_point_x,closest_point_y); 
            // set the velocity
            psi_N(i,j)=1;u_fixed(i,j)=0;v_fixed(i,j)=0;}}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
void WATER_INTO_BOX::Get_Source_Reseed_Mask(ARRAY<int,VECTOR<int,2> >& reseed_mask,double time){
    double dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    double padding=3*dx;
    copy(0,reseed_mask);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        double x=grid.x(i),y=grid.y(j);
        if((fabs(x-(xmin+xmax)/2) <= .1*(xmax-xmin)+padding) && y >= grid.y(n)-3.*dy-2*padding){
              (reseed_mask)(i,j)=1;}}}
//#####################################################################
// Function Delete_Particles_Inside_Geometry
//#####################################################################
void WATER_INTO_BOX::Delete_Particles_Inside_Geometry(){
    // stationary circle
    int k;
    double dx=grid.dx,dy=grid.dy,xmin=grid.xmin,xmax=grid.xmax;
    double epsilon=.0001*max(dx,dy);
    double x_center=xmin+.4*(xmax-xmin),y_center=ymin+.1*(ymax-ymin),radius=.1*(xmax-xmin);
    // delete negative particles inside object
    for(k=1;k<=particle_levelset.negative_particles.array_size;k++) if(particle_levelset.negative_particles.active(k)){
        double x=particle_levelset.negative_particles.x(k),y=particle_levelset.negative_particles.y(k);
        double distance=sqrt(sqr(x-x_center)+sqr(y-y_center)); 
        if(distance < radius-epsilon) particle_levelset.negative_particles.Delete_Particle(k);} // include buffer zone
        // delete positive particles inside object
    for(k=1;k<=particle_levelset.positive_particles.array_size;k++) if(particle_levelset.positive_particles.active(k)){
        double x=particle_levelset.positive_particles.x(k),y=particle_levelset.positive_particles.y(k);
        double distance=sqrt(sqr(x-x_center)+sqr(y-y_center)); 
        if(distance < radius-epsilon)  particle_levelset.positive_particles.Delete_Particle(k);}} // include buffer zone
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
