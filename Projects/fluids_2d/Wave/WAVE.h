//#####################################################################
// Copyright 2002, Doug Enright, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class WAVE
//##################################################################### 
//
//#####################################################################
// Enright - June 19, 2002
// Nguyen - June 24, 2002
//#####################################################################
#ifndef __WAVE__
#define __WAVE__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Particles/HEAVY_PARTICLES.h>

#include <Tools/Ordinary_Differential_Equations/EXAMPLE.h>

namespace PhysBAM{

class WAVE:public EXAMPLE
{
public:
    WAVE() {

        // grid
        m=400;n=70;
        xmin=-50;xmax=150;ymin=0;ymax=35;
        // time
        start_time=0;final_time=100;
        frame_rate=3*24;
        // particle
        use_deep_particles=0;
        // backdoor
        backdoor_id=2;
        //output
        matlab_directory_name="Wave/matlab";}

    ~WAVE() {}

//#####################################################################
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void WAVE::Initialize_Phi(){
    int i,j;
    int m=grid.m,n=grid.n;
    for(i=0;i<m;i++) for(j=0;j<n;j++){
        double x=grid.x(i),y=grid.y(j);
        double g=9.8,d=10,h=10,argument=sqrt(3*h/(4*cube(d)))*x;
        double sech=2/(exp(argument)+exp(-argument)),height=d+h*sqr(sech);
        (phi)(i,j)=y-height;}}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void WAVE::Initialize_Velocities(){
    int i,j;
    for(i=0;i<u_grid.m;i++) for(j=0;j<u_grid.n;j++){
        double x=u_grid.x(i),y=u_grid.y(j);
        double g=9.8,d=10,h=10,argument=sqrt(3*h/(4*cube(d)))*x;
        double sech=2/(exp(argument)+exp(-argument));
        (u)(i,j)=sqrt(g*d)*h/d*sqr(sech);}
    for(i=0;i<v_grid.m;i++) for(j=0;j<v_grid.n;j++){
        double x=v_grid.x(i),y=v_grid.y(j);
        double g=9.8,d=10,h=10,argument=sqrt(3*h/(4*cube(d)))*x;
        double sech=2/(exp(argument)+exp(-argument)),tanh=(exp(argument)-exp(-argument))/(exp(argument)+exp(-argument));
        (v)(i,j)=sqrt(3*g*d)*sqrt(cube(h/d))*y/d*sqr(sech)*tanh;}}
//#####################################################################
// Function Get_Objects
//#####################################################################
void WAVE::Get_Objects(double time){
    double epsilon=.0001*max(grid.dx,grid.dy);
    // sloped bottom - screws up the step, but ok if the water surface is never near there!
    for(int i=0;i<m;i++) for(int j=0;j<n;j++)
        if(grid.x(i) > 50 && 14*grid.y(j)-grid.x(i)+0 < 0){
            // extrapolate phi inward
            double length=sqrt(197),nx=-1/length,ny=14/length;
            double magnitude=-(14*grid.y(j)-grid.x(i)+0)/sqrt(197)+epsilon;
            double closest_point_x=grid.x(i)+magnitude*nx,closest_point_y=grid.y(j)+magnitude*ny;
            (phi)(i,j)=interpolate(grid,phi,closest_point_x,closest_point_y); 
            // set velocity
            (psi_N)(i,j)=1;(u_fixed)(i,j)=0;(v_fixed)(i,j)=0;}}
//#####################################################################
// Function Delete_Particles_Inside_Geometry
//#####################################################################
void WAVE::Delete_Particles_Inside_Geometry(PARTICLE_LEVELSET_2D& particle_levelset){
    double epsilon=.0001*max(grid.dx,grid.dy);
    // sloped bottom
    int k;
    // delete negative particles inside object
    for(k=0;k<particle_levelset.negative_particles.array_size;k++) if(particle_levelset.negative_particles.active(k)){
        double x=particle_levelset.negative_particles.x(k),y=particle_levelset.negative_particles.y(k);
        if(x > 50 && (14*y-x+0)/sqrt(197) < -epsilon) particle_levelset.negative_particles.Delete_Particle(k);} // include buffer zone
    // delete positive particles inside object
    for(k=0;k<particle_levelset.positive_particles.array_size;k++) if(particle_levelset.positive_particles.active(k)){
        double x=particle_levelset.positive_particles.x(k),y=particle_levelset.positive_particles.y(k);
        if(x > 50 && (14*y-x+0)/sqrt(197) < -epsilon) particle_levelset.positive_particles.Delete_Particle(k);}} // include buffer zone
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
