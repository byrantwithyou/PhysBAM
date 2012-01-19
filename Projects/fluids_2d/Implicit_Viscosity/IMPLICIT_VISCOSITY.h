//#####################################################################
// Copyright 2002, 2003, Doug Enright, Ronald Fedkiw, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_VISCOSITY
//##################################################################### 
//
//#####################################################################
// Enright - June 19, 2002
// Nguyen - June 27, 2002
// Fedkiw - July 19, 2003
//#####################################################################
#ifndef __IMPLICIT_VISCOSITY__
#define __IMPLICIT_VISCOSITY__

#include "../WATER_FREE_SURFACE_2D_EXAMPLE.h"
namespace PhysBAM{

template<class T,class RW=T>
class IMPLICIT_VISCOSITY:public WATER_FREE_SURFACE_2D_EXAMPLE<T,RW>
{
public:
    IMPLICIT_VISCOSITY()
    {
        frame_rate=200;
        last_frame=1000;//(int)(10*frame_rate);
        m=100;
        n=100;
        domain_walls[1][1]=true;domain_walls[1][2]=true;domain_walls[2][1]=true;domain_walls[2][2]=false;
        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
        //viscosity=1e6*.001137;
        //viscosity=100;
        implicit_viscosity=false;
        levelset_substeps=(T)1.9;
        cfl=2;
        incompressible_iterations=60;
        implicit_viscosity_iterations=60;
        write_levelset=write_velocity=write_particles=write_removed_positive_particles=write_removed_negative_particles=true;
        matlab_directory="Implicit_Viscosity/matlab1";data_directory="Implicit_Viscosity/output";
    }

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    for(int i=0;i<grid.m;i++) for(int j=0;j<grid.n;j++) phi(i,j)=minmag(grid.y(j)-(T).5+grid.dy/2,(T)sqrt(sqr(grid.x(i)-(T).5)+sqr(grid.y(j)-(T).75))-(T).125);
//    for(int i=0;i<m;i++) for(int j=0;j<n;j++) phi(i,j)=grid.y(j)-.5;
}
/*
//#####################################################################
// Function Specify_Variable_Viscosity_Values
//#####################################################################
virtual void Specify_Variable_Viscosity_Values(const GRID<TV>& grid,ARRAY<T,VECTOR<int,2> >& viscosity)
{
    for(int i=0;i<=grid.m+1;i++) for(int j=0;j<=grid.n+1;j++) 
        if(grid.x(i) > .5) viscosity(i,j)=1000;//(T)1e8*.001137;
        else viscosity(i,j)=0;//(T)(1e6*.001137*1e-6);
}*/
//#####################################################################
// Function Write_Matlab_Data_File
//#####################################################################
void Write_Matlab_Data_File(const int stepnumber)
{
    WATER_FREE_SURFACE_2D_EXAMPLE<T>::Write_Matlab_Data_File(stepnumber,grid,phi,p,u,v,particle_levelset);
                                                       
    int i,j;
    MATLAB_OUTPUT matlab_output;
    ARRAY<T,VECTOR<int,2> > output(1,m,1,n);
    char file_name[256];

    ARRAY<T,VECTOR<int,2> > phi_temp(phi);LEVELSET_2D<T> levelset_temp(grid,phi_temp);levelset_temp.Reinitialize();//Fast_Marching_Method();
    levelset_temp.Compute_Curvature();
    for(i=0;i<m;i++) for(j=0;j<n;j++) if(fabs(phi_temp(i,j)) > 5*grid.dx) (*levelset_temp.curvature)(i,j)=0; 
    
    for(i=0;i<m;i++) for(j=0;j<n;j++) output(i,j)=(*levelset_temp.curvature)(i,j);
    sprintf(file_name,"%s/curvature",matlab_directory);
    matlab_output.Write_Output_File(file_name,output,stepnumber);

    sprintf(file_name,"%s/phi_temp",matlab_directory);
    matlab_output.Write_Output_File(file_name,phi_temp,stepnumber);
}
//#####################################################################
};      
}
#endif
