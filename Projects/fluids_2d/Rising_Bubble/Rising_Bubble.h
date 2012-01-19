//#####################################################################
// Copyright 2002, Doug Enright, Duc Nguyen.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RISING_BUBBLE
//##################################################################### 
//
//#####################################################################
// Enright - June 19, 2002
// Nguyen - June 27, 2002
//#####################################################################
#ifndef __RISING_BUBBLE__
#define __RISING_BUBBLE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Particles/HEAVY_PARTICLES.h>

#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>

namespace PhysBAM{

class RISING_BUBBLE : public EXAMPLE
{
public:
    RISING_BUBBLE() {

        // grid
        m=100;n=100;
        xmin=0;xmax=1;ymin=0;ymax=1;
        // time
        start_time=0;final_time=10;
        frame_rate=3*24;
        // particle
        use_deep_particles=0;
        //output
        matlab_directory_name="rising_Bubble/matlab";}

    ~RISING_BUBBLE() {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void RISING_BUBBLE::Initialize_Phi(){
    int i,j;
    int m=grid.m,n=grid.n;
    for(i=0;i<m;i++) for(j=0;j<n;j++) phi(i,j)=minmag(grid.y(j)-.5,.125-sqrt(sqr(grid.x(i)-.5)+sqr(grid.y(j)-.25)));}
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
