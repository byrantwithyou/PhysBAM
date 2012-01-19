//#####################################################################
// Copyright 2003 Doug Enright.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHU_OSHER_ST
//#####################################################################
//
//#####################################################################
// Enright - September 9, 2003
//#####################################################################
#ifndef __SHU_OSHER_ST__
#define __SHU_OSHER_ST__

#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <fstream>
#include <iostream>
#include "../EULER_1D_EXAMPLE.h"
using namespace PhysBAM;

namespace PhysBAM{

template <class T>
class SHU_OSHER_ST:public EULER_1D_EXAMPLE<T>
{
public:    
    SHU_OSHER_ST()                                                                      
    {
        //grid
        m=101; xmin=-5.; xmax=5.;
        //time
        initial_time=0.;final_time=60*3.38e-2;frame_rate=60./final_time;
        cfl_number=1.;
        //custom stuff . . . 
        eos = new EOS_GAMMA<T>;
        conservation_method = new CONSERVATION_ENO_RF<T>;
        strcpy(output_directory,"Shu_Osher_ST/matlab");
    }
    
    ~SHU_OSHER_ST() {}

//#####################################################################
    void Initialize_U(const GRID<TV>& grid, ARRAY<T,VECTOR<int,1> >& u) PHYSBAM_OVERRIDE;
    void Write_Matlab_Data_File(const int stepnumber, const GRID<TV>& grid, const ARRAY<T,VECTOR<int,1> >& u) PHYSBAM_OVERRIDE;
//#####################################################################    
};
//#####################################################################
// Function Intialize_U
//#####################################################################
template<class T> void SHU_OSHER_ST<T>::
Initialize_U(const GRID<TV>& grid, ARRAY<T,VECTOR<int,1> >& u)
{
    EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(eos);

    //initialize grid variables
    //1 =- density, 2 == momentum, 3 == total energy
    for(int i=0;i<m;i++)
        if(grid.x(i) < 0) {u(1,i)=1.; u(2,i) = 0.0; u(3,i) = 1./(tmp_eos->gamma-1);}
        else {u(1,i)=0.125; u(2,i)=0.0; u(3,i)= .1/(tmp_eos->gamma-1.);}
}
//#####################################################################
// Function Write_Matlab_Data_File
//#####################################################################
template<class T> void SHU_OSHER_ST<T>::
Write_Matlab_Data_File(const int stepnumber, const GRID<TV>& grid, const ARRAY<T,VECTOR<int,1> >& u)
{
        int i,m=grid.m;
        ARRAY<T,VECTOR<int,1> > output(1,m),output1(1,m);
        char file_name[256];

        MATLAB_OUTPUT matlab_output;

        EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(eos);

        //output header
        sprintf(file_name,"%s/header",output_directory);
        matlab_output.Write_Header_File(file_name,grid,stepnumber);

        // output primitive variables . . . . 
        //density
        sprintf(file_name,"%s/density",output_directory);
        for(i=1;i<=m;i++) output(i) = u(1,i);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //velocity
        for(i=1;i<=m;i++) output(i) = u(2,i)/u(1,i);
        sprintf(file_name,"%s/velocity",output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //pressure
        for(i=1;i<=m;i++) output1(i) = (tmp_eos->gamma-1.0)*(u(3,i)-.5*u(2,i)*output(i));
        sprintf(file_name,"%s/pressure",output_directory);
        matlab_output.Write_Output_File(file_name,output1,stepnumber);
        //entropy
        for(i=1;i<=m;i++) output(i) = tmp_eos->Cv*log(output1(i)/pow(u(1,i),tmp_eos->gamma));
        sprintf(file_name,"%s/entropy",output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //speed of sound
        for(i=1;i<=m;i++) output(i) = sqrt(tmp_eos->gamma*output1(i)/u(1,i));
        sprintf(file_name,"%s/sound_speed",output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //Mach number
        for(i=1;i<=m;i++) output(i) = (u(2,i)/u(1,i))/output(i);
        sprintf(file_name,"%s/mach_number",output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
}
//#####################################################################    
}
#endif
