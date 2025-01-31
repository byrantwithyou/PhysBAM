//#####################################################################
// Copyright 2003 Doug Enright.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LANEY_ST1
//#####################################################################
//
//#####################################################################
// Enright - September 9, 2003
//#####################################################################
#ifndef __LANEY_ST1__
#define __LANEY_ST1__

#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <Compressible/Equations_Of_State/EOS_GAMMA.h>
#include "../EULER_1D_EXAMPLE.h"
#include <fstream>
#include <iostream>
using namespace PhysBAM;

namespace PhysBAM{

template <class T>
class LANEY_ST1:public EULER_1D_EXAMPLE<T>
{
public:    
    LANEY_ST1()                                                                      
    {
        //grid
        m=50; xmin=-10.; xmax=10.;
        //time
        final_time=24*4.276e-4;if(!this->user_frame_rate) frame_rate=24./final_time;
        cfl_number=1.;
        //custom stuff . . . 
        eos = new EOS_GAMMA<T>;
        //conservation_method = new CONSERVATION_ENO_RF<T>;
        strcpy(output_directory,"Laney_ST1/matlab");
    }
    
    ~LANEY_ST1() {}

//#####################################################################
    void Initialize_U(const GRID<TV>& grid, ARRAY<T,VECTOR<int,1> >& u) override;
    void Write_Matlab_Data_File(const int stepnumber, const GRID<TV>& grid, const ARRAY<T,VECTOR<int,1> >& u) override;
//#####################################################################    
};
//#####################################################################
// Function Intialize_U
//#####################################################################
template<class T> void LANEY_ST1<T>::
Initialize_U(const GRID<TV>& grid, ARRAY<T,VECTOR<int,1> >& u)
{
    EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(eos);

    //initialize grid variables
    //1 =- density, 2 == momentum, 3 == total energy
    for(int i=0;i<m;i++) 
        if(grid.x(i) < 0) {u(1,i)=1.0; u(2,i) = 0.0; u(3,i) = 100000/(tmp_eos->gamma,tmp_eos->gamma-1);} 
        else {u(1,i) = 0.125; u(2,i) = 0.0; u(3,i) = 10000/(tmp_eos->gamma-1);}

}
//#####################################################################
// Function Write_Matlab_Data_File
//#####################################################################
template<class T> void LANEY_ST1<T>::
Write_Matlab_Data_File(const int stepnumber, const GRID<TV>& grid, const ARRAY<T,VECTOR<int,1> >& u)
{
        int i,m=grid.m;
        ARRAY<T,VECTOR<int,1> > output(1,m),output1(1,m);
        char file_name[256];

        MATLAB_OUTPUT matlab_output;

        EOS_GAMMA<T> *tmp_eos = dynamic_cast<EOS_GAMMA<T>*>(eos);

        //output header
        sprintf(file_name,"%s/header",viewer_dir.output_directory);
        matlab_output.Write_Header_File(file_name,grid,stepnumber);

        // output primitive variables . . . . 
        //density
        sprintf(file_name,"%s/density",viewer_dir.output_directory);
        for(i=0;i<m;i++) output(i) = u(1,i);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //velocity
        for(i=0;i<m;i++) output(i) = u(2,i)/u(1,i);
        sprintf(file_name,"%s/velocity",viewer_dir.output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //pressure
        for(i=0;i<m;i++) output1(i) = (tmp_eos->gamma-1.0)*(u(3,i)-.5*u(2,i)*output(i));
        sprintf(file_name,"%s/pressure",viewer_dir.output_directory);
        matlab_output.Write_Output_File(file_name,output1,stepnumber);
        //entropy
        for(i=0;i<m;i++) output(i) = tmp_eos->Cv*log(output1(i)/pow(u(1,i),tmp_eos->gamma));
        sprintf(file_name,"%s/entropy",viewer_dir.output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //speed of sound
        for(i=0;i<m;i++) output(i) = sqrt(tmp_eos->gamma*output1(i)/u(1,i));
        sprintf(file_name,"%s/sound_speed",viewer_dir.output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
        //Mach number
        for(i=0;i<m;i++) output(i) = (u(2,i)/u(1,i))/output(i);
        sprintf(file_name,"%s/mach_number",viewer_dir.output_directory);
        matlab_output.Write_Output_File(file_name,output,stepnumber);
}
//#####################################################################    
}
#endif
