//#####################################################################
// Copyright 2003-2007 Doug Enright, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_1D_EXAMPLE
//#####################################################################
//
//#####################################################################
// Enright - August 29, 2003
//#####################################################################
#ifndef __EULER_1D_EXAMPLE__
#define __EULER_1D_EXAMPLE__

#include <fstream>
#include <iostream>

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Debug_Output/MATLAB_OUTPUT.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <Compressible/Equations_Of_State/EOS.h>

namespace PhysBAM{

template <class T>
class EULER_1D_EXAMPLE
{
public:
    //time
    T final_time;
    T frame_rate;                        
    
    int restart_step_number;
    
    T cfl_number;
    bool set_max_time_step;
    double max_time_step;
    
    int spatial_order;
    int rungekutta_order;
    
    //grid
    int m;
    T xmin, xmax;
    
    //allow for custom boundary conditions & eos's
    BOUNDARY<TV,T,TV::m+2>* u_boundary;
    EOS<T>* eos;
    CONSERVATION<T,3>* conservation_method;
    
    char output_directory[256];
    
    EULER_1D_EXAMPLE()
        :final_time((T)1),frame_rate((T)24),restart_step_number(0),cfl_number((T)1),
        set_max_time_step(false),max_time_step(1.e8),spatial_order(3),rungekutta_order(3),m(10),xmin((T)0),xmax((T)1),u_boundary(0),eos(0),conservation_method(0)
    {
        strcpy(output_directory,"matlab");
    }
    
    virtual ~EULER_1D_EXAMPLE() {}

//#####################################################################
    virtual void Initialize_U(const GRID<TV>& grid, ARRAY<T,3,VECTOR<int,1> >& u);
    virtual void Write_Matlab_Data_File(const int stepnumber, const GRID<TV>& grid, const ARRAY<T,3,VECTOR<int,1> >& u);
//#####################################################################    
};
//#####################################################################
// Function Intialize_U
//#####################################################################
template<class T> void EULER_1D_EXAMPLE<T>::
Initialize_U(const GRID<TV>& grid, ARRAY<T,3,VECTOR<int,1> >& u)
{
    // could be problematic to reinitialize if phi is always postive or negative, so maybe have this function call sources
    std::cout << "Careful - Need Initial Data For Euler" << std::endl;
    exit(1);
}
//#####################################################################
// Function Write_Matlab_Data_File
//#####################################################################
template<class T> void EULER_1D_EXAMPLE<T>::
Write_Matlab_Data_File(const int stepnumber, const GRID<TV>& grid, const ARRAY<T,3,VECTOR<int,1> >& u)
{
    int i,m=grid.m;
    ARRAY<T,VECTOR<int,1> > output(1,m);
    char file_name[256];

    MATLAB_OUTPUT matlab_output;

    //output header
    sprintf(file_name,"%s/header",viewer_dir.output_directory);
    matlab_output.Write_Header_File(file_name,grid,stepnumber);

    // output primitive variables . . . . 
    //density
    sprintf(file_name,"%s/u1",viewer_dir.output_directory);
    for(i=0;i<m;i++) output(i) = u(1,i);
    matlab_output.Write_Output_File(file_name,output,stepnumber);
    //momentum
    for(i=0;i<m;i++) output(i) = u(2,i);
    sprintf(file_name,"%s/u2",viewer_dir.output_directory);
    matlab_output.Write_Output_File(file_name,output,stepnumber);
    //internal energy
    for(i=0;i<m;i++) output(i) = u(3,i);
    sprintf(file_name,"%s/u3",viewer_dir.output_directory);
    matlab_output.Write_Output_File(file_name,output,stepnumber);
}
//#####################################################################    
}
#endif
