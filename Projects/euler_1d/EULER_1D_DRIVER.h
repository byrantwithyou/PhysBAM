//#####################################################################
// Copyright 2003-2007, Doug Enright, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_1D_DRIVER
//#####################################################################
//
//#####################################################################
// Enright - August 29, 2003
//#####################################################################
#ifndef __EULER_1D_DRIVER__
#define __EULER_1D_DRIVER__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/EULER_UNIFORM.h>

#include "EULER_1D_EXAMPLE.h"

namespace PhysBAM{

template<class T>
class EULER_1D_DRIVER
{
public:
    EULER_1D_EXAMPLE<T>& example;
    T time,time_per_frame;
    int total_steps;

    ARRAY<T,3,VECTOR<int,1> > U;
    EULER_UNIFORM<TV> euler;

public:
    EULER_1D_DRIVER(EULER_1D_EXAMPLE<T>& example_input):example(example_input),euler(U)
    {
        Initialize();
    }

//#####################################################################
// Function Initialize
//#####################################################################
void Initialize()
{
    // initialize time
    total_steps=(int)((example.final_time-example.initial_time)*example.frame_rate);
    time_per_frame=(example.final_time-example.initial_time)/total_steps;
    time=example.initial_time+example.restart_step_number*time_per_frame;
    
    std::cout << "total step = " << total_steps << " time_per_frame = " << time_per_frame << " time = " << time << std::endl;

    // initialize grids
    euler.grid.Initialize(example.m,example.xmin,example.xmax);
    int m=example.m;
    U.Resize(1,m);

    euler.eigensystems[1]=new EULER_1D_EIGENSYSTEM_F<T>();

    // set up specialized stuff . . .
    if(example.eos) euler.Set_Custom_Equation_Of_State(*example.eos);
    
    if(example.conservation_method) euler.Set_Custom_Conservation(*example.conservation_method);
    
    if(example.u_boundary) euler.Set_Custom_Boundary(*example.u_boundary);
    
    if(example.set_max_time_step) euler.Set_Max_Time_Step(example.max_time_step);
    euler.Set_CFL_Number(example.cfl_number);
        
    //set spatial order 
    euler.conservation->Set_Order(example.spatial_order);
    
    //initialize u
    example.Initialize_U(euler.grid,U);

    example.Write_Matlab_Data_File(0,euler.grid,U); //going to need the equation of state info for auxillary outputs
}
//#####################################################################
// Function Execute_Main_Program
//#####################################################################
void Execute_Main_Program()
{
    for(int k=example.restart_step_number+1;k<=total_steps;k++){
        std::cout << std::endl;
        Advance_One_Frame(example.initial_time+k*time_per_frame);
        std::cout << "Finished frame " << k << " - time " << example.initial_time+k*time_per_frame << std::endl;
        example.Write_Matlab_Data_File(k,euler.grid,U);}
}
//#####################################################################
// Function Advance_One_Frame
//#####################################################################
void Advance_One_Frame(const T stopping_time)
{
    int substep=0;bool done=false;
    while(!done){substep++;
        T dt=euler.cfl_number*euler.CFL();
        if(time+dt >= stopping_time){dt=stopping_time-time;done=true;}
        std::cout << "substep = " << substep << ", dt = " << dt << " stopping_time = " << stopping_time << std::endl;
        Advance_One_Time_Step(dt,true);}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
void Advance_One_Time_Step(const T dt,const bool verbose) 
{
    RUNGEKUTTA<T,ARRAY<T,3,VECTOR<int,1> >,GRID<TV> > rungekutta_u(U);
    rungekutta_u.Set_Grid_And_Boundary_Condition(euler.grid,*(euler.boundary));
    rungekutta_u.Set_Order(example.rungekutta_order);
    rungekutta_u.Set_Time(time);
    
    rungekutta_u.Start(dt);

    for (int kk=1; kk <= rungekutta_u.order; kk++) {
        euler.Euler_Step(dt);
        time=rungekutta_u.Main();
        //fix up bcs?
        euler.boundary->Apply_Boundary_Condition(euler.grid,U);}
}
//#####################################################################
};
}
#endif
