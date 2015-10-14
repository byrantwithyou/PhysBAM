//#####################################################################
// Copyright 2002, 2003, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_1D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_1D class.
// Input U as 3 by (1,m) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m).
// When psi=true, solve the equaitions. 
// When psi=false, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_1D__
#define __EULER_1D__

#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Euler_Equations/EULER_1D_EIGENSYSTEM_F.h>
namespace PhysBAM{

template<class T_input>
class EULER_1D:public EULER<VECTOR<T_input,1> >
{
    typedef T_input T;typedef VECTOR<T,1> TV;typedef VECTOR<T,3> TV_DIMENSION;
protected:
    using EULER<TV>::cut_out_grid;using EULER<TV>::max_time_step;
public:
    using EULER<TV>::boundary;using EULER<TV>::conservation;using EULER<TV>::eos;
    
    GRID<TV> grid;
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& U;         // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,1> >* psi_pointer; // defines cut out grid
    EULER_1D_EIGENSYSTEM_F<T> eigensystem_F;

    EULER_1D(ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_input)
        :U(U_input)
    {}
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,1> >& psi_input)
    {psi_pointer=&psi_input;cut_out_grid=true;}

    void Set_Custom_Equation_Of_State(EOS<T>& eos_input)
    {eigensystem_F.Set_Custom_Equation_Of_State(eos_input);EULER<TV>::Set_Custom_Equation_Of_State(eos_input);}
    
    void Initialize_Domain(const int m, const T xmin, const T xmax)
    {grid.Initialize(m,xmin,xmax);}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
}
#endif
