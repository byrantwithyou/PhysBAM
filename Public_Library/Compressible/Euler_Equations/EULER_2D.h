//#####################################################################
// Copyright 2002, 2003, Doug Enright, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_2D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_2D class.
// Input U as 4 by (0,m) by (0,n) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (0,m) by (0,n).
// When psi=true, solve the equaitions. 
// When psi=false, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_2D__
#define __EULER_2D__

#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_F.h>
#include <Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_G.h>
namespace PhysBAM{

template<class T_input>
class EULER_2D:public EULER<VECTOR<T_input,2> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,4> TV_DIMENSION;
protected:
    using EULER<TV>::cut_out_grid;using EULER<TV>::max_time_step;
public:
    using EULER<TV>::boundary;using EULER<TV>::conservation;using EULER<TV>::eos;

    GRID<TV> grid;
    ARRAY<TV_DIMENSION,VECTOR<int,2> >& U; // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,2> >* psi_pointer; // defines cut out grid
    EULER_2D_EIGENSYSTEM_F<T> eigensystem_F;
    EULER_2D_EIGENSYSTEM_G<T> eigensystem_G;

    EULER_2D(ARRAY<TV_DIMENSION,VECTOR<int,2> >& U_input)
        :U(U_input)
    {}
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,2> >& psi_input)
    {psi_pointer=&psi_input;cut_out_grid=true;}

    void Set_Custom_Equation_Of_State(EOS<T>& eos_input)
    {eigensystem_F.Set_Custom_Equation_Of_State(eos_input);eigensystem_G.Set_Custom_Equation_Of_State(eos_input);EULER<TV>::Set_Custom_Equation_Of_State(eos_input);}
    
    void Initialize_Domain(const int m, const int n, const T xmin, const T xmax, const T ymin, const T ymax)
    {grid.Initialize(TV_INT(m,n),RANGE<TV>(TV(xmin,ymin),TV(xmax,ymax)));}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
}
#endif
