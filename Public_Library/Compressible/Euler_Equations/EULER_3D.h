//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_3D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_3D class.
// Input U as 5 by (1,m) by (1,n) by (1,mn) for mass, momentum, and energy.
//
// Use Set_Up_Cut_Out_Grid(psi) to define a cut out grid with psi as (1,m) by (1,n) by (1,mn).
// When psi=1, solve the equaitions. 
// When psi=0, do NOT solve the equations.
//
//#####################################################################
#ifndef __EULER_3D__
#define __EULER_3D__

#include <Compressible/Euler_Equations/EULER.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
#include <Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{

template<class T_input>
class EULER_3D:public EULER<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,5> TV_DIMENSION;
protected:
    using EULER<TV>::boundary;using EULER<TV>::conservation;using EULER<TV>::eos;using EULER<TV>::e;using EULER<TV>::Set_Custom_Equation_Of_State;

    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,3> >& U;             // mass, momentum, and energy
    ARRAY<bool,VECTOR<int,3> >* psi_pointer; // defines cut out grid
    EULER_EIGENSYSTEM<TV> eigensystem_F;
    EULER_EIGENSYSTEM<TV> eigensystem_G;
    EULER_EIGENSYSTEM<TV> eigensystem_H;
    
public:
    EULER_3D(EOS<T>& eos_input,GRID<TV>& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,3> >& U_input)  
        :grid(grid_input),U(U_input),psi_pointer(0),eigensystem_F(&this->eos_default,0),eigensystem_G(&this->eos_default,1),eigensystem_H(&this->eos_default,2)
    {
        Set_Custom_Equation_Of_State(eos_input);
    }
    
    void Set_Up_Cut_Out_Grid(ARRAY<bool,VECTOR<int,3> >& psi_input)
    {psi_pointer=&psi_input;}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
}    
#endif
