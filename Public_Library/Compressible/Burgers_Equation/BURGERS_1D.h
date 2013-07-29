//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BURGERS_1D  
//##################################################################### 
//
// Input U as 1 by (1,m).
//
//#####################################################################
#ifndef __BURGERS_1D__
#define __BURGERS_1D__

#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Compressible/Burgers_Equation/BURGERS_1D_EIGENSYSTEM_F.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <cfloat>
namespace PhysBAM{

template<class T>
class BURGERS_1D
{
    typedef VECTOR<T,1> TV;
public:
    BOUNDARY<TV,TV>* boundary;
    CONSERVATION<TV,1>* conservation;
protected:
    GRID<TV>& grid;
    ARRAY<TV,VECTOR<int,1> >& U;
    BOUNDARY<TV,TV> boundary_default;
    CONSERVATION_ENO_LLF<TV,1> conservation_default;
    BURGERS_1D_EIGENSYSTEM_F<T> eigensystem_F;
private:
     T max_time_step;

public:
    BURGERS_1D(GRID<TV>& grid_input,ARRAY<TV,VECTOR<int,1> >& U_input)
        :grid(grid_input),U(U_input)
    {
        boundary=&boundary_default;
        conservation=&conservation_default;
        Set_Max_Time_Step();
    }
    
    void Set_Custom_Boundary(BOUNDARY<TV,TV>& boundary_input)
    {boundary=&boundary_input;}
    
    void Set_Custom_Conservation(CONSERVATION<TV,1>& conservation_input)
    {conservation=&conservation_input;}
    
    void Set_Max_Time_Step(const T max_time_step_input=FLT_MAX)
    {max_time_step=max_time_step_input;}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};   
}
#endif
