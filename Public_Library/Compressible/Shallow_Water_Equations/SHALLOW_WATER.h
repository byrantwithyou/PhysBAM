//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER
//##################################################################### 
//
// Inherited by SHALLOW_WATER_1D, SHALLOW_WATER_2D, and SHALLOW_WATER_3D.
//
//#####################################################################
#ifndef __SHALLOW_WATER__
#define __SHALLOW_WATER__    

#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
namespace PhysBAM{

template<class TV>
class SHALLOW_WATER
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,TV::m+1> TV_DIMENSION;
public:
    BOUNDARY<TV,TV_DIMENSION>* boundary;
    CONSERVATION<TV,TV::m+1>* conservation;
private:
    BOUNDARY<TV,TV_DIMENSION> boundary_default;
    CONSERVATION_ENO_LLF<TV,TV::m+1> conservation_default;

protected:
    SHALLOW_WATER()
    {
        boundary=&boundary_default;
        conservation=&conservation_default;
    }
public:

    void Set_Custom_Boundary(BOUNDARY<TV,TV_DIMENSION>& boundary_input)
    {boundary=&boundary_input;}
    
    void Set_Custom_Conservation(CONSERVATION<T,TV::m+1>& conservation_input)
    {conservation=&conservation_input;}

//#####################################################################
};
}    
#endif

