//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER
//##################################################################### 
#ifndef __SHALLOW_WATER__
#define __SHALLOW_WATER__    

#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
namespace PhysBAM{

template<class TV>
class SHALLOW_WATER
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    enum{d=TV::m+1};
    typedef VECTOR<T,d> TV_DIMENSION;
public:
    BOUNDARY<TV,TV_DIMENSION>* boundary;
    CONSERVATION<TV,TV::m+1>* conservation;
private:
    BOUNDARY<TV,TV_DIMENSION> boundary_default;
    CONSERVATION_ENO_LLF<TV,TV::m+1> conservation_default;
public:

    T gravity;
    T min_height;
    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,TV_INT>& U; // h, h*u, and h*v
    VECTOR<EIGENSYSTEM<T,d>*,TV::m> eigensystems;

    SHALLOW_WATER(GRID<TV>& grid,ARRAY<TV_DIMENSION,TV_INT>& U,T gravity=9.8,
        T min_height=1e-3);
    ~SHALLOW_WATER();

//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
}    
#endif

