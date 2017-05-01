//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL
//#####################################################################
//
// Assumes a MAC grid straddling the centerline on the left, applying reflection.
// Extrapolation is applied to the right, top and bottom.
//
//#####################################################################
#ifndef __BOUNDARY_EULER_EQUATIONS_CYLINDRICAL__
#define __BOUNDARY_EULER_EQUATIONS_CYLINDRICAL__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class T_input>
class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL:public BOUNDARY<VECTOR<T_input,2>,VECTOR<T_input,4> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,4> TV_DIMENSION;typedef VECTOR<int,2> TV_INT;
    enum {d=4};
public:
    BOUNDARY_EULER_EQUATIONS_CYLINDRICAL() = default;
    ~BOUNDARY_EULER_EQUATIONS_CYLINDRICAL() = default;

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u,ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<TV_DIMENSION,VECTOR<int,2> >& u,const T time) const override {} // do nothing
//#####################################################################
};
}
#endif
