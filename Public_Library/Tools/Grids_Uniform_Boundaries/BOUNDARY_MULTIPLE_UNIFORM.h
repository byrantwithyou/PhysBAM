//#####################################################################
// Copyright 2008, Jon Gretarsson, Nipun Kwatra
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MULTIPLE_UNIFORM 
//#####################################################################
#ifndef __BOUNDARY_MULTIPLE_UNIFORM__
#define __BOUNDARY_MULTIPLE_UNIFORM__

#include <Tools/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_MULTIPLE_UNIFORM:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_DIMENSION_T2;
    typedef VECTOR<BOUNDARY<TV,T2>*,2*TV::m> T_BOUNDARY_FACE_VECTOR;
    typedef BOUNDARY<TV,T2> BASE;

    T_BOUNDARY_FACE_VECTOR boundaries;

public:
    using BASE::Find_Ghost_Regions;
    BOUNDARY_MULTIPLE_UNIFORM(const T_BOUNDARY_FACE_VECTOR& boundaries_input)
        :boundaries(boundaries_input) {}

    ~BOUNDARY_MULTIPLE_UNIFORM()
    {for(int side=0;side<2*TV::m;side++) delete boundaries[side];}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_T2& u,T_ARRAYS_DIMENSION_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_DIMENSION_T2& u,const T time) const override;
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_MULTIPLE_UNIFORM<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const T_ARRAYS_DIMENSION_T2& u,T_ARRAYS_DIMENSION_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    T_ARRAYS_DIMENSION_T2::Put(u,u_ghost);
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=0;side<2*TV::m;side++) boundaries[side]->Fill_Single_Ghost_Region(grid,u_ghost,regions(side),side,dt,time,number_of_ghost_cells); 
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY_MULTIPLE_UNIFORM<TV,T2>::
Apply_Boundary_Condition(const GRID<TV>& grid,T_ARRAYS_DIMENSION_T2& u,const T time) const
{
    for(int side=0;side<2*TV::m;side++) boundaries[side]->Apply_Boundary_Condition(grid,u,time);
}
}
#endif
