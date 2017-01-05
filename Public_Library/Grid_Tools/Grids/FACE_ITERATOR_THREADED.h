//#####################################################################
// Copyright 2017, Ounan Ding, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ITERATOR_THREADED
//#####################################################################
#ifndef __FACE_ITERATOR_THREADED__
#define __FACE_ITERATOR_THREADED__

#include <Core/Utilities/PHYSBAM_ATTRIBUTE.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/GRID_ITERATOR_BASE.h>

namespace PhysBAM{

template<class TV>
class FACE_ITERATOR_THREADED:public GRID_ITERATOR_BASE<TV>
{
public:
    typedef typename GRID<TV>::REGION T_REGION;typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    typedef TV VECTOR_T;
    using GRID_ITERATOR_BASE<TV>::grid;using GRID_ITERATOR_BASE<TV>::index;using GRID_ITERATOR_BASE<TV>::region;using GRID_ITERATOR_BASE<TV>::valid;
    using GRID_ITERATOR_BASE<TV>::Reset;using GRID_ITERATOR_BASE<TV>::current_region;using GRID_ITERATOR_BASE<TV>::Add_Region;
    using GRID_ITERATOR_BASE<TV>::Reset_Regions;using GRID_ITERATOR_BASE<TV>::number_of_regions;

    int axis;
protected:
    int axes[TV::m];
    
public:
    // axis_input==0 means iterate through faces in all dimensions
    FACE_ITERATOR_THREADED(const GRID<TV>& grid_input,int number_of_ghost_cells_input=0,T_REGION region_type_input=GRID<TV>::WHOLE_REGION);

private:
    void Next_Helper()
    {GRID_ITERATOR_BASE<TV>::Next_Helper();if(valid) axis=axes[current_region];}

public:
    void Next() PHYSBAM_ALWAYS_INLINE // overloads GRID_ITERATOR_BASE::Next but we don't want that to be virtual to avoid virtual call overhead
    {if(index(TV::m-1)<region.max_corner(TV::m-1)-1) index(TV::m-1)++;else Next_Helper();}

    int Axis() const
    {return axis;}

    const TV_INT& Face_Index() const
    {return index;}

    FACE_INDEX<TV::m> Full_Index() const
    {return FACE_INDEX<TV::m>(axis,index);}

    TV Location() const
    {return grid.Face(Full_Index());}

    TV_INT First_Cell_Index() const
    {TV_INT i(index);i(axis)--;return i;}

    TV_INT Second_Cell_Index() const
    {return index;}

    RANGE<TV> Dual_Cell() const
    {return RANGE<TV>(Location()-(T).5*grid.dX,Location()+(T).5*grid.dX);}

    TV First_Cell_Center() const
    {return grid.Center(First_Cell_Index());}

    TV Second_Cell_Center() const
    {return grid.Center(Second_Cell_Index());}

    TV_INT Face_Node_Index(const int node) const // 1-based
    {return grid.Face_Node_Index(axis,index,node);}
//#####################################################################
};
}
#endif
