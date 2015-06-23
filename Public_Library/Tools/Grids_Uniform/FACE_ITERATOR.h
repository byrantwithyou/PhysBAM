//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ITERATOR
//#####################################################################
#ifndef __FACE_ITERATOR__
#define __FACE_ITERATOR__

#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Grids_Uniform/GRID_ITERATOR_BASE.h>
#include <Tools/Utilities/PHYSBAM_ATTRIBUTE.h>

namespace PhysBAM{

template<class TV>
class FACE_ITERATOR:public GRID_ITERATOR_BASE<TV>
{
public:
    typedef typename GRID<TV>::REGION T_REGION;typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    typedef TV VECTOR_T;
    using GRID_ITERATOR_BASE<TV>::grid;using GRID_ITERATOR_BASE<TV>::index;using GRID_ITERATOR_BASE<TV>::region;using GRID_ITERATOR_BASE<TV>::valid;
    using GRID_ITERATOR_BASE<TV>::Reset;using GRID_ITERATOR_BASE<TV>::current_region;using GRID_ITERATOR_BASE<TV>::Add_Region;
    using GRID_ITERATOR_BASE<TV>::Reset_Regions;

    int axis;
protected:
    T_REGION region_type;
    int side;
    bool single_axis;
    int number_of_ghost_cells;
    T face_size;

public:
    // axis_input==0 means iterate through faces in all dimensions
    FACE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells_input=0,const T_REGION& region_type_input=GRID<TV>::WHOLE_REGION,const int side_input=-1,
        int axis_input=-1);

    FACE_ITERATOR(const GRID<TV>& grid_input,const int axis_input,const TV_INT& face_index);

    FACE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,const int axis_input);

private:
    void Reset_Axis(const int axis_input);
    void Next_Helper();

public:
    void Next() PHYSBAM_ALWAYS_INLINE // overloads GRID_ITERATOR_BASE::Next but we don't want that to be virtual to avoid virtual call overhead
    {if(index(TV::dimension-1)<region.max_corner(TV::dimension-1)-1) index(TV::dimension-1)++;else Next_Helper();}

    int Axis() const
    {return axis;}

    const TV_INT& Face_Index() const
    {return index;}

    FACE_INDEX<TV::dimension> Full_Index() const
    {return FACE_INDEX<TV::dimension>(axis,index);}

    T Face_Size() const
    {return face_size;}

    TV Location() const
    {return grid.Face(Full_Index());}

    TV_INT First_Cell_Index() const
    {TV_INT i(index);i(axis)--;return i;}

    TV_INT Second_Cell_Index() const
    {return index;}

    void Unordered_Cell_Indices_Touching_Face(TV_INT& cell1,TV_INT& cell2) const
    {cell1=First_Cell_Index();cell2=Second_Cell_Index();}

    RANGE<TV> Dual_Cell() const
    {return RANGE<TV>(Location()-(T).5*grid.dX,Location()+(T).5*grid.dX);}

    TV First_Cell_Center() const
    {return grid.Center(First_Cell_Index());}

    TV Second_Cell_Center() const
    {return grid.Center(Second_Cell_Index());}

    bool First_Boundary() const // returns true if currently on left, bottom, or front boundary
    {assert(region_type==GRID<TV>::BOUNDARY_REGION);return (side<0 && current_region%2==0) || (side>=0 && side%2==0);}

    TV_INT Face_Node_Index(const int node) const // 1-based
    {return grid.Face_Node_Index(axis,index,node);}
//#####################################################################
};
}
#endif
