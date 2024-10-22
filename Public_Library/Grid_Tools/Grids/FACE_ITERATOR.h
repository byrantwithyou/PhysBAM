//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_ITERATOR
//#####################################################################
#ifndef __FACE_ITERATOR__
#define __FACE_ITERATOR__

#include <Core/Utilities/PHYSBAM_ATTRIBUTE.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>

namespace PhysBAM{

template<class TV>
class FACE_ITERATOR:public FACE_RANGE_ITERATOR<TV::m>
{
public:
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    typedef TV VECTOR_T;
    typedef FACE_RANGE_ITERATOR<TV::m> BASE;
    using BASE::Set_Range;using BASE::Initialize;using BASE::face;
    using BASE::side;

    const GRID<TV>& grid;

    // axis_input==0 means iterate through faces in all dimensions
    FACE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells_input=0,
        const typename GRID<TV>::REGION& region_type_input=GRID<TV>::WHOLE_REGION,
        const int side_input=-1,int axis_input=-1);

    FACE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,
        const int axis_input);

    int Axis() const
    {return face.axis;}

    const TV_INT& Face_Index() const
    {return face.index;}

    FACE_INDEX<TV::m> Full_Index() const
    {return face;}

    TV Location() const
    {return grid.Face(Full_Index());}

    TV_INT First_Cell_Index() const
    {TV_INT i(face.index);i(face.axis)--;return i;}

    TV_INT Second_Cell_Index() const
    {return face.index;}

    void Unordered_Cell_Indices_Touching_Face(TV_INT& cell1,TV_INT& cell2) const
    {cell1=First_Cell_Index();cell2=Second_Cell_Index();}

    RANGE<TV> Dual_Cell() const
    {return RANGE<TV>(Location()-(T).5*grid.dX,Location()+(T).5*grid.dX);}

    TV First_Cell_Center() const
    {return grid.Center(First_Cell_Index());}

    TV Second_Cell_Center() const
    {return grid.Center(Second_Cell_Index());}

    bool First_Boundary() const // returns true if currently on left, bottom, or front boundary
    {return side%2==0;}

    TV_INT Face_Node_Index(const int node) const // 1-based
    {return grid.Face_Node_Index(face.axis,face.index,node);}
//#####################################################################
};
}
#endif
