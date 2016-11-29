//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Avi Robinson-Mosher, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CELL_ITERATOR
//#####################################################################
#ifndef __CELL_ITERATOR__
#define __CELL_ITERATOR__

#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/GRID_ITERATOR_BASE.h>
namespace PhysBAM{

template<class TV>
class CELL_ITERATOR:public GRID_ITERATOR_BASE<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef typename GRID<TV>::REGION T_REGION;
    typedef TV VECTOR_T;
    using GRID_ITERATOR_BASE<TV>::grid;using GRID_ITERATOR_BASE<TV>::index;using GRID_ITERATOR_BASE<TV>::Add_Region;using GRID_ITERATOR_BASE<TV>::Reset;

    CELL_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells=0,const T_REGION& region_type=GRID<TV>::WHOLE_REGION,const int side=-1);

    CELL_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input)
        :GRID_ITERATOR_BASE<TV>(grid_input,region_input)
    {}

    const TV_INT& Cell_Index() const
    {return index;}

    TV Location() const
    {return grid.Center(index);}

    RANGE<TV> Bounding_Box() const
    {TV minimum_corner=grid.Node(index);return RANGE<TV>(minimum_corner,minimum_corner+grid.dX);}

    TV_INT Cell_Node_Index(const int node) const
    {return index+GRID<TV>::Binary_Counts(index)(node);}

    TV_INT First_Face_Index(const int axis) const
    {return grid.First_Face_Index_In_Cell(axis,index);}

    TV_INT Second_Face_Index(const int axis) const
    {return grid.Second_Face_Index_In_Cell(axis,index);}

    FACE_INDEX<TV::m> Full_First_Face_Index(const int axis) const
    {return FACE_INDEX<TV::m>(axis,grid.First_Face_Index_In_Cell(axis,index));}

    FACE_INDEX<TV::m> Full_Second_Face_Index(const int axis) const
    {return FACE_INDEX<TV::m>(axis,grid.Second_Face_Index_In_Cell(axis,index));}

    TV_INT Cell_Neighbor(const int i) const // 1 to 2*dimension
    {return GRID<TV>::Node_Neighbor(index,i);}
};
//#####################################################################
}
#endif
