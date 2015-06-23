//#####################################################################
// Copyright 2005-2008, Eran Guendelman, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EDGE_ITERATOR
//#####################################################################
#ifndef __EDGE_ITERATOR__
#define __EDGE_ITERATOR__

#include <Tools/Grids_Uniform/EDGE_INDEX.h>
#include <Tools/Grids_Uniform/GRID_ITERATOR_BASE.h>
#include <Tools/Utilities/PHYSBAM_ATTRIBUTE.h>

namespace PhysBAM{

template<class TV>
class EDGE_ITERATOR:public GRID_ITERATOR_BASE<TV>
{
public:
    typedef typename GRID<TV>::REGION T_REGION;typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    using GRID_ITERATOR_BASE<TV>::grid;using GRID_ITERATOR_BASE<TV>::index;using GRID_ITERATOR_BASE<TV>::region;using GRID_ITERATOR_BASE<TV>::valid;
    using GRID_ITERATOR_BASE<TV>::Reset;using GRID_ITERATOR_BASE<TV>::current_region;using GRID_ITERATOR_BASE<TV>::Add_Region;
    using GRID_ITERATOR_BASE<TV>::Reset_Regions;

protected:
    T_REGION region_type;
    int side;
    int axis;
    bool single_axis;
    int number_of_ghost_cells;

public:
    // axis_input==0 means iterate through faces in all dimensions
    EDGE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells_input=0,const T_REGION& region_type_input=GRID<TV>::WHOLE_REGION,const int side_input=-1,
        int axis_input=-1);

    EDGE_ITERATOR(const GRID<TV>& grid_input,const int axis_input,const TV_INT& face_index);

    EDGE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& explicit_region_input,const int axis_input);

private:
    void Reset_Axis(const int axis_input);
    void Next_Helper();

public:
    void Next() PHYSBAM_ALWAYS_INLINE // overloads GRID_ITERATOR_BASE::Next but we don't want that to be virtual to avoid virtual call overhead
    {if(index(TV::dimension-1)<region.max_corner(TV::dimension-1)-1) index(TV::dimension-1)++;else Next_Helper();}

    int Axis() const
    {return axis;}

    EDGE_INDEX<TV::dimension> Full_Index() const
    {return EDGE_INDEX<TV::dimension>(axis,index);}
//#####################################################################
};
}
#endif
