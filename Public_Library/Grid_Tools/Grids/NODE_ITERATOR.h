//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NODE_ITERATOR
//#####################################################################
#ifndef __NODE_ITERATOR__
#define __NODE_ITERATOR__

#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/GRID_ITERATOR_BASE.h>
namespace PhysBAM{

template<class TV>
class NODE_ITERATOR:public GRID_ITERATOR_BASE<TV>
{
public:
    typedef typename GRID<TV>::REGION T_REGION;typedef VECTOR<int,TV::m> TV_INT;
    typedef TV VECTOR_T;
    using GRID_ITERATOR_BASE<TV>::grid;using GRID_ITERATOR_BASE<TV>::index;using GRID_ITERATOR_BASE<TV>::Add_Region;using GRID_ITERATOR_BASE<TV>::Reset;

    NODE_ITERATOR(const GRID<TV>& grid_input,const int number_of_ghost_cells=0,const T_REGION& region_type=GRID<TV>::WHOLE_REGION,const int side=-1);
    NODE_ITERATOR(const GRID<TV>& grid_input,const RANGE<TV_INT>& region_input);

    const TV_INT& Node_Index() const
    {return index;}

    TV Location() const
    {return grid.Node(index);}

    TV_INT Node_Neighbor(const int i) const // i=1 to 6
    {return grid.Node_Neighbor(index,i);}

//#####################################################################
};
}
#endif
