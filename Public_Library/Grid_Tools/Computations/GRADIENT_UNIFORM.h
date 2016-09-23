//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace  GRADIENT_UNIFORM
//##################################################################### 
#ifndef __GRADIENT_UNIFORM__
#define __GRADIENT_UNIFORM__

#include <Core/Vectors/VECTOR_FORWARD.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
namespace PhysBAM{

namespace GRADIENT{
template<class T,class TV>
void Compute_Magnitude(const GRID<TV>& grid,const int number_of_ghost_cells,ARRAY<T,VECTOR<int,TV::m> >& values,
    ARRAY<T,VECTOR<int,TV::m> >& gradient);
}
}
#endif
