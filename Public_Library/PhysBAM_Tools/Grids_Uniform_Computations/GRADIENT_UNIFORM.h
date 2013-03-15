//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace  GRADIENT_UNIFORM
//##################################################################### 
#ifndef __GRADIENT_UNIFORM__
#define __GRADIENT_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

namespace GRADIENT{
template<class T,class TV>
void Compute_Magnitude(const GRID<TV>& grid,const int number_of_ghost_cells,ARRAY<T,VECTOR<int,TV::m> >& values,ARRAY<T,VECTOR<int,TV::m> >& gradient)
{
    for(CELL_ITERATOR<TV> iterator(grid,number_of_ghost_cells-1);iterator.Valid();iterator.Next()){
        VECTOR<int,TV::m> cell_index=iterator.Cell_Index();T sum_of_partials=0;
        for(int axis=0;axis<TV::m;axis++){
            T partial=(values(iterator.Cell_Neighbor(2*axis+1))-values(iterator.Cell_Neighbor(2*axis)))*(T).5*grid.one_over_dX[axis];
            sum_of_partials+=partial*partial;}
        gradient(cell_index)=sqrt(sum_of_partials);}
    for(CELL_ITERATOR<TV> iterator(grid,number_of_ghost_cells,GRID<TV>::BOUNDARY_INTERIOR_REGION);iterator.Valid();iterator.Next()) gradient(iterator.Cell_Index())=(T)0;
}
}
}
#endif
