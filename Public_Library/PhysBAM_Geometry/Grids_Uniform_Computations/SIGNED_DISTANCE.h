//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef __SIGNED_DISTANCE__
#define __SIGNED_DISTANCE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
namespace PhysBAM{
namespace SIGNED_DISTANCE{
//#####################################################################
template<class T,class TV,class TV_INT,class T_GEOMETRY> void Calculate(T_GEOMETRY& geometry,const GRID<TV>& grid,ARRAY<T,TV_INT>& phi,bool print_progress=false){
    for(UNIFORM_GRID_ITERATOR_NODE<TV> iterator(grid,grid.Domain_Indices());iterator.Valid();iterator.Next()){TV_INT index=iterator.Node_Index();
        phi(index)=geometry.Signed_Distance(grid.X(index));}
}
//#####################################################################
};
};
#endif
