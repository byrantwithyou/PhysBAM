//#####################################################################
// Copyright 2010, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace  SMOOTH_UNIFORM
//##################################################################### 
#ifndef __SMOOTH_UNIFORM__
#define __SMOOTH_UNIFORM__

#include <Grid_Tools/Grids/NODE_ITERATOR.h>

namespace PhysBAM{

namespace SMOOTH{

//#####################################################################
// Function Smooth
//#####################################################################
// explicit update with CFL=1/2 and all Neumann outer boundaries, no smoothing where phi is negative if defined
template<class T,int dim,class T_ARRAYS_T2>
void Smooth(T_ARRAYS_T2& d,const int steps,const ARRAY<T,VECTOR<int,dim> >* phi);
}
}
#endif
