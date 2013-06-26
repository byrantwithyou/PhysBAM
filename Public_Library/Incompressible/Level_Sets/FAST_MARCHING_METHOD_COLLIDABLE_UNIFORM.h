//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM  
//#####################################################################
#ifndef __FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM__
#define __FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Incompressible/Level_Sets/LEVELSET_COLLIDABLE.h>
namespace PhysBAM{

template<class T_GRID>
class FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM:public FAST_MARCHING_METHOD_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef FAST_MARCHING_METHOD_UNIFORM<T_GRID> BASE;
public:
    using BASE::Neighbor_Visible;using BASE::cell_grid;using BASE::dimension_start;using BASE::dimension_end;
    const LEVELSET_COLLIDABLE<TV>& levelset_collidable;

    FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM(const LEVELSET_COLLIDABLE<TV>& levelset_collidable,const int ghost_cells,THREAD_QUEUE* thread_queue_input=0);
    ~FAST_MARCHING_METHOD_COLLIDABLE_UNIFORM();

    void Initialize_Interface_Threaded(RANGE<TV_INT>& domain,ARRAY<T,TV_INT>& phi_ghost,ARRAY<T,TV_INT>& phi_new,ARRAY<bool,TV_INT>& done) PHYSBAM_OVERRIDE;

//#####################################################################
};
}
#endif
