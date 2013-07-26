//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __REFINEMENT_THREADS__
#define __REFINEMENT_THREADS__
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class TV> class PROJECTION_REFINEMENT_UNIFORM;

template<class TV>
class REFINEMENT_TASK:public THREAD_QUEUE::TASK
{
    typedef typename TV::SCALAR T;    
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    
public:
    PROJECTION_REFINEMENT_UNIFORM<TV>* projection;
    ARRAY<T,FACE_INDEX<TV::dimension> > *fine_face_velocities;
    TV_INT cell_index;
    T time,dt;

    void Run();
};
}
#endif
