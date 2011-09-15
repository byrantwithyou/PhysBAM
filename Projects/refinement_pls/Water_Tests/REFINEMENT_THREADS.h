//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __REFINEMENT_THREADS__
#define __REFINEMENT_THREADS__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class TV> class WATER_TESTS;

template<class TV>
class REFINEMENT_TASK:public THREAD_QUEUE::TASK
{
    typedef typename TV::SCALAR T;    
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    
public:
    WATER_TESTS<TV>* smoke_tests;
    TV_INT cell_index;
    T time,dt;

    void Run();
};
}
#endif
