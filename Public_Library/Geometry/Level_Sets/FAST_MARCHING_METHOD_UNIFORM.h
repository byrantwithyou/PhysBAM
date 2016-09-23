//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING_METHOD_UNIFORM  
//#####################################################################
#ifndef __FAST_MARCHING_METHOD_UNIFORM__
#define __FAST_MARCHING_METHOD_UNIFORM__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <functional>
namespace PhysBAM{

template<class TV>
class FAST_MARCHING_METHOD_UNIFORM:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    const LEVELSET<TV>& levelset;
protected:
    GRID<TV> cell_grid;
    TV_INT dimension_start,dimension_end;
    int ghost_cells;
public:
    THREAD_QUEUE* thread_queue;
    std::function<bool(const int axis,const TV_INT& current_index)> Neighbor_Visible;

    FAST_MARCHING_METHOD_UNIFORM(const LEVELSET<TV>& levelset,const int ghost_cells,THREAD_QUEUE* thread_queue_input=0);
    ~FAST_MARCHING_METHOD_UNIFORM();

//#####################################################################
    void Fast_Marching_Method_Threaded(RANGE<TV_INT>& domain,ARRAY<T,TV_INT>& phi_ghost,const T stopping_distance,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Fast_Marching_Method(ARRAY<T,TV_INT>& phi_ghost,const T stopping_distance=0,const ARRAY<TV_INT>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    virtual void Initialize_Interface_Threaded(RANGE<TV_INT>& domain,ARRAY<T,TV_INT>& phi_ghost,ARRAY<T,TV_INT>& phi_new,ARRAY<bool,TV_INT>& done);
private:
    void Update_Or_Add_Neighbor(ARRAY<T,TV_INT>& phi_ghost,ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,ARRAY<TV_INT>& heap,int& heap_length,const TV_INT& neighbor);
    void Initialize_Interface(RANGE<TV_INT>& domain,ARRAY<T,TV_INT>& phi_ghost,ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices,const bool add_seed_indices_for_ghost_cells=false);
    void Initialize_Interface(ARRAY<T,TV_INT>& phi_ghost,ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,ARRAY<TV_INT>& heap,int& heap_length,const ARRAY<TV_INT>* seed_indices=0,
        const bool add_seed_indices_for_ghost_cells=false);
    void Initialize_Interface(ARRAY<T,TV_INT>& phi_ghost,ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,ARRAY<TV_INT>& heap,int& heap_length,const bool add_seed_indices_for_ghost_cells=false);
    void Update_Close_Point(ARRAY<T,TV_INT>& phi_ghost,const ARRAY<bool,TV_INT>& done,const TV_INT& index);
    void Add_To_Initial(ARRAY<bool,TV_INT>& done,ARRAY<int,TV_INT>& close_k,const TV_INT& index);
//#####################################################################
};
}
#endif
