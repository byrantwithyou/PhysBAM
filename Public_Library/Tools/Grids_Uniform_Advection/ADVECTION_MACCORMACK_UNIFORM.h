//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_MACCORMACK_UNIFORM
//#####################################################################
#ifndef __ADVECTION_MACCORMACK_UNIFORM__
#define __ADVECTION_MACCORMACK_UNIFORM__

#include <Tools/Advection/ADVECTION.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class TV,class T2,class T_NESTED_ADVECTION>
class ADVECTION_MACCORMACK_UNIFORM:public ADVECTION<TV,T2>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef ARRAY<bool,FACE_INDEX<TV::m> > T_FACE_ARRAYS_BOOL;
private:
    // false to use 1st order advection
    const ARRAY<bool,TV_INT>* node_mask;
    const ARRAY<bool,TV_INT>* cell_mask;
    const T_FACE_ARRAYS_BOOL* face_mask;
    T_NESTED_ADVECTION& nested_advection;
public:
    bool clamp_extrema; // otherwise fall back to 1st order
    bool ensure_second_order; // sacrifice stability for accuracy.
    THREAD_QUEUE* thread_queue;

    ADVECTION_MACCORMACK_UNIFORM(T_NESTED_ADVECTION& nested_advection_input,const ARRAY<bool,TV_INT>* node_mask_input,const ARRAY<bool,TV_INT>* cell_mask_input,const T_FACE_ARRAYS_BOOL* face_mask_input,THREAD_QUEUE* thread_queue=0);

    void Update_Advection_Equation_Node(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost_input=0,const ARRAY<T2,TV_INT>* Z_max_ghost_input=0,ARRAY<T2,TV_INT>* Z_min_input=0,ARRAY<T2,TV_INT>* Z_max_input=0);
    void Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,const FACE_LOOKUP_UNIFORM<TV>& face_velocities,
        BOUNDARY<TV,T2>& boundary,const T dt,const T time,const ARRAY<T2,TV_INT>* Z_min_ghost_input,const ARRAY<T2,TV_INT>* Z_max_ghost_input,ARRAY<T2,TV_INT>* Z_min_input,
        ARRAY<T2,TV_INT>* Z_max_input);
    void Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& Z,const FACE_LOOKUP_UNIFORM<TV>& Z_ghost,
        const FACE_LOOKUP_UNIFORM<TV>& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,const FACE_LOOKUP_UNIFORM<TV>* Z_min_ghost_input,
        const FACE_LOOKUP_UNIFORM<TV>* Z_max_ghost_input,T_FACE_ARRAYS_SCALAR* Z_min_input,T_FACE_ARRAYS_SCALAR* Z_max_input);
private:
    void Apply_Clamped_Extrema_Limiter_Node_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_forward_ghost,const ARRAY<T2,TV_INT>& Z_backward_ghost,const ARRAY<T2,TV_INT>& Z_min,const ARRAY<T2,TV_INT>& Z_max);
    void Apply_Clamped_Extrema_Limiter_Cell_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_forward_ghost,const ARRAY<T2,TV_INT>& Z_backward_ghost,const ARRAY<T2,TV_INT>& Z_min,const ARRAY<T2,TV_INT>& Z_max);
    void Apply_Clamped_Extrema_Limiter_Face_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost,const T_FACE_ARRAYS_SCALAR& Z_min,const T_FACE_ARRAYS_SCALAR& Z_max);

    void Apply_Reversion_Limiter_Node_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_forward_ghost,const ARRAY<T2,TV_INT>& Z_backward_ghost,const ARRAY<T2,TV_INT>& Z_min,const ARRAY<T2,TV_INT>& Z_max);
    void Apply_Reversion_Limiter_Cell_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_forward_ghost,const ARRAY<T2,TV_INT>& Z_backward_ghost,const ARRAY<T2,TV_INT>& Z_min,const ARRAY<T2,TV_INT>& Z_max);
    void Apply_Reversion_Limiter_Face_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost,const T_FACE_ARRAYS_SCALAR& Z_min,const T_FACE_ARRAYS_SCALAR& Z_max);

    void Apply_Second_Order_Update_Face_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,int axis,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_forward_ghost,const T_FACE_ARRAYS_SCALAR& Z_backward_ghost);

//#####################################################################
};
}
#endif
