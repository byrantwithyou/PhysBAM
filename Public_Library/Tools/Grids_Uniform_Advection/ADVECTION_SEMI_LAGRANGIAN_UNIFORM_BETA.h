//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA__
#define __ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA__

#include <Tools/Advection/ADVECTION.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_UNIFORM_FORWARD.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Parallel_Computation/THREAD_QUEUE.h>
namespace PhysBAM{

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> //  T_AVERAGING=AVERAGING_UNIFORM<TV>, T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<TV,T2>
class ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA:public ADVECTION<TV,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
public:
    using ADVECTION<TV,T2,typename T_AVERAGING::FACE_LOOKUP>::Update_Advection_Equation_Cell;
    template<class T3> struct REBIND{typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T3,T_AVERAGING,typename T_INTERPOLATION::template REBIND<T3>::TYPE> TYPE;};
    template<class T_INTERPOLATION_2> struct REBIND_INTERPOLATION{typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION_2> TYPE;};

    THREAD_QUEUE* thread_queue;

//#####################################################################
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA(THREAD_QUEUE* thread_queue=0);
    void Update_Advection_Equation_Node(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max);
    void Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max);
    void Update_Advection_Equation_Node_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Cell_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Cell_Lookup_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max);
    void Update_Advection_Equation_Face_Lookup_Threaded(RANGE<TV_INT>& domain,int axis,const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max);
//#####################################################################
};
}
#endif
