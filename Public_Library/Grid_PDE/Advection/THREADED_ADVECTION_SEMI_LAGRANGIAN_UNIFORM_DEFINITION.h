//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Advection/THREADED_ADVECTION_SEMI_LAGRANGIAN_TASK.h>
#include <Grid_PDE/Advection/THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>

using namespace PhysBAM;

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T2,T_AVERAGING,T_INTERPOLATION>::
THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM()
    :thread_queue(0),row_jump(1)
{}

#ifdef USE_PTHREADS
template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Node(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
    const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    RANGE<TV_INT> domain(grid.Domain_Indices());
    int min_value=domain.min_corner.x,max_value=domain.max_corner.x;
    for(int i=min_value;i<max_value;i+=row_jump){
        domain.min_corner.x=min(i,max_value);domain.max_corner.x=min(i+row_jump-1,max_value);
        ADVECTION_SEMI_LAGRANGIAN_TASK_NODE<TV,T2,T_AVERAGING,T_INTERPOLATION>* task=
            new ADVECTION_SEMI_LAGRANGIAN_TASK_NODE<TV,T2,T_AVERAGING,T_INTERPOLATION>(grid,Z,Z_ghost,V,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max,domain);
        thread_queue->Queue(task);}
    thread_queue->Wait();
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    RANGE<TV_INT> domain(grid.Domain_Indices());
    int min_value=domain.min_corner.x,max_value=domain.max_corner.x;
    for(int i=min_value;i<max_value;i+=row_jump){
        domain.min_corner.x=min(i,max_value);domain.max_corner.x=min(i+row_jump-1,max_value);
        ADVECTION_SEMI_LAGRANGIAN_TASK_CELL<TV,T2,T_AVERAGING,T_INTERPOLATION>* task=
            new ADVECTION_SEMI_LAGRANGIAN_TASK_CELL<TV,T2,T_AVERAGING,T_INTERPOLATION>(grid,Z,Z_ghost,face_velocities,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max,domain);
        thread_queue->Queue(task);}
    thread_queue->Wait();
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T,FACE_INDEX<TV::m> >* Z_min,ARRAY<T,FACE_INDEX<TV::m> >* Z_max)
{
    for(int i=0;i<TV::dimension;i++){
        RANGE<TV_INT> domain(grid.Domain_Indices());
        int min_value=domain.min_corner(i),max_value=domain.max_corner(i)+1;
        for(int j=min_value;j<max_value;j+=row_jump){
            domain.min_corner(i)=j;domain.max_corner(i)=min(j+row_jump-1,max_value);
            ADVECTION_SEMI_LAGRANGIAN_TASK_FACE<TV,T2,T_AVERAGING,T_INTERPOLATION>* task=
                new ADVECTION_SEMI_LAGRANGIAN_TASK_FACE<TV,T2,T_AVERAGING,T_INTERPOLATION>(grid,Z,Z_ghost,face_velocities,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max,domain,i);
            thread_queue->Queue(task);}}
    thread_queue->Wait();
}
#endif
