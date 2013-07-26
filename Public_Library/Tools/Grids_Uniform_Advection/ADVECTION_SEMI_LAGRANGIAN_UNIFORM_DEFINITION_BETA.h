//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
#include <Tools/Vectors/VECTOR_3D.h>

using namespace PhysBAM;

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA(THREAD_QUEUE* thread_queue_input)
    :thread_queue(thread_queue_input)
{}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Node(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
    const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::All_Ones_Vector();
    DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>,TV>(domain,thread_queue).template Run<const GRID<TV>&,ARRAY<T2,TV_INT>&,const ARRAY<T2,TV_INT>&,const ARRAY<TV,TV_INT>&,BOUNDARY<TV,T2>&,T,T,const ARRAY<T2,TV_INT>*,const ARRAY<T2,TV_INT>*,ARRAY<T2,TV_INT>*,ARRAY<T2,TV_INT>*>(*this,&ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::Update_Advection_Equation_Node_Threaded,grid,Z,Z_ghost,V,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>,TV>(grid.Domain_Indices(),thread_queue).template Run<const GRID<TV>&,ARRAY<T2,TV_INT>&,const ARRAY<T2,TV_INT>&,const T_FACE_LOOKUP&,BOUNDARY<TV,T2>&,T,T,const ARRAY<T2,TV_INT>*,const ARRAY<T2,TV_INT>*,ARRAY<T2,TV_INT>*,ARRAY<T2,TV_INT>*>(*this,&ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::Update_Advection_Equation_Cell_Lookup_Threaded,grid,Z,Z_ghost,face_velocities,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    for(int axis=0;axis<TV::dimension;axis++){
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(axis);
        DOMAIN_ITERATOR_THREADED_ALPHA<ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>,TV>(domain,thread_queue).template Run<int,const GRID<TV>&,T_FACE_ARRAYS_SCALAR&,const T_FACE_LOOKUP&,const T_FACE_LOOKUP&,BOUNDARY<TV,T>&,T,T,const T_FACE_LOOKUP*,const T_FACE_LOOKUP*,T_FACE_ARRAYS_SCALAR*,T_FACE_ARRAYS_SCALAR*>(*this,&ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::Update_Advection_Equation_Face_Lookup_Threaded,axis,grid,Z,Z_ghost,face_velocities,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Node_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
    const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    T_INTERPOLATION interpolation;
    if(Z_min && Z_max) for(NODE_ITERATOR<TV> iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT node=iterator.Node_Index();TV X=iterator.Location()-dt*V(node);
        Z(node)=interpolation.Clamped_To_Array_Node(iterator.grid,Z_ghost,X);
        VECTOR<T2,2> extrema=interpolation.Extrema_Clamped_To_Array_Node(iterator.grid,*Z_min_ghost,*Z_max_ghost,X);
        (*Z_min)(node)=extrema.x;(*Z_max)(node)=extrema.y;}
    else for(NODE_ITERATOR<TV> iterator(grid,domain);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
        Z(node)=interpolation.Clamped_To_Array(iterator.grid,Z_ghost,iterator.Location()-dt*V(node));}
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
    const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    T_INTERPOLATION interpolation;
    if(Z_min && Z_max) for(CELL_ITERATOR<TV> iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();TV X=iterator.Location()-dt*V(cell);
        Z(cell)=interpolation.Clamped_To_Array(iterator.grid,Z_ghost,X);
        VECTOR<T2,2> extrema=interpolation.Extrema_Clamped_To_Array_Cell(iterator.grid,*Z_min_ghost,*Z_max_ghost,X);
        (*Z_min)(cell)=extrema.x;(*Z_max)(cell)=extrema.y;}
    else for(CELL_ITERATOR<TV> iterator(grid,domain);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        Z(cell)=interpolation.Clamped_To_Array(iterator.grid,Z_ghost,iterator.Location()-dt*V(cell));}
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell_Lookup_Threaded(RANGE<TV_INT>& domain,const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    T_INTERPOLATION interpolation;T_AVERAGING averaging;
    if(Z_min && Z_max) for(CELL_ITERATOR<TV> iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();
        TV X=iterator.Location()-dt*averaging.Face_To_Cell_Vector(iterator.grid,cell,face_velocities);
        Z(cell)=interpolation.Clamped_To_Array(iterator.grid,Z_ghost,X);
        VECTOR<T2,2> extrema=interpolation.Extrema_Clamped_To_Array_Cell(iterator.grid,*Z_min_ghost,*Z_max_ghost,X);
        (*Z_min)(cell)=extrema.x;(*Z_max)(cell)=extrema.y;}
    else for(CELL_ITERATOR<TV> iterator(grid,domain);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();
        TV cell_velocity(averaging.Face_To_Cell_Vector(iterator.grid,cell,face_velocities));
        cell_velocity*=dt; // TODO: gcc 4.1.2 compiler bug workaround
        Z(cell)=interpolation.Clamped_To_Array(iterator.grid,Z_ghost,iterator.Location()-cell_velocity);
    }
}

template<class TV,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup_Threaded(RANGE<TV_INT>& domain,int axis,const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    T_INTERPOLATION interpolation;T_AVERAGING averaging;
    if(Z_min && Z_max) for(FACE_ITERATOR<TV> iterator(grid,domain,axis);iterator.Valid();iterator.Next()){
        TV_INT face=iterator.Face_Index();int axis=iterator.Axis();TV X=iterator.Location()-dt*averaging.Face_To_Face_Vector(iterator.grid,axis,face,face_velocities);
        Z.Component(axis)(face)=interpolation.Clamped_To_Array_Face_Component(axis,iterator.grid,Z_ghost.Starting_Point_Face(axis,face),X);
        VECTOR<T,2> extrema=interpolation.Extrema_Clamped_To_Array_Face_Component(axis,iterator.grid,Z_min_ghost->Starting_Point_Face(axis,face),Z_max_ghost->Starting_Point_Face(axis,face),X);
        (*Z_min).Component(axis)(face)=extrema.x;(*Z_max).Component(axis)(face)=extrema.y;}
    else for(FACE_ITERATOR<TV> iterator(grid,domain,axis);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
        Z.Component(axis)(face)=interpolation.Clamped_To_Array_Face_Component(axis,iterator.grid,Z_ghost.Starting_Point_Face(axis,face),
            iterator.Location()-dt*averaging.Face_To_Face_Vector(iterator.grid,axis,face,face_velocities));}
}
