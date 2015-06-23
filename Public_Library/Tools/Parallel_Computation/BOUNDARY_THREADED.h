//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_THREADED
//#####################################################################
#ifndef __BOUNDARY_THREADED__
#define __BOUNDARY_THREADED__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
namespace PhysBAM{

template<class TV> struct BOUNDARY_POLICY;

template<class TV,class T2=typename TV::SCALAR>
class BOUNDARY_THREADED:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,TV::m> TV_SIDES;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<T2,FACE_INDEX<TV::m> > T_FACE_ARRAYS_T2;
public:
    THREAD_QUEUE& thread_queue;
    BOUNDARY<TV,T2>& boundary;

    BOUNDARY_THREADED(THREAD_QUEUE& thread_queue_input,BOUNDARY<TV,T2>& boundary_input)
        :thread_queue(thread_queue_input),boundary(boundary_input)
    {
        assert(&boundary);
    }

    void Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)))
    {boundary.Set_Constant_Extrapolation(constant_extrapolation_input);}

    bool Constant_Extrapolation(const int side) const override
    {return boundary.Constant_Extrapolation(side);}

    void Set_Fixed_Boundary(const bool use_fixed_boundary_input=true,const T2 fixed_boundary_value_input=T2())
    {boundary.Set_Fixed_Boundary(use_fixed_boundary_input,fixed_boundary_value_input);}

    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells_input) const override
    {
        DOMAIN_ITERATOR_THREADED_ALPHA<ARRAYS_ND_BASE<T2,TV_INT>,TV>(u.Domain_Indices(),&thread_queue).template Run<const ARRAYS_ND_BASE<T2,TV_INT>&,ARRAYS_ND_BASE<T2,TV_INT>&>(u_ghost,&ARRAYS_ND_BASE<T2,TV_INT>::Put_With_Range,u,u_ghost);
        VECTOR<RANGE<TV_INT>,2*TV::m> regions;boundary.Find_Ghost_Regions(grid,regions,number_of_ghost_cells_input);
        for(int side=0;side<GRID<TV>::number_of_faces_per_cell;side++){int axis=side/2;
            DOMAIN_ITERATOR_THREADED_ALPHA<BOUNDARY<TV,T2>,TV>(regions(side),&thread_queue,axis%TV::dimension+1).template Run<const GRID<TV>&,ARRAYS_ND_BASE<T2,TV_INT>&,int>(boundary,&BOUNDARY<TV,T2>::Fill_Single_Ghost_Region_Threaded,grid,u_ghost,side);}
    }

    void Fill_Ghost_Faces(const GRID<TV>& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells_input) const override
    {for(int axis=0;axis<TV::m;axis++)Fill_Ghost_Cells(grid.Get_Face_Grid(axis),u.Component(axis),u_ghost.Component(axis),0,time,number_of_ghost_cells_input);}

    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override
    {boundary.Apply_Boundary_Condition(grid,u,time);}

    void Apply_Boundary_Condition_Face(const GRID<TV>& grid,T_FACE_ARRAYS_T2& u,const T time) const override
    {boundary.Apply_Boundary_Condition_Face(grid,u,time);}

//#####################################################################
};
}
#endif
