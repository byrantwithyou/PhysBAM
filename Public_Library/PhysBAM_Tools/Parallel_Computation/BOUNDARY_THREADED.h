//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_THREADED
//#####################################################################
#ifndef __BOUNDARY_THREADED__
#define __BOUNDARY_THREADED__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/DOMAIN_ITERATOR_THREADED.h>
namespace PhysBAM{

template<class T_GRID> struct BOUNDARY_POLICY;

template<class T_GRID,class T2=typename T_GRID::SCALAR>
class BOUNDARY_THREADED:public BOUNDARY_UNIFORM<T_GRID,T2>
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef VECTOR<T,T_GRID::dimension> TV;typedef VECTOR<int,T_GRID::dimension> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<T2,FACE_INDEX<TV::m> > T_FACE_ARRAYS_T2;
public:
    THREAD_QUEUE& thread_queue;
    BOUNDARY_UNIFORM<T_GRID,T2>& boundary;

    BOUNDARY_THREADED(THREAD_QUEUE& thread_queue_input,BOUNDARY_UNIFORM<T_GRID,T2>& boundary_input)
        :thread_queue(thread_queue_input),boundary(boundary_input)
    {
        assert(&boundary);
    }

    void Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)))
    {boundary.Set_Constant_Extrapolation(constant_extrapolation_input);}

    bool Constant_Extrapolation(const int side) const PHYSBAM_OVERRIDE
    {return boundary.Constant_Extrapolation(side);}

    void Set_Fixed_Boundary(const bool use_fixed_boundary_input=true,const T2 fixed_boundary_value_input=T2())
    {boundary.Set_Fixed_Boundary(use_fixed_boundary_input,fixed_boundary_value_input);}

    void Fill_Ghost_Cells(const T_GRID& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells_input) PHYSBAM_OVERRIDE
    {
        DOMAIN_ITERATOR_THREADED_ALPHA<ARRAYS_ND_BASE<T2,TV_INT>,TV>(u.Domain_Indices(),&thread_queue).template Run<const ARRAYS_ND_BASE<T2,TV_INT>&,ARRAYS_ND_BASE<T2,TV_INT>&>(u_ghost,&ARRAYS_ND_BASE<T2,TV_INT>::Put_With_Range,u,u_ghost);
        ARRAY<RANGE<TV_INT> > regions;boundary.Find_Ghost_Regions(grid,regions,number_of_ghost_cells_input);
        for(int side=0;side<T_GRID::number_of_faces_per_cell;side++){int axis=side/2;
            DOMAIN_ITERATOR_THREADED_ALPHA<BOUNDARY_UNIFORM<T_GRID,T2>,TV>(regions(side),&thread_queue,axis%TV::dimension+1).template Run<const T_GRID&,ARRAYS_ND_BASE<T2,TV_INT>&,int>(boundary,&BOUNDARY_UNIFORM<T_GRID,T2>::Fill_Single_Ghost_Region_Threaded,grid,u_ghost,side);}
    }

    void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells_input) PHYSBAM_OVERRIDE
    {for(int axis=0;axis<T_GRID::dimension;axis++)Fill_Ghost_Cells(grid.Get_Face_Grid(axis),u.Component(axis),u_ghost.Component(axis),0,time,number_of_ghost_cells_input);}

    void Apply_Boundary_Condition(const T_GRID& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) PHYSBAM_OVERRIDE
    {boundary.Apply_Boundary_Condition(grid,u,time);}

    void Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE
    {boundary.Apply_Boundary_Condition_Face(grid,u,time);}

//#####################################################################
};
}
#endif
