//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY

//
// The ghost cell routines copy u into u_ghost and add a layer of three ghost cells around the boundary of u.
// The boundary condition routines modify u directly.
//
//#####################################################################
#ifndef __BOUNDARY__
#define __BOUNDARY__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <cassert>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,TV::dimension> TV_SIDES;typedef VECTOR<int,TV::m> TV_INT;
public:
    bool use_fixed_boundary,clamp_below,clamp_above;
    T2 fixed_boundary_value,lower_threshold,upper_threshold;
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;

    BOUNDARY();
    virtual ~BOUNDARY();

    virtual void Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input=TV_SIDES::Constant_Vector(VECTOR<bool,2>::Constant_Vector(true)))
    {constant_extrapolation=constant_extrapolation_input;}

    void Turn_Off_Constant_Extrapolation()
    {Set_Constant_Extrapolation(TV_SIDES());}

    virtual bool Constant_Extrapolation(const int side) const
    {assert((unsigned)side<6);int axis=side/2;return constant_extrapolation(axis)(side&1);}

    virtual void Set_Fixed_Boundary(const bool use_fixed_boundary_input=true,const T2 fixed_boundary_value_input=T2())
    {use_fixed_boundary=use_fixed_boundary_input;fixed_boundary_value=fixed_boundary_value_input;
    if(use_fixed_boundary) clamp_above=clamp_below=false;}

    void Limit_Minimum_Boundary_Value(const bool clamp_below_input=true,const T2 lower_threshold_input=T2())
    {clamp_below=clamp_below_input;lower_threshold=lower_threshold_input;}

    void Limit_Maximum_Boundary_Value(const bool clamp_above_input=true,const T2 upper_threshold_input=T2())
    {clamp_above=clamp_above_input;upper_threshold=upper_threshold_input;}

    void Fill_Ghost_Cells_Cell(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T time,const int number_of_ghost_cells)
    {Fill_Ghost_Cells(grid,u,u_ghost,0,time,number_of_ghost_cells);}

    int Boundary(const int side,const RANGE<TV_INT>& region) const
    {int axis=side/2;return side&1?region.Minimum_Corner()[axis]-1:region.Maximum_Corner()[axis];}

//#####################################################################
    virtual void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time);
    virtual void Apply_Boundary_Condition_Face(const GRID<TV>& grid,ARRAY<T2,FACE_INDEX<TV::m> >& u,const T time);
    virtual void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells);
    virtual void Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells);
    virtual void Apply_Boundary_Condition_Single_Side(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const int side,const T time) const;
    virtual void Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const;
    void Find_Ghost_Regions(const GRID<TV>& grid,ARRAY<RANGE<TV_INT> >& regions,const int number_of_ghost_cells) const;
    void Fill_Single_Ghost_Region(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const int side,const RANGE<TV_INT>& region) const;
    void Fill_Single_Ghost_Region_Threaded(RANGE<TV_INT>& region,const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const int side);
//#####################################################################
};
}
#endif
