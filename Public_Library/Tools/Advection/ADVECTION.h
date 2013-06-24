//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION
//#####################################################################
//
// Solves Z_t + u Z_x + v Z_y + w Z_z = 0 for one Euler step. Conservative schemes solve Z_t +(uZ)_x + (vZ)_y + (wZ)_z = 0.
// Z is updated using u, v, w, & Z_ghost. Input Z as (1,m) by (1,n) by (1,mn), and Z_ghost as (-2,m+3) by (-2,n+3) by (-2,mn+3).
// Input u, v, & w as (1,m) by (1,n) by (1,mn) for nonconservative schemes. Input u, v, & w as (-2,m+3) by (-2,n+3) by (-2,mn+3) for conservative schemes.
//
//#####################################################################
#ifndef __ADVECTION__
#define __ADVECTION__

#include <Tools/Advection/ADVECTION_FORWARD.h>
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T_GRID> struct BOUNDARY_POLICY;
template<class T_GRID> struct GRID_ARRAYS_POLICY;
template<class TV,class T2> class BOUNDARY;
template<class T,class ID> class ARRAY;

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=typename T_GRID::FACE_LOOKUP
class ADVECTION:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
public:

    ADVECTION()
    {}

    virtual ~ADVECTION()
    {}

    void Update_Advection_Equation_Cell(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_ARRAYS_SCALAR& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0)
    {Update_Advection_Equation_Cell_Lookup(grid,Z,Z_ghost,T_FACE_LOOKUP(face_velocities),boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

    void Update_Advection_Equation_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_ghost,
        const T_FACE_ARRAYS_SCALAR& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_ARRAYS_SCALAR* Z_min_ghost=0,const T_FACE_ARRAYS_SCALAR* Z_max_ghost=0,T_FACE_ARRAYS_SCALAR* Z_min=0,T_FACE_ARRAYS_SCALAR* Z_max=0)
    {if(Z_min_ghost && Z_max_ghost){
        T_FACE_LOOKUP Z_min_lookup(*Z_min_ghost),Z_max_lookup(*Z_max_ghost);
        Update_Advection_Equation_Face_Lookup(grid,Z,T_FACE_LOOKUP(Z_ghost),T_FACE_LOOKUP(face_velocities),boundary,dt,time,&Z_min_lookup,&Z_max_lookup,Z_min,Z_max);}
    else Update_Advection_Equation_Face_Lookup(grid,Z,T_FACE_LOOKUP(Z_ghost),T_FACE_LOOKUP(face_velocities),boundary,dt,time,0,0,0,0);}

    void Update_Advection_Equation_Cell(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0)
    {Update_Advection_Equation_Node(grid.Get_Regular_Grid_At_MAC_Positions(),Z,Z_ghost,V,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

//#####################################################################
    virtual void Update_Advection_Equation_Node(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost=0,const T_FACE_LOOKUP* Z_max_ghost=0,T_FACE_ARRAYS_SCALAR* Z_min=0,T_FACE_ARRAYS_SCALAR* Z_max=0)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
