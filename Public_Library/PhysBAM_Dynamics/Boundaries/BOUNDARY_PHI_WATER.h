//#####################################################################
// Copyright 2004-2006, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_PHI_WATER
//#####################################################################
#ifndef __BOUNDARY_PHI_WATER__
#define __BOUNDARY_PHI_WATER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_OPEN_CALLBACKS.h>
namespace PhysBAM{

template<class T_GRID>
class BOUNDARY_PHI_WATER:public BOUNDARY_UNIFORM<T_GRID,typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef BOUNDARY_UNIFORM<T_GRID,T> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Find_Ghost_Regions;using BASE::Boundary;

    bool use_extrapolation_mode;
    T tolerance;
    int sign; // used for multiphase. Set to -1 for multiphase air region

    bool use_open_boundary_mode;
    ARRAY<bool> open_boundary;
    const BOUNDARY_OPEN_CALLBACKS<T_GRID> *callbacks;
private:
    const T_FACE_ARRAYS_SCALAR* V;
public:

    BOUNDARY_PHI_WATER(const TV_SIDES& constant_extrapolation=TV_SIDES());
    ~BOUNDARY_PHI_WATER();

    void Use_Extrapolation_Mode(const bool use=true)
    {use_extrapolation_mode=use;}

    void Set_Tolerance(const T tolerance_input=(T)9.8/24)  // dt*gravity where dt=1/24 is based on the length of a frame
    {tolerance=tolerance_input;}

    void Set_Velocity_Pointer(const T_FACE_ARRAYS_SCALAR& V_input)
    {V=&V_input;}

    void Set_Open_Boundary(const bool left_open_boundary_input=false,const bool right_open_boundary_input=false,const bool bottom_open_boundary_input=false,const bool top_open_boundary_input=false,const bool front_open_boundary_input=false,const bool back_open_boundary_input=false)
    {
        open_boundary(0)=left_open_boundary_input;
        open_boundary(1)=right_open_boundary_input;
        open_boundary(2)=bottom_open_boundary_input;
        open_boundary(3)=top_open_boundary_input;
        open_boundary(4)=front_open_boundary_input;
        open_boundary(5)=back_open_boundary_input;
    }

    void Use_Open_Boundary_Mode(const bool use=true)
    {use_open_boundary_mode=use;}

    void Set_Boundary_Open_Callbacks(const BOUNDARY_OPEN_CALLBACKS<T_GRID>& callbacks_input)
    {callbacks=&callbacks_input;}

//#####################################################################
    void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_BASE& u,T_ARRAYS_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells); // uniform grids
    void Fill_Single_Ghost_Region_Threaded(RANGE<TV_INT>& region,const T_GRID& grid,T_ARRAYS_BASE& u_ghost,const int side);
//#####################################################################
};
}
#endif
