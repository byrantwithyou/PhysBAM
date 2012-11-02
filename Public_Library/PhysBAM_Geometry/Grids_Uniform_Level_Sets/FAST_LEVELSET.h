//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_LEVELSET
//#####################################################################
#ifndef __FAST_LEVELSET__
#define __FAST_LEVELSET__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class TV>
class FAST_LEVELSET:public LEVELSET<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef AVERAGING_UNIFORM<GRID<TV> > T_AVERAGING;
public:
    typedef LEVELSET<TV> BASE;
    using BASE::grid;using BASE::boundary;using BASE::phi;using BASE::curvature_motion;using BASE::sigma;using BASE::max_time_step;
    using BASE::small_number;using BASE::Set_Face_Velocities_Valid_Mask;

    T half_band_width;

    FAST_LEVELSET(GRID<TV>& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells=3);
    ~FAST_LEVELSET();

    void Set_Band_Width(const T number_of_cells=6)
    {half_band_width=number_of_cells*grid.dX.Max()/2;}

//#####################################################################
    T CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const;
    T CFL(const ARRAY<TV,TV_INT>& velocity) const;
public:
    void Fast_Marching_Method(const int local_advection_spatial_order, const T time=0,int process_sign=0);
    void Fast_Marching_Method_Outside_Band(const int local_advection_spatial_order, const T time=0,int process_sign=0);
//#####################################################################
};
}
#endif
