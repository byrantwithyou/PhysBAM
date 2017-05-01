//#####################################################################
// Copyright 2004-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW
//#####################################################################
#ifndef __BOUNDARY_SOLID_WALL_SLIP_OUTFLOW__
#define __BOUNDARY_SOLID_WALL_SLIP_OUTFLOW__
#include <Incompressible/Boundaries/BOUNDARY_MAC_GRID_SOLID_WALL_SLIP.h>
namespace PhysBAM{

template<class TV>
class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW:public BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV> BASE;
    using BASE::Constant_Extrapolation;using BASE::lower_threshold;using BASE::upper_threshold;using BASE::clamp_below;using BASE::clamp_above;

    BOUNDARY_SOLID_WALL_SLIP_OUTFLOW(const bool left_constant_extrapolation_input=false,const bool right_constant_extrapolation_input=false,
        const bool bottom_constant_extrapolation_input=false,const bool top_constant_extrapolation_input=false,
        const bool front_constant_extrapolation_input=false,const bool back_constant_extrapolation_input=false)
        :BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV>(left_constant_extrapolation_input,right_constant_extrapolation_input,bottom_constant_extrapolation_input,top_constant_extrapolation_input,
            front_constant_extrapolation_input,back_constant_extrapolation_input)
    {}

    ~BOUNDARY_SOLID_WALL_SLIP_OUTFLOW() = default;

//#####################################################################
    void Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<T,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells=3) const override;
//#####################################################################
};
}
#endif
