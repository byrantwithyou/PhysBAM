//#####################################################################
// Copyright 2002-2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_OPEN_WATER
//#####################################################################
#ifndef __BOUNDARY_OPEN_WATER__
#define __BOUNDARY_OPEN_WATER__

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_PDE/Boundaries/BOUNDARY.h>
#include <Incompressible/Boundaries/BOUNDARY_FORWARD.h>
namespace PhysBAM{

template<class TV>
class BOUNDARY_OPEN_WATER:public BOUNDARY<TV,typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    typedef BOUNDARY<TV,T> BASE;
    using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;using BASE::Find_Ghost_Regions;

    ARRAY<bool> open_boundary;
    T attenuate_inflow;
    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<TV> boundary_mac_grid_solid_wall_slip;

    BOUNDARY_OPEN_WATER(const T attenuate_inflow_input=T(1),const bool left_open_boundary_input=false,const bool right_open_boundary_input=false,const bool bottom_open_boundary_input=false,const bool top_open_boundary_input=false,
        const bool front_open_boundary_input=false,const bool back_open_boundary_input=false)
        :open_boundary(6,true),boundary_mac_grid_solid_wall_slip()
    {
        attenuate_inflow=attenuate_inflow_input;
        open_boundary(0)=left_open_boundary_input;
        open_boundary(1)=right_open_boundary_input;
        open_boundary(2)=bottom_open_boundary_input;
        open_boundary(3)=top_open_boundary_input;
        open_boundary(4)=front_open_boundary_input;
        open_boundary(5)=back_open_boundary_input;
    }

    ~BOUNDARY_OPEN_WATER(){}

public:

//#####################################################################
    void Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,ARRAY<T,FACE_INDEX<TV::m> >& u_ghost,const T time,const int number_of_ghost_cells=3) const override;
//#####################################################################
};
}
#endif
