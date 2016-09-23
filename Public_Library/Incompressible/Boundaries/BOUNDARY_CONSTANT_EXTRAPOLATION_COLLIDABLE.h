//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE
//#####################################################################
#ifndef __BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE__
#define __BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE__

#include <Grid_PDE/Boundaries/BOUNDARY.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename GRID<TV>::RIGID_BODY_LIST T_RIGID_BODY_LIST;
public:
    T_RIGID_BODY_LIST& body_list;

    BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE(T_RIGID_BODY_LIST& body_list_input)
        :body_list(body_list_input)
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAY<T2,TV_INT>& u,ARRAY<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAY<T2,TV_INT>& u,const T time) const override {} // do nothing
    void Collision_Aware_Extrapolate(const GRID<TV>& grid,ARRAY<T2,TV_INT>& u_ghost,const TV_INT& source_index,const TV_INT& ghost_index,const TV& direction,const T ray_length);
//#####################################################################
};
}
#endif
