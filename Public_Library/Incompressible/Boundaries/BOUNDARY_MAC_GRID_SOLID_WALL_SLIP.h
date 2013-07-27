//#####################################################################
// Copyright 2005-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP
//#####################################################################
#ifndef __BOUNDARY_MAC_GRID_SOLID_WALL_SLIP__
#define __BOUNDARY_MAC_GRID_SOLID_WALL_SLIP__

#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{

template<class TV>
class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP:public BOUNDARY<TV,typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,TV::m> TV_SIDES;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    typedef BOUNDARY<TV,T> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;using BASE::Boundary;using BASE::Find_Ghost_Regions;

    const ARRAY<T,TV_INT>* phi;

    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP(const TV_SIDES& constant_extrapolation=TV_SIDES());
    ~BOUNDARY_MAC_GRID_SOLID_WALL_SLIP();

    void Set_Phi(ARRAY<T,TV_INT>& phi_input)
    {phi=&phi_input;}

//#####################################################################
    void Fill_Ghost_Faces(const GRID<TV>& grid,const T_FACE_ARRAYS_SCALAR& u,T_FACE_ARRAYS_SCALAR& u_ghost,const T time,const int number_of_ghost_cells=3) const PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition_Face(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& u,const T time) const PHYSBAM_OVERRIDE;
    void Reflect_Single_Ghost_Region(const int face_axis,const GRID<TV>& face_grid,T_ARRAYS_BASE& u_ghost_component,const int side,const RANGE<TV_INT>& region) const;
protected:
    void Zero_Single_Boundary_Side(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& u,const int side) const;
//#####################################################################
};
}
#endif
