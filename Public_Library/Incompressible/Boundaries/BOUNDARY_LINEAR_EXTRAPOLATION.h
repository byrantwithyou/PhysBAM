//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_LINEAR_EXTRAPOLATION
//#####################################################################
#ifndef __BOUNDARY_LINEAR_EXTRAPOLATION__
#define __BOUNDARY_LINEAR_EXTRAPOLATION__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Boundaries/BOUNDARY.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY_LINEAR_EXTRAPOLATION:public BOUNDARY<TV,T2>
{
    typedef typename TV::SCALAR T;typedef BOUNDARY<TV,T2> BASE;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
public:
    using BASE::Find_Ghost_Regions;using BASE::Boundary;

    BOUNDARY_LINEAR_EXTRAPOLATION()
    {}

    virtual ~BOUNDARY_LINEAR_EXTRAPOLATION()
    {}

//#####################################################################
    void Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) const override;
    void Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const override {} // do nothing
//#####################################################################
};
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_LINEAR_EXTRAPOLATION<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    ARRAY<T2,TV_INT>::Put(u,u_ghost); // interior
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=0;axis<TV::m;axis++)for(int axis_side=0;axis_side<2;axis_side++){
        int side=2*axis+axis_side,outward_sign=axis_side?-1:1;
        int boundary=Boundary(side,regions(side));
        TV_INT inward_offset=-outward_sign*TV_INT::Axis_Vector(axis);
        for(NODE_ITERATOR<TV> iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            TV_INT boundary_node=node;boundary_node[axis]=boundary;int ghost_layer=outward_sign*(node[axis]-boundary);
            u_ghost(node)=u_ghost(boundary_node)+ghost_layer*(u_ghost(boundary_node)-u_ghost(boundary_node+inward_offset));}}
}
//#####################################################################
}
#endif
