//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_PERIODIC
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Grids_Uniform_Boundaries/BOUNDARY_PERIODIC.h>
using namespace PhysBAM;
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_PERIODIC<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,const int number_of_ghost_cells) const
{
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost); // interior
    TV_INT periods=grid.Domain_Indices().Maximum_Corner();
    VECTOR<RANGE<TV_INT>,2*TV::m> regions;
    Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int axis=0;axis<TV::m;axis++)
        for(int axis_side=0;axis_side<2;axis_side++){
            int side=2*axis+axis_side,sign=1-2*side;
            TV_INT period=sign*periods[axis]*TV_INT::Axis_Vector(axis);
            for(NODE_ITERATOR<TV> iterator(grid,regions(side));iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
                u_ghost(node)=u_ghost(node+period);}}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY_PERIODIC<TV,T2>::
Apply_Boundary_Condition(const GRID<TV>& grid,ARRAYS_ND_BASE<T2,TV_INT>& u,const T time) const
{
    assert(!grid.Is_MAC_Grid());
    for(int axis=0;axis<TV::m;axis++)
        for(NODE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,2*axis);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            TV_INT opposite_node=node;opposite_node[axis]=0;
            u(node)=u(opposite_node);}
}
namespace PhysBAM{
template class BOUNDARY_PERIODIC<VECTOR<float,2>,float>;
template class BOUNDARY_PERIODIC<VECTOR<float,3>,float>;
template class BOUNDARY_PERIODIC<VECTOR<double,2>,double>;
template class BOUNDARY_PERIODIC<VECTOR<double,3>,double>;
}
