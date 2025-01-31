//#####################################################################
// Copyright 2010, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Grid_Tools/Grids/NODE_ITERATOR.h>

namespace PhysBAM{

namespace SMOOTH{

//#####################################################################
// Function Smooth
//#####################################################################
// explicit update with CFL=1/2 and all Neumann outer boundaries, no smoothing where phi is negative if defined
template<class T,int dim,class T_ARRAYS_T2>
void Smooth(T_ARRAYS_T2& d,const int steps,const ARRAY<T,VECTOR<int,dim> >* phi){
    typedef typename T_ARRAYS_T2::ELEMENT T2;typedef VECTOR<int,dim> TV_INT;typedef VECTOR<T,dim> TV;

    if(!steps) return;
    int number_of_ghost_cells=-d.domain.min_corner.x;
    RANGE<TV_INT> domain_indices=d.Domain_Indices();
    GRID<TV> grid(domain_indices.Thickened(-number_of_ghost_cells).Maximum_Corner(),RANGE<TV>::Centered_Box());
    if(grid.Domain_Indices().Thickened(number_of_ghost_cells)!=d.Domain_Indices()) PHYSBAM_FATAL_ERROR();

    T_ARRAYS_T2 d_ghost(grid.Domain_Indices(number_of_ghost_cells+1),no_init);
    for(int step=0;step<steps;step++){
        T_ARRAYS_T2::Put(d,d_ghost);
        for(int axis=0;axis<TV::m;axis++)for(int axis_side=0;axis_side<2;axis_side++){
            int side=2*axis+axis_side;TV_INT ghost_offset=(axis_side==0?-1:1)*TV_INT::Axis_Vector(axis);
            for(NODE_ITERATOR<TV> iterator(grid,number_of_ghost_cells,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
                d_ghost(node+ghost_offset)=d_ghost(node);}}
        if(!phi) for(NODE_ITERATOR<TV> iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
            T2 sum=(T)GRID<TV>::number_of_neighbors_per_node*d_ghost(node);
            for(int n=0;n<GRID<TV>::number_of_neighbors_per_node;n++)sum+=d_ghost(GRID<TV>::Node_Neighbor(node,n));
            d(node)=(T)1/(2*GRID<TV>::number_of_neighbors_per_node)*sum;}
        else for(NODE_ITERATOR<TV> iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();if((*phi)(node) > 0){
            T2 sum=(T)GRID<TV>::number_of_neighbors_per_node*d_ghost(node);
            for(int n=0;n<GRID<TV>::number_of_neighbors_per_node;n++)sum+=d_ghost(GRID<TV>::Node_Neighbor(node,n));
            d(node)=(T)1/(2*GRID<TV>::number_of_neighbors_per_node)*sum;}}}
}
template void Smooth<float,2,ARRAY<float,VECTOR<int,2> > >(ARRAY<float,VECTOR<int,2> >&,int,
    ARRAY<float,VECTOR<int,2> > const*);
template void Smooth<float,3,ARRAY<float,VECTOR<int,3> > >(ARRAY<float,VECTOR<int,3> >&,int,
    ARRAY<float,VECTOR<int,3> > const*);
template void Smooth<double,3,ARRAY<VECTOR<double,3>,VECTOR<int,3> > >(
    ARRAY<VECTOR<double,3>,VECTOR<int,3> >&,int,ARRAY<double,VECTOR<int,3> > const*);
template void Smooth<float,3,ARRAY<VECTOR<float,3>,VECTOR<int,3> > >(
    ARRAY<VECTOR<float,3>,VECTOR<int,3> >&,int,ARRAY<float,VECTOR<int,3> > const*);
template void Smooth<double,2,ARRAY<double,VECTOR<int,2> > >(
    ARRAY<double,VECTOR<int,2> >&,int,ARRAY<double,VECTOR<int,2> > const*);
}
}
