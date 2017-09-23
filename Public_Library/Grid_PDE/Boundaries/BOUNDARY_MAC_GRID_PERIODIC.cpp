//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_MAC_GRID_PERIODIC
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_Tools/Grids/NODE_ITERATOR.h>
#include <Grid_PDE/Boundaries/BOUNDARY_MAC_GRID_PERIODIC.h>
namespace PhysBAM{
//#####################################################################
// Function Apply_Boundary_Condition_Face_Periodic
//#####################################################################
template<class T2,class TV> void
Apply_Boundary_Condition_Face_Periodic(const GRID<TV>& grid,ARRAY<T2,FACE_INDEX<TV::m> >& u,
    const VECTOR<bool,TV::m>& is_periodic)
{
    typedef VECTOR<int,TV::m> TV_INT;
    TV_INT periods=grid.numbers_of_cells;
    for(int axis=0;axis<TV::m;axis++){
        if(!is_periodic(axis)) continue;
        for(FACE_ITERATOR<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,2*axis,axis);it.Valid();it.Next()){
            FACE_INDEX<TV::m> face(it.Full_Index());
            face.index(axis)+=periods(axis);
            u(face)=u(it.Full_Index());}}
}
//#####################################################################
// Function Fill_Ghost_Cells_Periodic
//#####################################################################
template<class T2,class TV> void
Fill_Ghost_Cells_Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u,
    ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells)
{
    typedef VECTOR<int,TV::m> TV_INT;
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost); // interior
    TV_INT periods=grid.numbers_of_cells;
    for(CELL_ITERATOR<TV> it(grid,number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next()){
        TV_INT ind=it.index;
        for(int i=0;i<TV::m;i++)
            if(is_periodic(i))
                ind(i)=wrap(ind(i),0,grid.numbers_of_cells(i));
        u_ghost(it.index)=u(ind);}
}
//#####################################################################
// Function Fill_Ghost_Faces_Periodic
//#####################################################################
template<class T2,class TV> void
Fill_Ghost_Faces_Periodic(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,
    ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells)
{
    typedef VECTOR<int,TV::m> TV_INT;
    assert(grid.Is_MAC_Grid());
    for(int axis=0;axis<TV::m;axis++){
        if(!is_periodic(axis)) continue;
        const ARRAYS_ND_BASE<T2,TV_INT>& u_axis=u.Component(axis);
        ARRAYS_ND_BASE<T2,TV_INT>& u_ghost_axis=u_ghost.Component(axis);
        ARRAYS_ND_BASE<T2,TV_INT>::Put(u_axis,u_ghost_axis); // interior
        GRID<TV> face_grid=grid.Get_Face_Grid(axis);
        TV_INT periods=grid.numbers_of_cells;
        VECTOR<RANGE<TV_INT>,2*TV::m> regions;
        BOUNDARY_MAC_GRID_PERIODIC<TV,typename TV::SCALAR>::Find_Ghost_Regions(face_grid,regions,number_of_ghost_cells);
        for(int face_axis=0;face_axis<TV::m;face_axis++)
            for(int axis_side=0;axis_side<2;axis_side++){
                int side=2*face_axis+axis_side;
                TV_INT period=(axis_side==0?1:-1)*periods[face_axis]*TV_INT::Axis_Vector(face_axis);
                for(NODE_ITERATOR<TV> iterator(face_grid,regions(side));iterator.Valid();iterator.Next()){
                    TV_INT node=iterator.Node_Index();
                    u_ghost_axis(node)=u_ghost_axis(node+period);}}}
}
//#####################################################################
// Function Fill_Ghost_Cells_Periodic
//#####################################################################
template<class T2,class TV> void
Fill_Ghost_Cells_Periodic_Accum(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u,
    ARRAYS_ND_BASE<T2,VECTOR<int,TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells)
{
    typedef VECTOR<int,TV::m> TV_INT;
    ARRAYS_ND_BASE<T2,TV_INT>::Put(u,u_ghost); // interior
    for(CELL_ITERATOR<TV> it(grid,number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next()){
        TV_INT ind(wrap(it.index,TV_INT(),grid.numbers_of_cells));
        u_ghost(ind)+=u(it.index);}
    Fill_Ghost_Cells_Periodic(grid,u_ghost,u_ghost,is_periodic,number_of_ghost_cells);
}
//#####################################################################
// Function Fill_Ghost_Faces_Periodic
//#####################################################################
template<class T2,class TV> void
Fill_Ghost_Faces_Periodic_Accum(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,
    ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const VECTOR<bool,TV::m>& is_periodic,
    const int number_of_ghost_cells)
{
    typedef VECTOR<int,TV::m> TV_INT;
    assert(grid.Is_MAC_Grid());
    ARRAY<T2,FACE_INDEX<TV::m> >::Put(u,u_ghost); // interior
    TV_INT periods=grid.numbers_of_cells;
    for(FACE_ITERATOR<TV> it(grid,number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next()){
        FACE_INDEX<TV::m> fit(it.Full_Index());
        for(int i=0;i<TV::m;i++)
            if(is_periodic(i))
                fit.index(i)=wrap(fit.index(i),0,grid.numbers_of_cells(i));
        if(fit.index!=it.face.index)
            u_ghost(fit)+=u(it.Full_Index());}
    for(int axis=0;axis<TV::m;axis++){
        if(!is_periodic(axis)) continue;
        for(FACE_ITERATOR<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,2*axis,axis);it.Valid();it.Next()){
            FACE_INDEX<TV::m> face(it.Full_Index());
            face.index(axis)+=periods(axis);
            u_ghost(it.Full_Index())+=u(face);
            u_ghost(face)=u_ghost(it.Full_Index());}}
    for(FACE_ITERATOR<TV> it(grid,number_of_ghost_cells,GRID<TV>::GHOST_REGION);it.Valid();it.Next()){
        FACE_INDEX<TV::m> fit(it.Full_Index());
        for(int i=0;i<TV::m;i++)
            if(is_periodic(i))
                fit.index(i)=wrap(fit.index(i),0,grid.numbers_of_cells(i));
        u_ghost(it.Full_Index())=u_ghost(fit);}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class TV,class T2> void BOUNDARY_MAC_GRID_PERIODIC<TV,T2>::
Apply_Boundary_Condition_Face(const GRID<TV>& grid,ARRAY<T2,FACE_INDEX<TV::m> >& u,
    const T time) const
{
    Apply_Boundary_Condition_Face_Periodic(grid,u,is_periodic);
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class TV,class T2> void BOUNDARY_MAC_GRID_PERIODIC<TV,T2>::
Fill_Ghost_Cells(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,
    ARRAYS_ND_BASE<T2,TV_INT>& u_ghost,const T dt,const T time,
    const int number_of_ghost_cells) const
{
    Fill_Ghost_Cells_Periodic(grid,u,u_ghost,is_periodic,number_of_ghost_cells);
}
//#####################################################################
// Function Fill_Ghost_Faces
//#####################################################################
template<class TV,class T2> void BOUNDARY_MAC_GRID_PERIODIC<TV,T2>::
Fill_Ghost_Faces(const GRID<TV>& grid,const ARRAY<T2,FACE_INDEX<TV::m> >& u,
    ARRAY<T2,FACE_INDEX<TV::m> >& u_ghost,const T time,
    const int number_of_ghost_cells) const
{
    Fill_Ghost_Faces_Periodic(grid,u,u_ghost,is_periodic,number_of_ghost_cells);
}
//#####################################################################
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,2>,VECTOR<float,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,2>,MATRIX<double,2> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,2>,VECTOR<double,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,2>,double>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,2>,int>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,3>,MATRIX<double,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,3>,double>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<double,3>,int>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,2>,MATRIX<float,2> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,2>,float>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,2>,int>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,3>,MATRIX<float,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3> >;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,3>,float>;
template class BOUNDARY_MAC_GRID_PERIODIC<VECTOR<float,3>,int>;
template void Apply_Boundary_Condition_Face_Periodic<bool,VECTOR<double,2> >(
    GRID<VECTOR<double,2> > const&,ARRAY<bool,FACE_INDEX<VECTOR<double,2>::m> >&,
    VECTOR<bool,VECTOR<double,2>::m> const&);
template void Apply_Boundary_Condition_Face_Periodic<bool,VECTOR<double,3> >(
    GRID<VECTOR<double,3> > const&,ARRAY<bool,FACE_INDEX<VECTOR<double,3>::m> >&,
    VECTOR<bool,VECTOR<double,3>::m> const&);
template void Apply_Boundary_Condition_Face_Periodic<bool,VECTOR<float,2> >(
    GRID<VECTOR<float,2> > const&,ARRAY<bool,FACE_INDEX<VECTOR<float,2>::m> >&,
    VECTOR<bool,VECTOR<float,2>::m> const&);
template void Apply_Boundary_Condition_Face_Periodic<bool,VECTOR<float,3> >(
    GRID<VECTOR<float,3> > const&,ARRAY<bool,FACE_INDEX<VECTOR<float,3>::m> >&,
    VECTOR<bool,VECTOR<float,3>::m> const&);
template void Fill_Ghost_Faces_Periodic<bool,VECTOR<double,2> >(
    GRID<VECTOR<double,2> > const&,ARRAY<bool,FACE_INDEX<VECTOR<double,2>::m> > const&,
    ARRAY<bool,FACE_INDEX<VECTOR<double,2>::m> >&,VECTOR<bool,VECTOR<double,2>::m> const&,int);
template void Fill_Ghost_Faces_Periodic<bool,VECTOR<double,3> >(
    GRID<VECTOR<double,3> > const&,ARRAY<bool,FACE_INDEX<VECTOR<double,3>::m> > const&,
    ARRAY<bool,FACE_INDEX<VECTOR<double,3>::m> >&,VECTOR<bool,VECTOR<double,3>::m> const&,int);
template void Fill_Ghost_Faces_Periodic<bool,VECTOR<float,2> >(
    GRID<VECTOR<float,2> > const&,ARRAY<bool,FACE_INDEX<VECTOR<float,2>::m> > const&,
    ARRAY<bool,FACE_INDEX<VECTOR<float,2>::m> >&,VECTOR<bool,VECTOR<float,2>::m> const&,int);
template void Fill_Ghost_Faces_Periodic<bool,VECTOR<float,3> >(
    GRID<VECTOR<float,3> > const&,ARRAY<bool,FACE_INDEX<VECTOR<float,3>::m> > const&,
    ARRAY<bool,FACE_INDEX<VECTOR<float,3>::m> >&,VECTOR<bool,VECTOR<float,3>::m> const&,int);
template void Fill_Ghost_Cells_Periodic<int,VECTOR<double,2> >(GRID<VECTOR<double,2> > const&,
    ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<double,2>::m> > const&,
    ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<double,2>::m> >&,
    VECTOR<bool,VECTOR<double,2>::m> const&,int);
template void Fill_Ghost_Cells_Periodic<int,VECTOR<float,2> >(GRID<VECTOR<float,2> > const&,
    ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<float,2>::m> > const&,
    ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<float,2>::m> >&,
    VECTOR<bool,VECTOR<float,2>::m> const&,int);
template void Fill_Ghost_Faces_Periodic_Accum<double,VECTOR<double,2> >(
    GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<VECTOR<double,2>::m> > const&,
    ARRAY<double,FACE_INDEX<VECTOR<double,2>::m> >&,
    VECTOR<bool,VECTOR<double,2>::m> const&,int);
template void Fill_Ghost_Faces_Periodic_Accum<double,VECTOR<double,3> >(
    GRID<VECTOR<double,3> > const&,ARRAY<double,FACE_INDEX<VECTOR<double,3>::m> > const&,
    ARRAY<double,FACE_INDEX<VECTOR<double,3>::m> >&,
    VECTOR<bool,VECTOR<double,3>::m> const&,int);
template void Fill_Ghost_Faces_Periodic_Accum<float,VECTOR<float,2> >(
    GRID<VECTOR<float,2> > const&,ARRAY<float,FACE_INDEX<VECTOR<float,2>::m> > const&,
    ARRAY<float,FACE_INDEX<VECTOR<float,2>::m> >&,
    VECTOR<bool,VECTOR<float,2>::m> const&,int);
template void Fill_Ghost_Faces_Periodic_Accum<float,VECTOR<float,3> >(
    GRID<VECTOR<float,3> > const&,ARRAY<float,FACE_INDEX<VECTOR<float,3>::m> > const&,
    ARRAY<float,FACE_INDEX<VECTOR<float,3>::m> >&,
    VECTOR<bool,VECTOR<float,3>::m> const&,int);
template void Fill_Ghost_Faces_Periodic<float,VECTOR<float,2> >(
    GRID<VECTOR<float,2> > const&,ARRAY<float,FACE_INDEX<VECTOR<float,2>::m> > const&,
    ARRAY<float,FACE_INDEX<VECTOR<float,2>::m> >&,VECTOR<bool,VECTOR<float,2>::m> const&,
    int);
template void Fill_Ghost_Faces_Periodic<double,VECTOR<double,2> >(
    GRID<VECTOR<double,2> > const&,ARRAY<double,FACE_INDEX<VECTOR<double,2>::m> > const&,
    ARRAY<double,FACE_INDEX<VECTOR<double,2>::m> >&,VECTOR<bool,VECTOR<double,2>::m> const&,
    int);
template void Fill_Ghost_Cells_Periodic<int,VECTOR<double,3> >(
    GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<double,3>::m> > const&,
    ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<double,3>::m> >&,VECTOR<bool,VECTOR<double,3>::m> const&,
    int);
template void Fill_Ghost_Cells_Periodic<int,VECTOR<float,3> >(
    GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<float,3>::m> > const&,
    ARRAYS_ND_BASE<int,VECTOR<int,VECTOR<float,3>::m> >&,VECTOR<bool,VECTOR<float,3>::m> const&,
    int);
}
