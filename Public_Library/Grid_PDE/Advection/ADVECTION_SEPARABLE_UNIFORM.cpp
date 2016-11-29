//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEPARABLE_UNIFORM
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_PDE/Advection/ADVECTION_SEPARABLE_UNIFORM.h>
using namespace PhysBAM;
namespace{
// The following construct is necessary over straight function overloading to
// keep MSVC9 from getting confused :(
template<int d> struct UPDATE_ADVECTION_EQUATION_HELPER;
//#####################################################################
// Struct UPDATE_ADVECTION_EQUATION_HELPER<1>
//#####################################################################
template<> struct UPDATE_ADVECTION_EQUATION_HELPER<1>
{
    template<class T_ADVECTION_SEPARABLE_UNIFORM,class T,class T2>
    static void Apply(T_ADVECTION_SEPARABLE_UNIFORM& advection,const GRID<VECTOR<T,1> >& grid,ARRAYS_ND_BASE<T2,VECTOR<int,1> >& Z,
        const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& Z_ghost,const ARRAYS_ND_BASE<VECTOR<T,1>,VECTOR<int,1> >& V,const T dt,const T time)
    {
        int i;int m=grid.counts.x;T dx=grid.dX.x;ARRAY<T2,VECTOR<int,1> > rhs(0,m);

        ARRAY<T2,VECTOR<int,1> > Z_1d_x(-3,m+3);ARRAY<T,VECTOR<int,1> > u_1d(V.domain.min_corner.x,V.domain.max_corner.x);
        for(i=-3;i<m+3;i++) Z_1d_x(i)=Z_ghost(i);
        for(i=V.domain.min_corner.x;i<V.domain.max_corner.x;i++) u_1d(i)=V(i).x;
        advection.Advection_Solver(m,dx,Z_1d_x,u_1d,rhs);

        for(i=0;i<m;i++) Z(i)-=dt*rhs(i);
    }
};
//#####################################################################
// Struct UPDATE_ADVECTION_EQUATION_HELPER<2>
//#####################################################################
template<> struct UPDATE_ADVECTION_EQUATION_HELPER<2>
{
    template<class T_ADVECTION_SEPARABLE_UNIFORM,class T,class T2>
    static void Apply(T_ADVECTION_SEPARABLE_UNIFORM& advection,const GRID<VECTOR<T,2> >& grid,ARRAYS_ND_BASE<T2,VECTOR<int,2> >& Z,
        const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& Z_ghost,const ARRAYS_ND_BASE<VECTOR<T,2>,VECTOR<int,2> >& V,const T dt,const T time)
    {
        int i,j;int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;ARRAY<T2,VECTOR<int,2> > rhs(0,m,0,n);

        {ARRAY<T2,VECTOR<int,1> > Z_1d_x(-3,m+3),u_Zx_1d(0,m);ARRAY<T,VECTOR<int,1> > u_1d(V.domain.min_corner.x,V.domain.max_corner.x);
        for(j=0;j<n;j++){
            for(i=-3;i<m+3;i++) Z_1d_x(i)=Z_ghost(i,j);
            for(i=V.domain.min_corner.x;i<V.domain.max_corner.x;i++) u_1d(i)=V(i,j).x;
            advection.Advection_Solver(m,dx,Z_1d_x,u_1d,u_Zx_1d);
            for(i=0;i<m;i++) rhs(i,j)=u_Zx_1d(i);}}

        {ARRAY<T2,VECTOR<int,1> > Z_1d_y(-3,n+3),v_Zy_1d(0,n);ARRAY<T,VECTOR<int,1> > v_1d(V.domain.min_corner.y,V.domain.max_corner.y);
        for(i=0;i<m;i++){
            for(j=-3;j<n+3;j++) Z_1d_y(j)=Z_ghost(i,j);
            for(j=V.domain.min_corner.y;j<V.domain.max_corner.y;j++) v_1d(j)=V(i,j).y;
            advection.Advection_Solver(n,dy,Z_1d_y,v_1d,v_Zy_1d);
            for(j=0;j<n;j++) rhs(i,j)+=v_Zy_1d(j);}}

        for(i=0;i<m;i++) for(j=0;j<n;j++) Z(i,j)-=dt*rhs(i,j);
    }
};
//#####################################################################
// Struct UPDATE_ADVECTION_EQUATION_HELPER<3>
//#####################################################################
template<> struct UPDATE_ADVECTION_EQUATION_HELPER<3>
{
    template<class T_ADVECTION_SEPARABLE_UNIFORM,class T,class T2>
    static void Apply(T_ADVECTION_SEPARABLE_UNIFORM& advection,const GRID<VECTOR<T,3> >& grid,ARRAYS_ND_BASE<T2,VECTOR<int,3> >& Z,
        const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& Z_ghost,const ARRAYS_ND_BASE<VECTOR<T,3>,VECTOR<int,3> >& V,const T dt,const T time)
    {
        int i,j,ij;int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;ARRAY<T2,VECTOR<int,3> > rhs(0,m,0,n,0,mn);

        {ARRAY<T2,VECTOR<int,1> > Z_1d_x(-3,m+3),u_Zx_1d(0,m);ARRAY<T,VECTOR<int,1> > u_1d(V.domain.min_corner.x,V.domain.max_corner.x);
        for(j=0;j<n;j++) for(ij=0;ij<mn;ij++){
            for(i=-3;i<m+3;i++) Z_1d_x(i)=Z_ghost(i,j,ij);
            for(i=V.domain.min_corner.x;i<V.domain.max_corner.x;i++) u_1d(i)=V(i,j,ij).x;
            advection.Advection_Solver(m,dx,Z_1d_x,u_1d,u_Zx_1d);
            for(i=0;i<m;i++) rhs(i,j,ij)=u_Zx_1d(i);}}

        {ARRAY<T2,VECTOR<int,1> > Z_1d_y(-3,n+3),v_Zy_1d(0,n);ARRAY<T,VECTOR<int,1> > v_1d(V.domain.min_corner.y,V.domain.max_corner.y);
        for(i=0;i<m;i++) for(ij=0;ij<mn;ij++){
            for(j=-3;j<n+3;j++) Z_1d_y(j)=Z_ghost(i,j,ij);
            for(j=V.domain.min_corner.y;j<V.domain.max_corner.y;j++) v_1d(j)=V(i,j,ij).y;
            advection.Advection_Solver(n,dy,Z_1d_y,v_1d,v_Zy_1d);
            for(j=0;j<n;j++) rhs(i,j,ij)+=v_Zy_1d(j);}}

        {ARRAY<T2,VECTOR<int,1> > Z_1d_z(-3,mn+3),w_Zz_1d(0,mn);ARRAY<T,VECTOR<int,1> > w_1d(V.domain.min_corner.z,V.domain.max_corner.z);
        for(i=0;i<m;i++) for(j=0;j<n;j++){
            for(ij=-3;ij<mn+3;ij++) Z_1d_z(ij)=Z_ghost(i,j,ij);
            for(ij=V.domain.min_corner.z;ij<V.domain.max_corner.z;ij++) w_1d(ij)=V(i,j,ij).z;
            advection.Advection_Solver(mn,dz,Z_1d_z,w_1d,w_Zz_1d);
            for(ij=0;ij<mn;ij++) rhs(i,j,ij)+=w_Zz_1d(ij);}}

        for(RANGE_ITERATOR<3> it(grid.Domain_Indices());it.Valid();it.Next()) Z(it.index)-=dt*rhs(it.index);
    }
};
}
//#####################################################################
// Function Update_Advection_Equation_Node
//#####################################################################
template<class TV,class T2,class T_AVERAGING> void ADVECTION_SEPARABLE_UNIFORM<TV,T2,T_AVERAGING>::
Update_Advection_Equation_Node(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    assert(!Z_min && !Z_max);
    UPDATE_ADVECTION_EQUATION_HELPER<TV::m>::Apply(*this,grid,Z,Z_ghost,V,dt,time);
}
//#####################################################################
// Function Update_Advection_Equation_Cell_Lookup
//#####################################################################
template<class TV,class T2,class T_AVERAGING> void ADVECTION_SEPARABLE_UNIFORM<TV,T2,T_AVERAGING>::
Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,const T_FACE_LOOKUP& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
    const ARRAY<T2,TV_INT>* Z_min_ghost,const ARRAY<T2,TV_INT>* Z_max_ghost,ARRAY<T2,TV_INT>* Z_min,ARRAY<T2,TV_INT>* Z_max)
{
    assert(!Z_min && !Z_max);

    ARRAY<TV,TV_INT> V_cell(grid.Domain_Indices(3));T_AVERAGING averaging;

    for(CELL_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) V_cell(iterator.Cell_Index())=averaging.Face_To_Cell_Vector(grid,iterator.Cell_Index(),V);

    UPDATE_ADVECTION_EQUATION_HELPER<TV::m>::Apply(*this,grid,Z,Z_ghost,V_cell,dt,time);
}
//#####################################################################
// Function Update_Advection_Equation_Cell
//#####################################################################
template<class TV,class T2,class T_AVERAGING> void ADVECTION_SEPARABLE_UNIFORM<TV,T2,T_AVERAGING>::
Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const T_FACE_LOOKUP& Z_ghost,const T_FACE_LOOKUP& V,BOUNDARY<TV,T>& boundary,
    const T dt,const T time,const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,ARRAY<T,FACE_INDEX<TV::m> >* Z_min,ARRAY<T,FACE_INDEX<TV::m> >* Z_max)
{
    assert(!Z_min && !Z_max);

    ARRAY<TV,FACE_INDEX<TV::m> > V_face(grid,0);T_AVERAGING averaging;

    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()) V_face(iterator.Full_Index())=averaging.Face_To_Face_Vector(grid,iterator.Full_Index(),V);

    for(int i=0;i<TV::m;i++){
        GRID<TV> node_grid(grid.Get_Face_Grid(i));
        UPDATE_ADVECTION_EQUATION_HELPER<TV::m>::Apply(*this,node_grid,Z.Component(i),Z_ghost.V_face.Component(i),V_face.Component(i),dt,time);}
}
//#####################################################################
namespace PhysBAM{
template class ADVECTION_SEPARABLE_UNIFORM<VECTOR<float,1>,float>;
template class ADVECTION_SEPARABLE_UNIFORM<VECTOR<float,2>,float>;
template class ADVECTION_SEPARABLE_UNIFORM<VECTOR<float,3>,float>;
template class ADVECTION_SEPARABLE_UNIFORM<VECTOR<double,1>,double>;
template class ADVECTION_SEPARABLE_UNIFORM<VECTOR<double,2>,double>;
template class ADVECTION_SEPARABLE_UNIFORM<VECTOR<double,3>,double>;
}
