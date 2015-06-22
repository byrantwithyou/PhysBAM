//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Log/DEBUG_SUBSTEPS.h>
#include <Tools/Math_Tools/INTERVAL.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Geometry/Level_Sets/REINITIALIZATION.h>
#include <climits>
using namespace PhysBAM;
static void set_if_abs_less(float a,float b,float& x,float y,float z){if(abs(a)<abs(b)) x=y;else x=z;}
static void set_if_abs_less(double a,double b,double& x,double y,double z){if(abs(a)<abs(b)) x=y;else x=z;}
template<class T>
static void set_if_abs_less(const SYMMETRIC_MATRIX<T,3>& a,const SYMMETRIC_MATRIX<T,3>& b,SYMMETRIC_MATRIX<T,3>& x,const SYMMETRIC_MATRIX<T,3>& y,const SYMMETRIC_MATRIX<T,3>& z)
{
    set_if_abs_less(a.x00,b.x00,x.x00,y.x00,z.x00);
    set_if_abs_less(a.x10,b.x10,x.x10,y.x10,z.x10);
    set_if_abs_less(a.x20,b.x20,x.x20,y.x20,z.x20);
    set_if_abs_less(a.x11,b.x11,x.x11,y.x11,z.x11);
    set_if_abs_less(a.x21,b.x21,x.x21,y.x21,z.x21);
    set_if_abs_less(a.x22,b.x22,x.x22,y.x22,z.x22);
}
template<class T>
static void set_if_abs_less(const SYMMETRIC_MATRIX<T,2>& a,const SYMMETRIC_MATRIX<T,2>& b,SYMMETRIC_MATRIX<T,2>& x,const SYMMETRIC_MATRIX<T,2>& y,const SYMMETRIC_MATRIX<T,2>& z)
{
    set_if_abs_less(a.x00,b.x00,x.x00,y.x00,z.x00);
    set_if_abs_less(a.x10,b.x10,x.x10,y.x10,z.x10);
    set_if_abs_less(a.x11,b.x11,x.x11,y.x11,z.x11);
}
template<class T>
static void set_if_abs_less(const SYMMETRIC_MATRIX<T,1>& a,const SYMMETRIC_MATRIX<T,1>& b,SYMMETRIC_MATRIX<T,1>& x,const SYMMETRIC_MATRIX<T,1>& y,const SYMMETRIC_MATRIX<T,1>& z)
{
    set_if_abs_less(a.x00,b.x00,x.x00,y.x00,z.x00);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
EXTRAPOLATION_HIGHER_ORDER(const GRID<TV>& grid,const LEVELSET<TV>& phi,int iterations,int order,int fill_width)
    :grid(grid),phi(phi),iterations(iterations),order(order),fill_width(fill_width)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
~EXTRAPOLATION_HIGHER_ORDER()
{
}
//#####################################################################
// Function Add_Neighbors
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Add_Neighbors(ARRAY<TV_INT>& next,const ARRAY<TV_INT>& neighbors,const TV_INT& index,int unregistered,int registered)
{
    for(int i=0;i<neighbors.m;i++){
        TV_INT neighbor_index=index+neighbors(i);
        Periodic_Index(neighbor_index);
        int& n=node_to_index(neighbor_index);
        if(n!=unregistered) continue;
        next.Append(neighbor_index);
        n=registered;}
}
//#####################################################################
// Function Register_Nodes
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Register_Nodes(std::function<bool(const TV_INT& index)> inside_mask,ARRAY<VECTOR<STENCIL,TV::m> >& stencil)
{
    node_to_index.Resize(grid.Domain_Indices(periodic?0:fill_width+1)); // Need an extra ring for the sentinals
    index_to_node.Append(TV_INT()+INT_MAX); // First index is the "outside" index.
    normal.Append(TV()+FLT_MAX);

    ARRAY<TV_INT> next_inside,next_outside,current;
    ARRAY<TV_INT> neighbors_outside,neighbors_inside;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT()-1,TV_INT()+2));it.Valid();it.Next())
        neighbors_outside.Append(it.index);
    for(int d=0;d<TV::m;d++)
        for(int s=-1;s<=1;s+=2){
            TV_INT index;
            index(d)=s;
            neighbors_inside.Append(index);}

    if(!periodic)
        for(NODE_ITERATOR<TV> it(grid,fill_width,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
            node_to_index(it.index)=-1;

    for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        bool b=inside_mask(it.index);
        node_to_index(it.index)=-1-b;}

    if(!periodic)
        for(NODE_ITERATOR<TV> it(grid,1,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
            Add_Neighbors(next_inside,neighbors_inside,it.index,-2,-4);

    for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
        bool b=inside_mask(it.index);
        if(b) Add_Neighbors(next_outside,neighbors_outside,it.index,-1,-3);
        else Add_Neighbors(next_inside,neighbors_inside,it.index,-2,-4);}

    for(int o=0;o<fill_width;o++){
        current.Exchange(next_outside);
        for(int j=0;j<current.m;j++){
            Register_Index(current(j));
            Add_Neighbors(next_outside,neighbors_outside,current(j),-1,-3);}
        current.Remove_All();}
    next_outside.Remove_All();

    for(int i=1;i<index_to_node.m;i++){
        TV N=phi.Normal(grid.X(index_to_node(i)));
        normal.Append(N);
        for(int d=0;d<TV::m;d++){
            TV_INT index=index_to_node(i);
            int s=N(d)<0?1:-1;
            for(int j=0;j<3;j++){
                Register_Index(index,-3);
                index(d)+=s;}}}

    for(int o=0;o<order;o++){
        solve_indices(o)=INTERVAL<int>(1,index_to_node.m);
        fill_indices(o).min_corner=index_to_node.m;
        current.Exchange(next_inside);
        for(int j=0;j<current.m;j++){
            Register_Index(current(j));
            Add_Neighbors(next_inside,neighbors_inside,current(j),-2,-4);}
        current.Remove_All();}

    for(int o=order-1;o>=0;o--){
        current.Exchange(next_inside);
        for(int j=0;j<current.m;j++){
            Register_Index(current(j));
            Add_Neighbors(next_inside,neighbors_inside,current(j),-2,-4);}
        current.Remove_All();
        fill_indices(o).max_corner=index_to_node.m;}
    next_inside.Remove_All();

    for(int i=normal.m;i<fill_indices(1).max_corner;i++)
        normal.Append(phi.Normal(grid.X(index_to_node(i))));

    VECTOR<STENCIL,TV::m> st;
    stencil.Append(st);
    for(int i=solve_indices(order-1).min_corner;i<solve_indices(order-1).max_corner;i++){
        TV N=normal(i);
        for(int d=0;d<TV::m;d++){
            TV_INT index=index_to_node(i);
            int s=N(d)<0?1:-1;
            for(int j=0;j<3;j++){
                st(d).nodes(j+1)=Lookup_Index(index);
                st(d).scale=s*N(d)*grid.one_over_dX(d);
                index(d)+=s;}
            index(d)-=4*s;
            st(d).nodes(0)=Lookup_Index(index);}
        stencil.Append(st);}
}
//#####################################################################
// Function Fill_un
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Fill_un(const TV& one_over_dx,const ARRAY<T2>& x,ARRAY<T2>& xn,int o,int mo)
{
    xn.Resize(fill_indices(o).max_corner);
    for(int i=fill_indices(o).min_corner;i<fill_indices(o).max_corner;i++){
        const TV_INT& index=index_to_node(i);
        T2 v=T2();
        for(int d=0;d<TV::m;d++){
            TV_INT a=index,b=index;
            a(d)--;
            b(d)++;
            Periodic_Index(a);
            Periodic_Index(b);
            v+=(x(node_to_index(b))-x(node_to_index(a)))*(T).5*one_over_dx(d)*normal(i)(d);}
        xn(i)=v;}
}
//#####################################################################
// Function Constant_Extrapolate_FE
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_FE(const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,const ARRAY<T2>& u,ARRAY<T2>& y,const ARRAY<T2>* z,int o,T dt,T alpha)
{
    for(int i=solve_indices(o).min_corner;i<solve_indices(o).max_corner;i++){
        T2 dot=T2(),zi=z?(*z)(i):T2(),a=T2();
        for(int d=0;d<TV::m;d++){
            // Second order ENO.
            const STENCIL& s=stencil(i)(d);
            VECTOR<T2,4> f(u(s.nodes(0)),u(s.nodes(1)),u(s.nodes(2)),u(s.nodes(3)));
            VECTOR<T2,3> df(f(1)-f(0),f(2)-f(1),f(3)-f(2));
            VECTOR<T2,2> ddf(df(1)-df(0),df(2)-df(1));
            set_if_abs_less(ddf(0),ddf(1),a,(T).5*(f(2)-f(0)),df(1)-(T).5*ddf(1));
            dot+=a*s.scale;}
        y(i)+=alpha*(u(i)-y(i)-dt*(dot-zi));}
}
//#####################################################################
// Function Constant_Extrapolate_RK2
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_RK2(const ARRAY<VECTOR<STENCIL,TV::m> >& stencil,ARRAY<T2>& u,const ARRAY<T2>* z,ARRAY<T2>& tmp,int o,T dt)
{
    tmp.Resize(u.m);
    for(int i=fill_indices(o).min_corner;i<fill_indices(o).max_corner;i++) tmp(i)=u(i);
    Extrapolate_FE(stencil,u,tmp,z,o,dt,1);
    Extrapolate_FE(stencil,tmp,u,z,o,dt,(T).5);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Node(std::function<bool(const TV_INT& index)> inside_mask,ARRAYS_ND_BASE<T2,TV_INT>& u)
{
    PHYSBAM_ASSERT(order>=1 && order<=3);
    PHYSBAM_ASSERT(!grid.Is_MAC_Grid());
    T dt=grid.dX.Max()/(TV::m+1);
    ARRAY<VECTOR<STENCIL,TV::m> > stencil;
    Register_Nodes(inside_mask,stencil);
    int num_indices=fill_indices(0).max_corner;
    ARRAY<T2> du[3],tmp(num_indices);
    du[0].Resize(num_indices);
    for(int i=fill_indices(0).min_corner;i<fill_indices(0).max_corner;i++) du[0](i)=u(index_to_node(i));
    for(int o=1;o<order;o++) Fill_un(grid.one_over_dX,du[o-1],du[o],o,order);
    for(int i=0;i<order;i++) du[i](0)=T2()+FLT_MAX/100; // Sentinal values for ENO.
    tmp(0)=T2()+FLT_MAX/100;
    for(int o=order-1;o>=0;o--) for(int i=0;i<iterations;i++) Extrapolate_RK2(stencil,du[o],(o!=order-1?&du[o+1]:0),tmp,o,dt);
    for(int i=solve_indices(0).min_corner;i<solve_indices(0).max_corner;i++) u(index_to_node(i))=du[0](i);
    if(periodic)
        for(int a=0;a<TV::m;a++)
            for(NODE_ITERATOR<TV> it(grid,0,GRID<TV>::BOUNDARY_REGION,2*a+1);it.Valid();it.Next()){
                TV_INT itindex=it.index;
                itindex(a)-=combine_ends(a);
                TV_INT index=itindex;
                Periodic_Index(index);
                int k=node_to_index(index);
                if(k>0 && k<fill_indices(0).min_corner)
                    u(itindex)=u(index);}
    node_to_index.Clean_Memory();
    index_to_node.Remove_All();
    normal.Remove_All();
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Cell(std::function<bool(const TV_INT& index)> inside_mask,ARRAYS_ND_BASE<T2,TV_INT>& u)
{
    grid=grid.Get_Regular_Grid_At_MAC_Positions();
    Extrapolate_Node(inside_mask,u);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Extrapolate_Face(std::function<bool(const FACE_INDEX<TV::m>& index)> inside_mask,ARRAY<T2,FACE_INDEX<TV::m> >& u)
{
    GRID<TV> mac_grid(grid);
    for(int i=0;i<TV::m;i++){
        grid=mac_grid.Get_Face_Grid(i);
        if(periodic) combine_ends=VECTOR<bool,TV::m>::Axis_Vector(i);
        Extrapolate_Node([&](const TV_INT& index){return inside_mask(FACE_INDEX<TV::m>(i,index));},u.Component(i));}
}
//#####################################################################
// Function Lookup_Index
//#####################################################################
template<class TV,class T2> int EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Lookup_Index(TV_INT index) const
{
    Periodic_Index(index);
    if(!node_to_index.Valid_Index(index)) return 0;
    int i=node_to_index(index);
    return i>=0?i:0;
}
//#####################################################################
// Function Register_Index
//#####################################################################
template<class TV,class T2> int EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Register_Index(TV_INT index,int only_neg)
{
    Periodic_Index(index);
    if(!node_to_index.Valid_Index(index)) return 0;
    int& n=node_to_index(index);
    if(n<0 && (only_neg>=0 || n==only_neg)) n=index_to_node.Append(index);
    return n<0?0:n;
}
//#####################################################################
// Function Periodic_Index
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER<TV,T2>::
Periodic_Index(TV_INT& index) const
{
    if(!periodic) return;
    TV_INT a=index;
    TV_INT size(grid.counts-TV_INT(combine_ends));
    for(int i=0;i<TV::m;i++){
        if(index(i)<0) index(i)+=size(i);
        else if(index(i)>=size(i)) index(i)-=size(i);}
}
namespace PhysBAM{
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,1>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,2>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,3>,float>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,1>,double>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,2>,double>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,3>,double>;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,1>,SYMMETRIC_MATRIX<float,1> >;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2> >;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3> >;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,1>,SYMMETRIC_MATRIX<double,1> >;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2> >;
template class EXTRAPOLATION_HIGHER_ORDER<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3> >;
}
