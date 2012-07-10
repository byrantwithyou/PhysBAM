//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION_HIGHER_ORDER_POLY.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T2>::
EXTRAPOLATION_HIGHER_ORDER_POLY()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2> EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T2>::
~EXTRAPOLATION_HIGHER_ORDER_POLY()
{
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T2>::
Extrapolate_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<bool,TV::m> >& inside_mask,int ghost,ARRAYS_ND_BASE<VECTOR<T2,TV::m> >& u,int order,int fill_width,T order_reduction_penalty)
{
    PHYSBAM_ASSERT(fill_width<=ghost);
    ARRAY<int,TV_INT> distance(grid.Domain_Indices(ghost+1));
    ARRAY<TV_INT,TV_INT> nearest(grid.Domain_Indices(ghost+1));
    ARRAY<ARRAY<TV_INT> > todo;
    int inf=INT_MAX/10,ignore=-inf;

    ARRAY<TV_INT> neighbors;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT()-1,TV_INT()+2));it.Valid();it.Next())
        neighbors.Append(it.index);

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        if(inside_mask(it.index)) continue;
        distance(it.index)|=1;
        for(int d=0;d<TV::m;d++){
            TV_INT ia(it.index),ib(it.index);
            ia(d)--;ib(d)++;
            distance(ia)|=2;
            distance(ib)|=2;}
        for(int d=0;d<TV::m;d++)
            for(int e=d+1;e<TV::m;e++){
                TV_INT ia(it.index),ib(it.index),ic(it.index),id(it.index);
                ia(d)--;ib(d)--;ic(d)++;id(d)++;ia(e)--;ib(e)++;ic(e)--;id(e)++;
                distance(ia)|=4;
                distance(ib)|=4;
                distance(ic)|=4;
                distance(id)|=4;}}

    VECTOR<ARRAY<TV_INT>,3> seed;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        int& a=distance(it.index),oo=-1;
        for(int o=0,m=1;o<order;o++,m|=1<<o)
            if((a&m)==0)
                oo=o;
        if(oo==-1){
            a=ignore;
            continue;}
        a=~oo;
        seed(oo).Append(it.index);}

    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost+1,GRID<TV>::GHOST_REGION);it.Valid();it.Next())
        distance(it.index)=ignore;

    ARRAY<TV_INT> next_outside,current,solve;
    for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid);it.Valid();it.Next()){
        int d=distance(it.index);
        if(d==ignore || d==inf) continue;
        for(int i=0;i<neighbors.m;i++){
            TV_INT neighbor_index=it.index+neighbors(i);
            int& a=distance(neighbor_index);
            if(a==ignore){
                a=inf;
                next_outside.Append(neighbor_index);}}}

    for(int o=1;o<fill_width;o++){
        solve.Append_Elements(next_outside);
        current.Exchange(next_outside);
        for(int j=0;j<current.m;j++)
            for(int i=0;i<neighbors.m;i++){
                TV_INT neighbor_index=current(j)+neighbors(i);
                int& a=distance(neighbor_index);
                if(a==ignore){
                    a=inf;
                    next_outside.Append(neighbor_index);}}
        current.Remove_All();}
    solve.Append_Elements(next_outside);
    next_outside.Remove_All();

    for(int o=0;o<order;o++){
        todo.Resize(1);
        for(int i=0;i<seed(o).m;i++){
            distance(seed(o)(i))=0;
            nearest(seed(o)(i))=seed(o)(i);
            todo(0).Append(seed(o)(i));}

        for(int l=0;l<todo.m;l++){
            for(int i=0;i<todo(l).m;i++){
                TV_INT index=todo(l)(i),pt=nearest(index);
                if((index-pt).Magnitude_Squared()!=distance(index)) continue;
                for(int j=0;j<neighbors.m;j++){
                    TV_INT index2=neighbors(j)+index;
                    int nd=(index2-pt).Magnitude_Squared();
                    if(nd<distance(index2)){
                        distance(index2)=nd;
                        nearest(index2)=pt;
                        if(todo.m<=nd) todo.Resize(nd+1);
                        todo(nd).Append(index2);}}}
            todo(l).Clean_Memory();}
        todo.Clean_Memory();

        if(o!=order-1)
            for(UNIFORM_GRID_ITERATOR_NODE<TV> it(grid,ghost+1);it.Valid();it.Next()){
                int& a=distance(it.index);
                if(a>=0 && a<inf) a=(int)sqr(sqrt(a)+order_reduction_penalty);}}

    for(int i=0;i<solve.m;i++){
        TV_INT index=solve(i);
        TV_INT pt=nearest(index);
        bool boundary=false;
        T& uu=u(index),u0=u(pt);
        uu=u0;
        if(order<=1) continue;
        for(int d=0;d<TV::m;d++)
            if(pt(d)==0 || pt(d)==grid.counts(d)-1){
                boundary=true;
                break;}
        if(boundary) continue;
        TV_INT dx=index-pt;
        for(int d=0;d<TV::m;d++){
            TV_INT ia(pt),ib(pt);
            ia(d)--;ib(d)++;
            if(!inside_mask(ia) || !inside_mask(ib)) continue;
            T ua=u(ia),ub=u(ib);
            uu+=dx(d)*(ub-ua)/2;
            if(order>=3) uu+=sqr(dx(d))*(ub-2*u0+ua)/2;}
        if(order<=2) continue;
        for(int d=0;d<TV::m;d++)
            for(int e=d+1;e<TV::m;e++){
                TV_INT ia(pt),ib(pt),ic(pt),id(pt);
                ia(d)--;ib(d)--;ic(d)++;id(d)++;ia(e)--;ib(e)++;ic(e)--;id(e)++;
                if(!inside_mask(ia) || !inside_mask(ib) || !inside_mask(ic) || !inside_mask(id)) continue;
                uu+=dx(d)*dx(e)*(u(ia)-u(ib)-u(ic)+u(id))/4;}}
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T2>::
Extrapolate_Cell(const GRID<TV>& grid,const ARRAYS_ND_BASE<VECTOR<bool,TV::m> >& inside_mask,
    int ghost,ARRAYS_ND_BASE<VECTOR<T2,TV::m> >& u,int order,int fill_width,T order_reduction_penalty)
{
    GRID<TV> node_grid(grid.Get_Regular_Grid_At_MAC_Positions());
    Extrapolate_Node(node_grid,inside_mask,ghost,u,order,fill_width,order_reduction_penalty);
}
//#####################################################################
// Function Quadratic_Extrapolate
//#####################################################################
template<class TV,class T2> void EXTRAPOLATION_HIGHER_ORDER_POLY<TV,T2>::
Extrapolate_Face(const GRID<TV>& grid,const ARRAY<bool,FACE_INDEX<TV::m> >& inside_mask,
    int ghost,ARRAY<T2,FACE_INDEX<TV::m> >& u,int order,int fill_width,T order_reduction_penalty)
{
    for(int i=0;i<TV::m;i++){
        GRID<TV> node_grid(grid.Get_Face_Grid(i));
        Extrapolate_Node(node_grid,inside_mask.Component(i),ghost,u.Component(i),order,fill_width,order_reduction_penalty);}
}
template class EXTRAPOLATION_HIGHER_ORDER_POLY<VECTOR<float,1>,float>;
template class EXTRAPOLATION_HIGHER_ORDER_POLY<VECTOR<float,2>,float>;
template class EXTRAPOLATION_HIGHER_ORDER_POLY<VECTOR<float,3>,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXTRAPOLATION_HIGHER_ORDER_POLY<VECTOR<double,1>,double>;
template class EXTRAPOLATION_HIGHER_ORDER_POLY<VECTOR<double,2>,double>;
template class EXTRAPOLATION_HIGHER_ORDER_POLY<VECTOR<double,3>,double>;
#endif
