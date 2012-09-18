//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include "MARCHING_TETRAHEDRA_CUTTING.h"
using namespace PhysBAM;
const int edge_table[2][6][2]=
{
    {{1,2},{0,2},{0,1}},
    {{0,2},{0,1},{1,2},{0,3},{2,3},{1,3}}
};
//#####################################################################
// Function Case_Table
//#####################################################################
template<class TV> const ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m> >& MARCHING_TETRAHEDRA_CUTTING<TV>::
Case_Table()
{
    static ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m> > table;
    static bool filled=false;
    if(!filled){
        Initialize_Case_Table(table);
        filled=true;}
    return table;
}
#define A(i,j) cc.elements[i][j]%16
#define B(i,j) (cc.elements[i][j]>>4)%16
#define C(i,j) (cc.elements[i][j]>>8)%16
#define D(i,j) (cc.elements[i][j]>>12)
int depth=0;
static void Fill_Helper(ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<2> >& table,const MARCHING_TETRAHEDRA_CUTTING_CASE<2>& cc,int cs)
{
    if(table(cs).elements[0][0] || table(cs).elements[1][0]){
        printf("%*scase filled %i%i%i\n",2*depth,"",cs/9,(cs/3)%3,cs%3);
        return;}
    table(cs)=cc;
    printf("%*sfill %i%i%i { %03x %03x } { %03x %03x }\n",2*depth,"",cs/9,(cs/3)%3,cs%3,cc.elements[0][0],cc.elements[0][1],cc.elements[1][0],cc.elements[1][1]);
    depth++;

    int a=cs/9,b=(cs/3)%3,c=cs%3;
    // sign flip
    MARCHING_TETRAHEDRA_CUTTING_CASE<2> cc1={{0,0,0,0},{{cc.elements[1][0],cc.elements[1][1]},
                                                        {cc.elements[0][0],cc.elements[0][1]}}};
    Fill_Helper(table,cc1,(a^(1-a/2))*9+(b^(1-b/2))*3+(c^(1-c/2)));
#define F(i,j,m) (cc.elements[i][j]?(m[B(i,j)]+(m[A(i,j)]<<4)+(m[C(i,j)]<<8)):0)
    // flip
    const int fl[6]={1,0,2,4,3,5};
    printf("%03x %03x\n", cc.elements[0][0], F(0,0,fl));
    MARCHING_TETRAHEDRA_CUTTING_CASE<2> cc2={{0,0,0,0},{{F(0,0,fl),F(0,1,fl)},{F(1,0,fl),F(1,1,fl)}}};
    Fill_Helper(table,cc2,b*9+a*3+c);
#undef F
#define R(i,j,m) (cc.elements[i][j]?(m[A(i,j)]+(m[B(i,j)]<<4)+(m[C(i,j)]<<8)):0)
    // rotate
    const int ro[6]={1,2,0,4,5,3};
    MARCHING_TETRAHEDRA_CUTTING_CASE<2> cc0={{0,0,0,0},{{R(0,0,ro),R(0,1,ro)},{R(1,0,ro),R(1,1,ro)}}};
    Fill_Helper(table,cc0,c*9+a*3+b);
#undef R
    depth--;
}
template<class TV> void MARCHING_TETRAHEDRA_CUTTING<TV>::
Initialize_Case_Table(ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<2> >& table)
{
    typedef MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m> CASE;
#define EE(x,y,z) (((z)<<8)+((y)<<4)+(x))
    CASE c0={{0,0,0,0},{{EE(0,1,2),0},{0,0}}};
    CASE c1={{0,0,0,0},{{EE(0,1,3),EE(3,4,0)},{EE(4,3,2),0}}};
    CASE c5={{0,0,0,0},{{EE(0,5,2),0},{EE(5,1,2),0}}};
#undef EE
    table.Resize(27);
    Fill_Helper(table,c0,0);
    Fill_Helper(table,c1,1);
    Fill_Helper(table,c0,2);
    Fill_Helper(table,c5,5);
    Fill_Helper(table,c0,8);
    Fill_Helper(table,c0,26);
}
static void Fill_Helper(ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<3> >& table,const MARCHING_TETRAHEDRA_CUTTING_CASE<3>& cc,int cs)
{
    if(table(cs).elements[0][0] || table(cs).elements[1][0]) return;
    table(cs)=cc;

    int a=cs/27,b=(cs/9)%3,c=(cs/3)%3,d=cs%3;
    // sign flip
    MARCHING_TETRAHEDRA_CUTTING_CASE<3> cc1={
        {cc.comp[0],cc.comp[1],cc.comp[2],cc.comp[3]},
        {{cc.elements[1][0],cc.elements[1][1]},{cc.elements[0][0],cc.elements[0][1]}}};
    Fill_Helper(table,cc1,(a^(1-a/2))*27+(b^(1-b/2))*9+(c^(1-c/2))*3+(d^(1-d/2)));
#define F(i,j,m) (cc.elements[i][j]?(m[B(i,j)]+(m[A(i,j)]<<4)+(m[C(i,j)]<<8)+(m[D(i,j)]<<12)):0)
    // flip
    const int fl[10]={1,0,2,3,6,5,4,9,8,7};
    MARCHING_TETRAHEDRA_CUTTING_CASE<3> cc2={
        {fl[cc.comp[0]],fl[cc.comp[1]],fl[cc.comp[2]],fl[cc.comp[3]]},
        {{F(0,0,fl),F(0,1,fl)},{F(1,0,fl),F(1,1,fl)}}};
    Fill_Helper(table,cc2,b*27+a*9+c*3+d);
    // rotate
    const int ro[10]={1,2,3,0,9,6,8,5,7,4};
    MARCHING_TETRAHEDRA_CUTTING_CASE<3> cc0={
        {ro[cc.comp[0]],ro[cc.comp[1]],ro[cc.comp[2]],ro[cc.comp[3]]},
        {{F(0,0,ro),F(0,1,ro)},{F(1,0,ro),F(1,1,ro)}}};
    Fill_Helper(table,cc0,d*27+a*9+b*3+c);
#undef R
}
template<class TV> void MARCHING_TETRAHEDRA_CUTTING<TV>::
Initialize_Case_Table(ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<3> >& table)
{
    typedef MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m> CASE;
#define EE(w,x,y,z) (((z)<<12)+((y)<<8)+((x)<<4)+(w))
    CASE c0={{0,0,0,0},{{EE(0,1,2,3),0,0},{0,0,0}}};
    CASE c1={{1,2,2,0},{{EE(8,1,9,7),EE(8,2,1,7),EE(0,1,2,7)},{EE(7,9,8,3),0,0}}};
    CASE c4={{2,0,3,1},{{EE(0,5,2,7),EE(6,5,7,2),EE(6,7,8,2)},{EE(7,6,8,3),EE(7,5,6,3),EE(3,5,6,1)}}};
    CASE c5={{2,0,0,0},{{EE(0,5,2,3),EE(2,5,6,3),0},{EE(5,1,6,3),0,0}}};
    CASE c17={{0,0,0,0},{{EE(4,1,2,3),0,0},{EE(0,1,4,3),0,0}}};
#undef EE
    table.Resize(81);
    Fill_Helper(table,c0,0);
    Fill_Helper(table,c1,1);
    Fill_Helper(table,c0,2);
    Fill_Helper(table,c4,4);
    Fill_Helper(table,c5,5);
    Fill_Helper(table,c0,8);
    Fill_Helper(table,c17,17);
    Fill_Helper(table,c0,26);
    Fill_Helper(table,c0,80);
}
//#####################################################################
// Function Query_Case
//#####################################################################
template<class TV> void MARCHING_TETRAHEDRA_CUTTING<TV>::
Query_Case(ARRAY<E>& parents,ARRAY<E>& children,ARRAY<E>& split_parents,const ARRAY<T>& phi,ARRAY<PAIR<S,T> >& weights)
{
    ARRAY<E> new_parents,new_children;
    ARRAY<bool> side;
    const ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m> >& table=Case_Table();
    HASHTABLE<S,int> edge_hash;
    HASHTABLE<PAIR<TV_INT,int>,int> neighbor;
    UNION_FIND<> union_find(2*max_side_elements*(TV::m+1)*parents.m);
    int next_parent=0;
    for(int i=0;i<children.m;i++){
        E e=children(i);
        VECTOR<T,TV::m+1> p(phi.Subset(e));
        int cs=0;
        for(int j=0;j<TV::m+1;j++) cs=cs*3+(p(j)==0?2:p(j)>0);
        printf("case %i\n", cs);
        const MARCHING_TETRAHEDRA_CUTTING_CASE<TV::m>& mtcc=table(cs);
        int inv=0;
        if(e(mtcc.comp[0])>e(mtcc.comp[1])){
            exchange(e(mtcc.comp[0]),e(mtcc.comp[1]));
            inv^=1;}
        if(e(mtcc.comp[2])>e(mtcc.comp[3])){
            exchange(e(mtcc.comp[2]),e(mtcc.comp[3]));
            inv^=1;}
        if(e(mtcc.comp[0])>e(mtcc.comp[1])){
            exchange(e(mtcc.comp[0]),e(mtcc.comp[1]));
            inv^=1;}
        for(int s=0;s<2;s++)
            for(int t=0;t<max_side_elements;t++){
                int mask=mtcc.elements[s][t];
                if(!mask) break;
                printf("mask %03x\n", mask);
                E elem,par;
                for(int j=0;j<TV::m+1;j++){
                    par(j)=next_parent++;
                    int ed=(mask>>(4*j))%16;
                    if(ed>=TV::m+1){
                        ed-=TV::m+1;
                        printf("%i -> %i %i\n", ed, edge_table[TV::m-2][ed][0], edge_table[TV::m-2][ed][1]);
                        S seg(e(edge_table[TV::m-2][ed][0]),e(edge_table[TV::m-2][ed][1]));
                        if(!edge_hash.Get(seg,ed)){
                            ed=weights.Append(PAIR<S,T>(seg,phi(seg.x)/(phi(seg.x)-phi(seg.y))))+phi.m;
                            edge_hash.Set(seg,ed);}}
                    else ed=e(ed);
                    elem(j)=ed;}
                if(inv) exchange(elem(0),elem(1));
                LOG::cout<<e<<elem<<std::endl;
                new_parents.Append(parents(i));
                split_parents.Append(par);
                new_children.Append(elem);

                for(int j=0;j<TV::m+1;j++){
                    int other=-1;
                    PAIR<TV_INT,int> pr(elem.Remove_Index(j).Sorted(),s);
                    LOG::cout<<pr<<std::endl;
                    if(neighbor.Get(pr,other)){
                        LOG::cout<<parents(i)<<"  "<<new_parents(other)<<"  "<<par<<"  "<<split_parents(other)<<std::endl;
                        for(int a=0;a<TV::m+1;a++)
                            for(int b=0;b<TV::m+1;b++)
                                if(parents(i)(a)==new_parents(other)(b))
                                    union_find.Union(par(a),split_parents(other)(b));}
                    else neighbor.Set(pr,split_parents.m-1);}}}

    ARRAY<int> index_map(split_parents.m*(TV::m+1));
    index_map.Fill(-1);
    int condense=0;
    for(int i=0;i<index_map.m;i++){
        int p=union_find.Find(i);
        if(index_map(p)<0) index_map(p)=index_map(i)=condense++;}
    LOG::cout<<index_map<<std::endl;
    split_parents.Flattened()=index_map.Subset(split_parents.Flattened());
    new_parents.Exchange(parents);
    new_children.Exchange(children);
}
template const ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<2> >& MARCHING_TETRAHEDRA_CUTTING<VECTOR<float,2> >::Case_Table();
template const ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<3> >& MARCHING_TETRAHEDRA_CUTTING<VECTOR<float,3> >::Case_Table();
template void MARCHING_TETRAHEDRA_CUTTING<VECTOR<float,3> >::Query_Case(ARRAY<VECTOR<int,4>,int>&,
    ARRAY<VECTOR<int,4>,int>&,ARRAY<VECTOR<int,4>,int>&,ARRAY<float,int> const&,ARRAY<PAIR<VECTOR<int,2>,float>,int>&);
template void MARCHING_TETRAHEDRA_CUTTING<VECTOR<float,2> >::Query_Case(ARRAY<VECTOR<int,3>,int>&,
    ARRAY<VECTOR<int,3>,int>&,ARRAY<VECTOR<int,3>,int>&,ARRAY<float,int> const&,ARRAY<PAIR<VECTOR<int,2>,float>,int>&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template const ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<2> >& MARCHING_TETRAHEDRA_CUTTING<VECTOR<double,2> >::Case_Table();
template const ARRAY<MARCHING_TETRAHEDRA_CUTTING_CASE<3> >& MARCHING_TETRAHEDRA_CUTTING<VECTOR<double,3> >::Case_Table();
template void MARCHING_TETRAHEDRA_CUTTING<VECTOR<double,3> >::Query_Case(ARRAY<VECTOR<int,4>,int>&,
    ARRAY<VECTOR<int,4>,int>&,ARRAY<VECTOR<int,4>,int>&,ARRAY<double,int> const&,ARRAY<PAIR<VECTOR<int,2>,double>,int>&);
template void MARCHING_TETRAHEDRA_CUTTING<VECTOR<double,2> >::Query_Case(ARRAY<VECTOR<int,3>,int>&,
    ARRAY<VECTOR<int,3>,int>&,ARRAY<VECTOR<int,3>,int>&,ARRAY<double,int> const&,ARRAY<PAIR<VECTOR<int,2>,double>,int>&);
#endif
