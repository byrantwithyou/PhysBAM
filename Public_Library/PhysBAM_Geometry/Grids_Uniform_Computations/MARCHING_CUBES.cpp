//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
#define TRI(a,b,c,s) (((s)<<15) | ((c)<<10) | ((b)<<5) | (a))
#define TRI_ORIENT_MAP(x,f) (x?((x&0x8000) | (f[(x>>10)&31]<<10) | (f[x&31]<<5) | f[(x>>5)&31]):0)
template<> int MARCHING_CUBES_CASE<2>::edge_lookup[4][4]={{-1,0,2,-1},{0,-1,-1,3},{2,-1,-1,1},{-1,3,1,-1}};
template<> int MARCHING_CUBES_CASE<2>::vertex_lookup[4][2]={{0,1},{2,3},{0,2},{1,3}};
template<> int MARCHING_CUBES_CASE<2>::permute_map[2][8]={{2,3,0,1,4,6,5,7},{0,1,3,2,5,4,7,6}};
template<> int MARCHING_CUBES_CASE<3>::edge_lookup[8][8]={};
template<> int MARCHING_CUBES_CASE<3>::vertex_lookup[12][2]={};
template<> int MARCHING_CUBES_CASE<3>::permute_map[6][20]={};
//#####################################################################
// Function Case_Table
//#####################################################################
template<class TV> const ARRAY<MARCHING_CUBES_CASE<TV::m> >& MARCHING_CUBES<TV>::
Case_Table()
{
    static ARRAY<MARCHING_CUBES_CASE<TV::m> > table;
    static bool filled=false;
    if(!filled){
        Initialize_Case_Table(table);
        filled=true;}
    return table;
}
//#####################################################################
// Function has_ambig
//#####################################################################
static int has_ambig(int n)
{
    if((n&0x0f)==0x09 || (n&0x0f)==0x06) return true;
    if((n&0xf0)==0x90 || (n&0xf0)==0x60) return true;
    if((n&0x33)==0x21 || (n&0x33)==0x12) return true;
    if((n&0xaa)==0x82 || (n&0xaa)==0x28) return true;
    if((n&0x55)==0x41 || (n&0x55)==0x14) return true;
    if((n&0xcc)==0x84 || (n&0xcc)==0x48) return true;
    return false;
}
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Initialize_Neighbor_Cases(ARRAY<MARCHING_CUBES_CASE<TV::m> >& table, int c)
{
    for(int p=0;p<4*TV::m-6;p++){
        int b=0;
        for(int i=0;i<num_corners;i++)
            if(c&(1<<i))
                b|=1<<(MARCHING_CUBES_CASE<TV::m>::permute_map[p][i+num_edges]-num_edges);
        if(!table(b).surface[0]){
            for(int i=0;i<MARCHING_CUBES_CASE<TV::m>::max_surface;i++) table(b).surface[i]=TRI_ORIENT_MAP(table(c).surface[i],MARCHING_CUBES_CASE<TV::m>::permute_map[p]);
            for(int i=0;i<MARCHING_CUBES_CASE<TV::m>::max_boundary;i++) table(b).boundary[i]=TRI_ORIENT_MAP(table(c).boundary[i],MARCHING_CUBES_CASE<TV::m>::permute_map[p]);
            if(TV::m==3 && p>=3 && table(c).proj_dir!=p-3) table(b).proj_dir=6-table(c).proj_dir-p;
            else if(TV::m==2 && p==0) table(b).proj_dir=1-table(c).proj_dir;
            else table(b).proj_dir=table(c).proj_dir;
            table(b).enclose_inside=table(c).enclose_inside;
            Initialize_Neighbor_Cases(table, b);}}

    int b=(1<<num_corners)-c-1;
    if((TV::m==2 || !has_ambig(c)) && !table(b).surface[0]){
        table(b)=table(c);
        table(b).enclose_inside=1-table(b).enclose_inside;
        Initialize_Neighbor_Cases(table, b);}
}
static MARCHING_CUBES_CASE<2> d0 = {{}, {}, 0, 1};
static MARCHING_CUBES_CASE<2> d1 = {{TRI(0,2,0,1)}, {TRI(2,4,0,1)}, 0, 1};
static MARCHING_CUBES_CASE<2> d5 = {{TRI(0,1,0,1)}, {TRI(6,4,0,1)}, 0, 1};
static MARCHING_CUBES_CASE<2> d9 = {{TRI(0,2,0,1), TRI(1,3,0,1)}, {TRI(2,4,0,1), TRI(3,7,0,1)}, 0, 1};
static MARCHING_CUBES_CASE<2> d15 = {{}, {}, 0, 0};
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<2> >& table)
{
    table.Resize(16);
    table(0)=d0;
    table(1)=d1;
    table(5)=d5;
    table(9)=d9;
    table(15)=d15;

    Initialize_Neighbor_Cases(table,1);
    Initialize_Neighbor_Cases(table,5);
    Initialize_Neighbor_Cases(table,9);
}
static MARCHING_CUBES_CASE<3> c0 = {{}, {}, 0, 1};
static MARCHING_CUBES_CASE<3> c1 = {{TRI(8,0,4,1)}, {TRI(4,12,8,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c3 = {{TRI(9,4,8,1), TRI(4,9,5,0)}, {TRI(4,12,8,1), TRI(13,5,9,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c6 = {{TRI(9,5,0,1), TRI(1,10,4,1)}, {TRI(13,5,9,1), TRI(14,4,10,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c22 = {{TRI(9,5,0,1), TRI(6,2,8,1), TRI(1,10,4,1)}, {TRI(13,5,9,1), TRI(16,6,8,1), TRI(4,10,14,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c23 = {{TRI(1,10,6,1), TRI(2,1,6,0), TRI(5,1,2,0), TRI(9,5,2,0)}, {TRI(6,10,14,1), TRI(16,6,14,0), TRI(12,16,14,0), TRI(13,5,9,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c24 = {{TRI(6,2,8,1), TRI(1,5,11,1)}, {TRI(16,6,8,1), TRI(11,5,15,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c25 = {{TRI(6,2,4,1), TRI(4,2,0,0), TRI(1,5,11,1)}, {TRI(16,6,4,1), TRI(12,16,4,0), TRI(5,15,11,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c27 = {{TRI(6,2,4,1), TRI(4,2,1,0), TRI(1,2,9,0), TRI(11,1,9,0)}, {TRI(12,6,4,1), TRI(16,6,12,0), TRI(13,11,9,0), TRI(15,11,13,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c60 = {{TRI(6,7,8,1), TRI(8,7,9,0), TRI(10,4,5,1), TRI(10,5,11,0)}, {TRI(16,6,8,1), TRI(7,17,9,0), TRI(14,4,10,1), TRI(5,15,11,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c61 = {{TRI(6,7,9,1), TRI(0,6,9,0), TRI(5,6,0,0), TRI(11,6,5,0), TRI(10,6,11,0)}, {TRI(6,10,14,1), TRI(6,14,16,0), TRI(14,12,16,0), TRI(15,11,5,0), TRI(9,7,17,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c69 = {{TRI(6,1,3,1), TRI(8,1,6,0), TRI(0,1,8,0)}, {TRI(18,14,6,1), TRI(6,14,8,0), TRI(8,14,12,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c85 = {{TRI(0,1,2,1), TRI(2,1,3,0)}, {TRI(14,12,16,1), TRI(18,14,16,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c101 = {{TRI(6,1,3,1), TRI(8,1,6,0), TRI(0,1,8,0), TRI(2,7,9,1)}, {TRI(18,14,6,1), TRI(6,14,8,0), TRI(8,14,12,0), TRI(7,17,9,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c105 = {{TRI(4,8,0,1), TRI(7,9,2,1), TRI(11,1,5,1), TRI(10,3,6,1)}, {TRI(12,8,4,1), TRI(17,9,7,1), TRI(11,5,15,1), TRI(10,6,18,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c125 = {{TRI(5,11,3,1), TRI(5,3,0,0), TRI(0,3,9,0), TRI(9,3,7,0)}, {TRI(7,17,9,1), TRI(15,11,5,0), TRI(14,12,16,0), TRI(18,14,16,0)}, 0, 1};
static MARCHING_CUBES_CASE<3> c151 = {{TRI(1,10,6,1), TRI(2,1,6,0), TRI(5,1,2,0), TRI(9,5,2,0), TRI(3,11,7,1)}, {TRI(6,10,14,1), TRI(16,6,14,0), TRI(12,16,14,0), TRI(13,5,9,0), TRI(7,11,19,1)}, 0, 1};
static MARCHING_CUBES_CASE<3> c255 = {{}, {}, 0, 0};
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Initialize_Case_Table(ARRAY<MARCHING_CUBES_CASE<3> >& table)
{
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++)
            if(!(v&mask)){
                MARCHING_CUBES_CASE<TV::m>::vertex_lookup[k][0]=v;
                MARCHING_CUBES_CASE<TV::m>::vertex_lookup[k][1]=v|mask;
                MARCHING_CUBES_CASE<TV::m>::edge_lookup[v][v|mask]=MARCHING_CUBES_CASE<TV::m>::edge_lookup[v|mask][v]=k++;}}

    for(int a=0;a<3;a++) for(int v=0;v<8;v++) MARCHING_CUBES_CASE<TV::m>::permute_map[a][v+12]=(v^(1<<a))+12;
    for(int v=0;v<8;v++) MARCHING_CUBES_CASE<TV::m>::permute_map[3][v+12]=((v&1)|(v&2)*2|(v&4)/2)+12;
    for(int v=0;v<8;v++) MARCHING_CUBES_CASE<TV::m>::permute_map[4][v+12]=((v&1)*4|(v&2)|(v&4)/4)+12;
    for(int v=0;v<8;v++) MARCHING_CUBES_CASE<TV::m>::permute_map[5][v+12]=((v&1)*2|(v&2)/2|(v&4))+12;

    for(int i=0;i<6;i++)
        for(int e=0;e<12;e++){
            int v0=MARCHING_CUBES_CASE<TV::m>::permute_map[i][MARCHING_CUBES_CASE<TV::m>::vertex_lookup[e][0]+12]-12;
            int v1=MARCHING_CUBES_CASE<TV::m>::permute_map[i][MARCHING_CUBES_CASE<TV::m>::vertex_lookup[e][1]+12]-12;
            MARCHING_CUBES_CASE<TV::m>::permute_map[i][e]=MARCHING_CUBES_CASE<TV::m>::edge_lookup[v0][v1];}

    table.Resize(256);
    table(0)=c0;
    table(1)=c1;
    table(3)=c3;
    table(6)=c6;
    table(22)=c22;
    table(23)=c23;
    table(24)=c24;
    table(25)=c25;
    table(27)=c27;
    table(60)=c60;
    table(61)=c61;
    table(69)=c69;
    table(85)=c85;
    table(101)=c101;
    table(105)=c105;
    table(125)=c125;
    table(151)=c151;
    table(255)=c255;

    Initialize_Neighbor_Cases(table,1);
    Initialize_Neighbor_Cases(table,3);
    Initialize_Neighbor_Cases(table,6);
    Initialize_Neighbor_Cases(table,22);
    Initialize_Neighbor_Cases(table,23);
    Initialize_Neighbor_Cases(table,24);
    Initialize_Neighbor_Cases(table,25);
    Initialize_Neighbor_Cases(table,27);
    Initialize_Neighbor_Cases(table,60);
    Initialize_Neighbor_Cases(table,61);
    Initialize_Neighbor_Cases(table,69);
    Initialize_Neighbor_Cases(table,85);
    Initialize_Neighbor_Cases(table,101);
    Initialize_Neighbor_Cases(table,105);
    Initialize_Neighbor_Cases(table,125);
    Initialize_Neighbor_Cases(table,151);
}
//#####################################################################
// Function Create_Surface
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Create_Surface(T_SURFACE& surface,const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi)
{
    VECTOR<TV_INT,num_corners> bits=GRID<TV>::Binary_Counts(TV_INT());
    VECTOR<T,num_corners> phis;
    TV pts[num_pts];
    const ARRAY<MARCHING_CUBES_CASE<TV::m> >& table=Case_Table();
    HASHTABLE<FACE_INDEX<TV::m>,int> ht;

    for(UNIFORM_ARRAY_ITERATOR<TV::m> it(phi.domain.To_Closed());it.Valid();it.Next()){
        int c=0;
        for(int i=0;i<bits.m;i++){
            TV_INT index=it.index+bits(i);
            pts[i+num_edges]=grid.Node(index);
            phis(i)=phi(index);
            c|=(phis(i)<0)<<i;}
        for(int i=0;i<num_edges;i++){
            int v0=MARCHING_CUBES_CASE<TV::m>::vertex_lookup[i][0];
            int v1=MARCHING_CUBES_CASE<TV::m>::vertex_lookup[i][1];
            if(((c>>v0)^(c>>v1))&1){
                T t=phis(v0)/(phis(v0)-phis(v1));
                pts[i]=pts[v0+num_edges]+t*(pts[v1+num_edges]-pts[v0+num_edges]);}}

        const MARCHING_CUBES_CASE<TV::m> & cs=table(c);

        for(int i=0;i<MARCHING_CUBES_CASE<TV::m>::max_surface && cs.surface[i];i++){
            TV_INT face;
            for(int j=0;j<TV::m;j++){
                int e=(cs.surface[i]>>5*j)&31;
                FACE_INDEX<TV::m> fi(e/(1<<(TV::m-1)),it.index+bits(MARCHING_CUBES_CASE<TV::m>::vertex_lookup[e][0]));
                if(!ht.Get(fi,face(j))){
                    int index=surface.particles.Add_Element();
                    face(j)=index;
                    ht.Set(fi,index);
                    surface.particles.X(index)=pts[e];}}
            if(!cs.enclose_inside) exchange(face.x,face.y);
            surface.mesh.elements.Append(face);}}

    surface.Update_Number_Nodes();
}
//#####################################################################
// Function Get_Triangles_For_Cell
//#####################################################################
template<class TV> void MARCHING_CUBES<TV>::
Get_Elements_For_Cell(ARRAY<T_FACE>& surface,ARRAY<T_FACE>& boundary,int& direction,bool& enclose_inside,const ARRAY<T,TV_INT>& phi,const TV_INT& cell)
{
    VECTOR<TV_INT,num_corners> bits=GRID<TV>::Binary_Counts(TV_INT());
    VECTOR<T,num_corners> phis;
    TV pts[num_pts];
    const ARRAY<MARCHING_CUBES_CASE<TV::m> >& table=Case_Table();

    int c=0;
    for(int i=0;i<bits.m;i++){
        pts[i+num_edges]=TV(bits(i));
        phis(i)=phi(cell+bits(i));
        c|=(phis(i)<0)<<i;}

    for(int i=0;i<num_edges;i++){
        int v0=MARCHING_CUBES_CASE<TV::m>::vertex_lookup[i][0];
        int v1=MARCHING_CUBES_CASE<TV::m>::vertex_lookup[i][1];
        if(((c>>v0)^(c>>v1))&1){
            T t=phis(v0)/(phis(v0)-phis(v1));
            pts[i]=pts[v0+num_edges]+t*(pts[v1+num_edges]-pts[v0+num_edges]);}}

    const MARCHING_CUBES_CASE<TV::m> & cs=table(c);

    int len[2]={MARCHING_CUBES_CASE<TV::m>::max_surface,MARCHING_CUBES_CASE<TV::m>::max_boundary};
    const unsigned short* elements[2]={cs.surface,cs.boundary};
    ARRAY<T_FACE>* list[2]={&surface,&boundary};
    for(int s=0;s<2;s++){
        for(int i=0;i<len[s] && elements[s][i];i++){
            T_FACE face;
            for(int j=0;j<TV::m;j++){
                int e=(elements[s][i]>>5*j)&31;
                face.X(j)=pts[e];}
            list[s]->Append(face);}}
    direction=cs.proj_dir;
    enclose_inside=cs.enclose_inside;
}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void MARCHING_CUBES<VECTOR<float,2> >::Get_Elements_For_Cell(ARRAY<SEGMENT_2D<float>,int>&,ARRAY<SEGMENT_2D<float>,int>&,int&,bool&,ARRAY<float,VECTOR<int,2> > const&,VECTOR<int,2> const&);
template const ARRAY<MARCHING_CUBES_CASE<2> >& MARCHING_CUBES<VECTOR<float,2> >::Case_Table();
template const ARRAY<MARCHING_CUBES_CASE<3> >& MARCHING_CUBES<VECTOR<float,3> >::Case_Table();
template void MARCHING_CUBES<VECTOR<float,3> >::Create_Surface(TRIANGULATED_SURFACE<float>&,GRID<VECTOR<float,3> > const&,ARRAY<float,VECTOR<int,3> > const&);
template void MARCHING_CUBES<VECTOR<float,3> >::Get_Elements_For_Cell(ARRAY<TRIANGLE_3D<float>,int>&,ARRAY<TRIANGLE_3D<float>,int>&,int&,bool&,ARRAY<float,VECTOR<int,3> > const&,VECTOR<int,3> const&);
template void MARCHING_CUBES<VECTOR<float,2> >::Create_Surface(SEGMENTED_CURVE_2D<float>&,GRID<VECTOR<float,2> > const&,ARRAY<float,VECTOR<int,2> > const&);
#endif
template void MARCHING_CUBES<VECTOR<double,2> >::Get_Elements_For_Cell(ARRAY<SEGMENT_2D<double>,int>&,ARRAY<SEGMENT_2D<double>,int>&,int&,bool&,ARRAY<double,VECTOR<int,2> > const&,VECTOR<int,2> const&);
template const ARRAY<MARCHING_CUBES_CASE<2> >& MARCHING_CUBES<VECTOR<double,2> >::Case_Table();
template const ARRAY<MARCHING_CUBES_CASE<3> >& MARCHING_CUBES<VECTOR<double,3> >::Case_Table();
template void MARCHING_CUBES<VECTOR<double,3> >::Create_Surface(TRIANGULATED_SURFACE<double>&,GRID<VECTOR<double,3> > const&,ARRAY<double,VECTOR<int,3> > const&);
template void MARCHING_CUBES<VECTOR<double,3> >::Get_Elements_For_Cell(ARRAY<TRIANGLE_3D<double>,int>&,ARRAY<TRIANGLE_3D<double>,int>&,int&,bool&,ARRAY<double,VECTOR<int,3> > const&,VECTOR<int,3> const&);
template void MARCHING_CUBES<VECTOR<double,2> >::Create_Surface(SEGMENTED_CURVE_2D<double>&,GRID<VECTOR<double,2> > const&,ARRAY<double,VECTOR<int,2> > const&);
