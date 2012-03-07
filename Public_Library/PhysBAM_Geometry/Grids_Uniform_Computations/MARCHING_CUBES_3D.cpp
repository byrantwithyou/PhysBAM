//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/UNIFORM_ARRAY_ITERATOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
static int edge_lookup[8][8];
static int vertex_lookup[12][2];
static int permute_map[6][20];
#define TRI(a,b,c,s) (((s)<<15) | ((c)<<10) | ((b)<<5) | (a))
#define TRI_ORIENT_MAP(x,f) (x?((x&0x8000) | (f[(x>>10)&31]<<10) | (f[x&31]<<5) | f[(x>>5)&31]):0)
//#####################################################################
// Function Case_Table
//#####################################################################
template<class T> const ARRAY<MARCHING_CUBES_3D_CASE>& MARCHING_CUBES_3D<T>::
Case_Table()
{
    static ARRAY<MARCHING_CUBES_3D_CASE> table;
    static bool filled=false;
    if(!filled){
        Initialize_Case_Table(table);
        filled=true;}
    return table;
}
//#####################################################################
// Function flip
//#####################################################################
static int permute_case(int c,int* f)
{
    int r=0;
    for(int i=0;i<8;i++)
        if(c&(1<<i))
            r|=1<<(f[i+12]-12);
    return r;
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
template<class T> void MARCHING_CUBES_3D<T>::
Initialize_Neighbor_Cases(ARRAY<MARCHING_CUBES_3D_CASE>& table, int c)
{
    for(int p=0;p<6;p++){
        int b=permute_case(c, permute_map[p]);
        if(!table(b).elements[0]){
            for(int i=0;i<MARCHING_CUBES_3D_CASE::max_elements;i++) table(b).elements[i]=TRI_ORIENT_MAP(table(c).elements[i],permute_map[p]);
            for(int i=0;i<MARCHING_CUBES_3D_CASE::sheet_elements;i++) table(b).boundary[i]=TRI_ORIENT_MAP(table(c).boundary[i],permute_map[p]);
            if(p>=3 && table(c).proj_dir!=p-3) table(b).proj_dir=6-table(c).proj_dir-p;
            else table(b).proj_dir=table(c).proj_dir;
            table(b).enclose_inside=table(c).enclose_inside;
            Initialize_Neighbor_Cases(table, b);}}

    if(!has_ambig(c)){
        int b=255-c;
        if(!table(b).elements[0]){
            table(b)=table(c);
            table(b).enclose_inside=1-table(b).enclose_inside;
            Initialize_Neighbor_Cases(table, b);}}
}
static MARCHING_CUBES_3D_CASE c0 = {{}, {}, 0, 1};
static MARCHING_CUBES_3D_CASE c1 = {{TRI(8,0,4,1)}, {TRI(4,12,8,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c3 = {{TRI(9,4,8,1), TRI(4,9,5,0)}, {TRI(4,12,8,1), TRI(13,5,9,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c6 = {{TRI(9,5,0,1), TRI(1,10,4,1)}, {TRI(13,5,9,1), TRI(14,4,10,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c22 = {{TRI(9,5,0,1), TRI(6,2,8,1), TRI(1,10,4,1)}, {TRI(13,5,9,1), TRI(16,6,8,1), TRI(4,10,14,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c23 = {{TRI(1,10,6,1), TRI(2,1,6,0), TRI(5,1,2,0), TRI(9,5,2,0)}, {TRI(6,10,14,1), TRI(16,6,14,0), TRI(12,16,14,0), TRI(13,5,9,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c24 = {{TRI(6,2,8,1), TRI(1,5,11,1)}, {TRI(16,6,8,1), TRI(11,5,15,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c25 = {{TRI(6,2,4,1), TRI(4,2,0,0), TRI(1,5,11,1)}, {TRI(16,6,4,1), TRI(12,16,4,0), TRI(5,15,11,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c27 = {{TRI(6,2,4,1), TRI(4,2,1,0), TRI(1,2,9,0), TRI(11,1,9,0)}, {TRI(12,6,4,1), TRI(16,6,12,0), TRI(13,11,9,0), TRI(15,11,13,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c60 = {{TRI(6,7,8,1), TRI(8,7,9,0), TRI(10,4,5,1), TRI(10,5,11,0)}, {TRI(16,6,8,1), TRI(7,17,9,0), TRI(14,4,10,1), TRI(5,15,11,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c61 = {{TRI(6,7,9,1), TRI(0,6,9,0), TRI(5,6,0,0), TRI(11,6,5,0), TRI(10,6,11,0)}, {TRI(6,10,14,1), TRI(6,14,16,0), TRI(14,12,16,0), TRI(15,11,5,0), TRI(9,7,17,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c69 = {{TRI(6,1,3,1), TRI(8,1,6,0), TRI(0,1,8,0)}, {TRI(18,14,6,1), TRI(6,14,8,0), TRI(8,14,12,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c85 = {{TRI(0,1,2,1), TRI(2,1,3,0)}, {TRI(14,12,16,1), TRI(18,14,16,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c101 = {{TRI(6,1,3,1), TRI(8,1,6,0), TRI(0,1,8,0), TRI(2,7,9,1)}, {TRI(18,14,6,1), TRI(6,14,8,0), TRI(8,14,12,0), TRI(7,17,9,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c105 = {{TRI(4,8,0,1), TRI(7,9,2,1), TRI(11,1,5,1), TRI(10,3,6,1)}, {TRI(12,8,4,1), TRI(17,9,7,1), TRI(11,5,15,1), TRI(10,6,18,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c125 = {{TRI(5,11,3,1), TRI(5,3,0,0), TRI(0,3,9,0), TRI(9,3,7,0)}, {TRI(7,17,9,1), TRI(15,11,5,0), TRI(14,12,16,0), TRI(18,14,16,0)}, 0, 1};
static MARCHING_CUBES_3D_CASE c151 = {{TRI(1,10,6,1), TRI(2,1,6,0), TRI(5,1,2,0), TRI(9,5,2,0), TRI(3,11,7,1)}, {TRI(6,10,14,1), TRI(16,6,14,0), TRI(12,16,14,0), TRI(13,5,9,0), TRI(7,11,19,1)}, 0, 1};
static MARCHING_CUBES_3D_CASE c255 = {{}, {}, 0, 0};
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<class T> void MARCHING_CUBES_3D<T>::
Initialize_Case_Table(ARRAY<MARCHING_CUBES_3D_CASE>& table)
{
    for(int a=0,k=0;a<3;a++){
        int mask=1<<a;
        for(int v=0;v<8;v++)
            if(!(v&mask)){
                vertex_lookup[k][0]=v;
                vertex_lookup[k][1]=v|mask;
                edge_lookup[v][v|mask]=edge_lookup[v|mask][v]=k++;}}

    for(int a=0;a<3;a++) for(int v=0;v<8;v++) permute_map[a][v+12]=(v^(1<<a))+12;
    for(int v=0;v<8;v++) permute_map[3][v+12]=((v&1)|(v&2)*2|(v&4)/2)+12;
    for(int v=0;v<8;v++) permute_map[4][v+12]=((v&1)*4|(v&2)|(v&4)/4)+12;
    for(int v=0;v<8;v++) permute_map[5][v+12]=((v&1)*2|(v&2)/2|(v&4))+12;

    for(int i=0;i<6;i++)
        for(int e=0;e<12;e++){
            int v0=permute_map[i][vertex_lookup[e][0]+12]-12;
            int v1=permute_map[i][vertex_lookup[e][1]+12]-12;
            permute_map[i][e]=edge_lookup[v0][v1];}

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
template<class T> void MARCHING_CUBES_3D<T>::
Create_Surface(TRIANGULATED_SURFACE<T>& surface,const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi)
{
    VECTOR<TV_INT,8> bits=GRID<TV>::Binary_Counts(TV_INT());

    VECTOR<T,8> phis;
    TV pts[20];

    const ARRAY<MARCHING_CUBES_3D_CASE>& table=Case_Table();
    
    HASHTABLE<FACE_INDEX<TV::m>,int> ht;

    for(UNIFORM_ARRAY_ITERATOR<TV::m> it(phi.domain.To_Closed());it.Valid();it.Next()){
        int c=0;
        for(int i=0;i<bits.m;i++){
            TV_INT ind=it.index+bits(i);
            pts[i+12]=grid.Node(ind);
            phis(i)=phi(ind);
            c|=(phis(i)<0)<<i;}
        for(int i=0;i<12;i++){
            int v0=vertex_lookup[i][0];
            int v1=vertex_lookup[i][1];
            if(((c>>v0)^(c>>v1))&1){
                T t=phis(v0)/(phis(v0)-phis(v1));
                pts[i]=pts[v0+12]+t*(pts[v1+12]-pts[v0+12]);}}

        const MARCHING_CUBES_3D_CASE& cs=table(c);
        printf("case %i\n", c);

        for(int i=0;i<MARCHING_CUBES_3D_CASE::max_elements && cs.elements[i];i++){
            TV_INT face;
            for(int j=0;j<3;j++){
                int e=(cs.elements[i]>>5*j)&31;
                FACE_INDEX<TV::m> fi(e/4,it.index+bits(vertex_lookup[e][0]));
                if(!ht.Get(fi,face(j))){
                    int index=surface.particles.array_collection->Add_Element();
                    face(j)=index;
                    ht.Set(fi,index);
                    surface.particles.X(index)=pts[e];}}
            if(!cs.enclose_inside) exchange(face.x,face.y);
            surface.mesh.elements.Append(face);
        }
    }
    surface.Update_Number_Nodes();
}
template class MARCHING_CUBES_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHING_CUBES_3D<double>;
#endif
