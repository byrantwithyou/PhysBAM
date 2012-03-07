//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_3D.h>
using namespace PhysBAM;
static int edge_lookup[8][8];
static int vertex_lookup[12][2];
static int permute_map[6][20];
#define TRI(a,b,c,s) (((s)<<15) | ((c)<<10) | ((b)<<5) | (a))
#define TRI_ORIENT_MAP(x,f) (x?((x&0x8000) | (f[(x>>10)&31]<<10) | (f[x&31]<<5) | f[(x>>5)&31]):0)
//#####################################################################
// Function Case_Table
//#####################################################################
template<class T> const ARRAY<typename MARCHING_CUBES_3D<T>::CASE>& MARCHING_CUBES_3D<T>::
Case_Table()
{
    static ARRAY<CASE> table;
    static bool filled=false;
    if(!filled) Initialize_Case_Table(table);
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
Initialize_Neighbor_Cases(ARRAY<typename MARCHING_CUBES_3D<T>::CASE>& table, int c)
{
    for(int p=0;p<6;p++){
        int b=permute_case(c, permute_map[p]);
        if(!table(b).elements[0]){
            printf("permute[%i] %i -> %i\n", p, c, b);
            for(int i=0;i<max_elements;i++) table(b).elements[i]=TRI_ORIENT_MAP(table(c).elements[i],permute_map[p]);
            for(int i=0;i<sheet_elements;i++) table(b).boundary[i]=TRI_ORIENT_MAP(table(c).boundary[i],permute_map[p]);
            if(p>=3 && table(c).proj_dir!=p-3) table(b).proj_dir=6-table(c).proj_dir-p;
            else table(b).proj_dir=table(c).proj_dir;
            table(b).enclose_inside=table(c).enclose_inside;
            Initialize_Neighbor_Cases(table, b);}}

    if(!has_ambig(c)){
        int b=255-c;
        if(!table(b).elements[0]){
            printf("inv %i -> %i\n", c, b);
            table(b)=table(c);
            table(b).enclose_inside=1-table(b).enclose_inside;
            Initialize_Neighbor_Cases(table, b);}}
}
//#####################################################################
// Function Initialize_Case_Table
//#####################################################################
template<class T> void MARCHING_CUBES_3D<T>::
Initialize_Case_Table(ARRAY<typename MARCHING_CUBES_3D<T>::CASE>& table)
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
    static CASE c0 = {{}, {}, 0, 1};
    static CASE c1 = {{TRI(0,8,4,1)}, {TRI(12,4,8,1)}, 0, 1};
    static CASE c3 = {{TRI(4,9,8,1), TRI(9,4,5,0)}, {TRI(12,4,8,1), TRI(5,13,9,0)}, 0, 1};
    static CASE c6 = {{TRI(5,9,0,1), TRI(10,1,4,1)}, {TRI(5,13,9,1), TRI(4,14,10,1)}, 0, 1};
    static CASE c22 = {{TRI(5,9,0,1), TRI(2,6,8,1), TRI(10,1,4,1)}, {TRI(5,13,9,1), TRI(6,16,8,1), TRI(10,4,14,1)}, 0, 1};
    static CASE c23 = {{TRI(10,1,6,1), TRI(1,2,6,0), TRI(1,5,2,0), TRI(5,9,2,0)}, {TRI(10,6,14,1), TRI(6,16,14,0), TRI(16,12,14,0), TRI(5,13,9,0)}, 0, 1};
    static CASE c24 = {{TRI(2,6,8,1), TRI(5,1,11,1)}, {TRI(6,16,8,1), TRI(5,11,15,1)}, 0, 1};
    static CASE c25 = {{TRI(2,6,4,1), TRI(2,4,0,0), TRI(5,1,11,1)}, {TRI(6,16,4,1), TRI(16,12,4,0), TRI(15,5,11,1)}, 0, 1};
    static CASE c27 = {{TRI(2,6,4,1), TRI(2,4,1,0), TRI(2,1,9,0), TRI(1,11,9,0)}, {TRI(6,12,4,1), TRI(6,16,12,0), TRI(11,13,9,0), TRI(11,15,13,0)}, 0, 1};
    static CASE c60 = {{TRI(7,6,8,1), TRI(7,8,9,0), TRI(4,10,5,1), TRI(5,10,11,0)}, {TRI(6,16,8,1), TRI(17,7,9,0), TRI(4,14,10,1), TRI(15,5,11,0)}, 0, 1};
    static CASE c61 = {{TRI(7,6,9,1), TRI(6,0,9,0), TRI(6,5,0,0), TRI(6,11,5,0), TRI(6,10,11,0)}, {TRI(10,6,14,1), TRI(14,6,16,0), TRI(12,14,16,0), TRI(11,15,5,0), TRI(7,9,17,0)}, 0, 1};
    static CASE c69 = {{TRI(1,6,3,1), TRI(1,8,6,0), TRI(1,0,8,0)}, {TRI(14,18,6,1), TRI(14,6,8,0), TRI(14,8,12,0)}, 0, 1};
    static CASE c85 = {{TRI(1,0,2,1), TRI(1,2,3,0)}, {TRI(12,14,16,1), TRI(14,18,16,0)}, 0, 1};
    static CASE c101 = {{TRI(1,6,3,1), TRI(1,8,6,0), TRI(1,0,8,0), TRI(7,2,9,1)}, {TRI(14,18,6,1), TRI(14,6,8,0), TRI(14,8,12,0), TRI(17,7,9,1)}, 0, 1};
    static CASE c105 = {{TRI(8,4,0,1), TRI(9,7,2,1), TRI(1,11,5,1), TRI(3,10,6,1)}, {TRI(8,12,4,1), TRI(9,17,7,1), TRI(5,11,15,1), TRI(6,10,18,1)}, 0, 1};
    static CASE c125 = {{TRI(11,5,3,1), TRI(3,5,0,0), TRI(3,0,9,0), TRI(3,9,7,0)}, {TRI(17,7,9,1), TRI(11,15,5,0), TRI(12,14,16,0), TRI(14,18,16,0)}, 0, 1};
    static CASE c151 = {{TRI(10,1,6,1), TRI(1,2,6,0), TRI(1,5,2,0), TRI(5,9,2,0), TRI(11,3,7,1)}, {TRI(10,6,14,1), TRI(6,16,14,0), TRI(16,12,14,0), TRI(11,7,19,1)}, 0, 1};
    static CASE c255 = {{}, {}, 0, 0};

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

    for(int c=1;c<256;c++)
        if(!table(c).elements[0])
            printf("missing: %i\n", c);
}
template class MARCHING_CUBES_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHING_CUBES_3D<double>;
#endif
