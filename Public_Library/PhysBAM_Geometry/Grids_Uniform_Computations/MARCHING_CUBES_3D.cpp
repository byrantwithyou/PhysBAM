//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/MARCHING_CUBES_3D.h>
using namespace PhysBAM;
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
int flip(int c,int m)
{
    int r=0;
    for(int i=0;i<8;i++)
        if(c&(1<<i))
            r|=1<<(i^m);
    return r;
}
//#####################################################################
// Function swapyz
//#####################################################################
int swapyz(int c)
{
    int r=0;
    for(int i=0;i<8;i++)
        if(c&(1<<i))
            r|=1<<((i&1)|(i&2)*2|(i&4)/2);
    return r;
}
//#####################################################################
// Function has_ambig
//#####################################################################
int has_ambig(int n)
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
#define TRI(a,b,c,s) (((s)<<15) | ((c)<<10) | ((b)<<5) | (a))
template<class T> void MARCHING_CUBES_3D<T>::
Initialize_Case_Table(ARRAY<typename MARCHING_CUBES_3D<T>::CASE>& table)
{
    table.Resize(256);
    printf("CS %x\n", TRI(8,9,16,0));
    static CASE c0 = {{0}, {0}};
    static CASE c1 = {{TRI(8,16,12,0)}, {TRI(8,4,16,0), TRI(8,1,4,0), TRI(1,5,4,0)}};
    static CASE c3 = {{TRI(12,17,16,0), TRI(17,12,13,0)}, {TRI(16,17,5,0), TRI(5,4,16,0)}};
    static CASE c5 = {{TRI(9,8,16,0), TRI(9,16,18,0)}, {TRI(2,9,18,0), TRI(8,4,16,0), TRI(8,1,4,0), TRI(1,5,4,0)}};
    static CASE c6 = {{TRI(13,17,8,0), TRI(18,9,12,1)}, {TRI(17,5,8,0), TRI(5,4,8,0), TRI(8,4,0,0), TRI(18,6,9,1), TRI(9,6,3,1), TRI(3,6,7,1)}};
    static CASE c7 = {{TRI(16,18,9,0), TRI(17,16,13,0), TRI(16,9,13,0)}, {TRI(2,9,18,0), TRI(16,17,5,0), TRI(5,4,16,0)}};
    static CASE c15 = {{TRI(16,18,19,0), TRI(19,18,17,0)}, {TRI(2,3,18,0), TRI(3,19,18,0), TRI(16,17,5,0), TRI(5,4,16,0)}};
    static CASE c20 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c21 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c22 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c23 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c24 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c25 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c27 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c28 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c29 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c30 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c60 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c61 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c85 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c86 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c90 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c91 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c105 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c107 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c111 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};
    static CASE c125 = {{TRI(1,2,3,0)}, {TRI(4,5,6,0)}};

    table(0)=c0;
    table(1)=c1;
    table(3)=c3;
    table(5)=c5;
    table(6)=c6;
    table(7)=c7;
    table(15)=c15;
    table(20)=c20;
    table(21)=c21;
    table(22)=c22;
    table(23)=c23;
    table(24)=c24;
    table(25)=c25;
    table(27)=c27;
    table(28)=c28;
    table(29)=c29;
    table(30)=c30;
    table(60)=c60;
    table(61)=c61;
    table(85)=c85;
    table(86)=c86;
    table(90)=c90;
    table(91)=c91;
    table(105)=c105;
    table(107)=c107;
    table(111)=c111;
    table(125)=c125;

    UNION_FIND<> uf(256);
    int base[256];
    for(int c=0;c<256;c++){
        base[c]=c;
        uf.Union(c,flip(c,1));
        uf.Union(c,flip(c,2));
        uf.Union(c,flip(c,4));
        uf.Union(c,swapyz(c));
        if(!has_ambig(c)) uf.Union(c,255-c);}

    for(int c=0;c<256;c++){
        int p=uf.Find(c);
        if(base[p]>c) base[p]=c;
        base[c]=base[p];}
}
template class MARCHING_CUBES_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MARCHING_CUBES_3D<double>;
#endif
