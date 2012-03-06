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
    static CASE c0 = {{}, {}, 0, 1};
    static CASE c1 = {{TRI(0,8,4,1)}, {TRI(12,4,8,1)}, 0, 1};
    static CASE c3 = {{TRI(4,9,8,1), TRI(9,4,5,0)}, {TRI(12,4,8,1), TRI(5,13,9,0)}, 0, 1};
    static CASE c6 = {{TRI(5,9,0,1), TRI(10,1,4,1)}, {TRI(5,13,9,1), TRI(4,14,10,1)}, 0, 1};
    static CASE c22 = {{TRI(5,9,0,1), TRI(2,6,8,1), TRI(10,1,4,1)}, {TRI(5,13,9,1), TRI(6,16,8,1), TRI(10,4,14,1)}, 0, 1};
    static CASE c23 = {{TRI(10,1,6,1), TRI(1,2,6,0), TRI(1,5,2,0), TRI(5,9,2,0)}, {TRI(10,6,14,1), TRI(6,16,14,0), TRI(16,12,14,0)}, 0, 1};
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
