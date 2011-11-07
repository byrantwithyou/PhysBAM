//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_ORIGIN_AREAS
//##################################################################### 
#ifndef __TRIANGLE_ORIGIN_AREAS__
#define __TRIANGLE_ORIGIN_AREAS__

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
namespace TRIANGLE_ORIGIN_AREAS
{

template<class T,int n>
struct VOL_DATA
{
    T V;
    VECTOR<T,3> G[n];
    MATRIX<T,3> H[n][n];
};

template<class T>
struct PT_DATA
{
    int n;
    int index[5];
    VECTOR<T,3> V;
    MATRIX<T,3> G[5];
    MATRIX<T,3> H[3][5][5];
};

template<class T,class TV> void Volume_From_Triangles(VOL_DATA<T,6>& data,TV A,TV B,TV C,TV D,TV E,TV F);
}
}
#endif

