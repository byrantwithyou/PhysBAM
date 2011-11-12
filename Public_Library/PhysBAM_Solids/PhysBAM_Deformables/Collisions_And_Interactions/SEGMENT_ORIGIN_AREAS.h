//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_ORIGIN_AREAS
//##################################################################### 
#ifndef __SEGMENT_ORIGIN_AREAS__
#define __SEGMENT_ORIGIN_AREAS__

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
namespace ORIGIN_AREAS
{

template<class T,int m,int n>
struct VOL_DATA
{
    T V;
    VECTOR<T,m> G[n];
    MATRIX<T,m> H[n][n];

    VOL_DATA(){V=0;}
};

template<class T,int m,int n> void Clear(VOL_DATA<T,m,n>& data);
template<class T,class TV> void Volume_From_Simplices(VOL_DATA<T,2,4>& data,TV A[4]);
}
}
#endif

