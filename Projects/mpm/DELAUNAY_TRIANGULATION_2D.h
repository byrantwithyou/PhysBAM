//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DELAUNAY_TRIANGULATION_2D
//#####################################################################
#ifndef __DELAUNAY_TRIANGULATION_2D__
#define __DELAUNAY_TRIANGULATION_2D__
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
namespace PhysBAM{
template<class T>
class DELAUNAY_TRIANGULATION_2D
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,3> E;
public:
    DELAUNAY_TRIANGULATION_2D(){}
    ~DELAUNAY_TRIANGULATION_2D(){}

    static void Triangulate(const ARRAY_VIEW<TV>& X,ARRAY<E>& elements,const T threshold=(T)99999);
//#####################################################################
};
}
#endif
