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
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
namespace PhysBAM{
template<class T>
class DELAUNAY_TRIANGULATION_2D
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,3> E;
public:
    DELAUNAY_TRIANGULATION_2D(){}
    ~DELAUNAY_TRIANGULATION_2D(){}

    static void Triangulate(const ARRAY_VIEW<TV>& X,TRIANGULATED_AREA<T>& ta,const T length_max=(T)99999,const T theta_min=-(T)1.0);
    static void Test(const TRIANGULATED_AREA<T>& ta);
//#####################################################################
};
}
#endif
