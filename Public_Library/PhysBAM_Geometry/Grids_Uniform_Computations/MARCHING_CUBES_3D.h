//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MARCHING_CUBES_3D
//#####################################################################
#ifndef __MARCHING_CUBES_3D__
#define __MARCHING_CUBES_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> class TRIANGLE_3D;

/*
  scccccbbbbbaaaaa
  aaaaa, bbbbb, ccccc = vertex index (0-7 cube vertices, 8-19 cube edges)
  s = sheet index (0 or 1)
 */
template<class T>
class MARCHING_CUBES_3D
{
public:
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    enum WORKAROUND {max_elements=5,sheet_elements=5};

    struct CASE
    {
        unsigned short elements[max_elements];
        unsigned short boundary[sheet_elements];
        unsigned short proj_dir;
        unsigned short enclose_inside;
    };

    MARCHING_CUBES_3D() {}
    ~MARCHING_CUBES_3D() {}

    static const ARRAY<CASE>& Case_Table();
    static void Initialize_Case_Table(ARRAY<CASE>& table);
    static void Initialize_Neighbor_Cases(ARRAY<typename MARCHING_CUBES_3D<T>::CASE>& table, int c);
//#####################################################################
};
}
#endif
