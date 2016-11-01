//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHOTON__
#define __PHOTON__

#include <Core/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class PHOTON
{
    typedef VECTOR<T,3> TV;
public:
    typedef enum {KDTREE_LEAF=0,KDTREE_SPLIT_X=1,KDTREE_SPLIT_Y=2,KDTREE_SPLIT_Z=3} PHOTON_KDTREE_SPLIT_AXIS;

    TV location;
    TV direction;
    TV power;
    PHOTON_KDTREE_SPLIT_AXIS kdtree_split_axis;

//#####################################################################
};
}
#endif
