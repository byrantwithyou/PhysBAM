
//
//  Created by Yuting Wang on 5/23/12.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.
//

#ifndef _levelset_mesh_cutting_h
#define _levelset_mesh_cutting_h

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Data_Structures/UNION_FIND.h>

using namespace PhysBAM;

class LEVELSET_MESH_CUTTING_3D
{
public:
    typedef double T;
    typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<int,4> TV_INT4;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,4> TV4;

    struct TET
    {
        TV_INT4 parent,indices;
        VECTOR<TV4,4> weights;
    };

    static void Subdivide(const ARRAY<TV_INT4>& mesh,ARRAY<T>& phi0,ARRAY<T>& phi1,ARRAY<TET>& cut_mesh);
};

#endif
