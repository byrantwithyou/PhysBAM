
//
//  Created by Yuting Wang on 5/23/12.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.
//

#ifndef _levelset_mesh_cutting_h
#define _levelset_mesh_cutting_h

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>

using namespace PhysBAM;

class LEVELSET_MESH_CUTTING_3D
{
public:
    typedef float T;
    typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<int,4> TV_INT4;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<T,4> TV4;


    struct TET
    {
        TV_INT4 parent,indices;
    };

    static void Subdivide(const ARRAY<TV_INT4>& mesh,ARRAY<T>& phi0,ARRAY<T>& phi1,ARRAY<TET>& cut_mesh);
};

#endif
