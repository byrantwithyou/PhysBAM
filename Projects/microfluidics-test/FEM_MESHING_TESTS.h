//#####################################################################
// Copyright 2018.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FEM_MESHING_TESTS__
#define __FEM_MESHING_TESTS__
#include "PARSE_DATA_FEM.h"

namespace PhysBAM{

template<typename TV>
void Test_Degree2_Joint(JOINT_TYPE jt,typename TV::SCALAR a0,typename TV::SCALAR a1,typename TV::SCALAR da);

template<typename TV>
void Test_Degree2_Circle(JOINT_TYPE jt,typename TV::SCALAR h0,typename TV::SCALAR h1,typename TV::SCALAR dh);

template<typename TV>
void Test_Degree3_Joint(JOINT_TYPE jt,typename TV::SCALAR h,int n,int seed);
}
#endif

