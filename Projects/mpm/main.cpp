//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include "RECONSTRUCTION_PARTICLES.h"
using namespace PhysBAM;

int main(int argc,char *argv[])
{
    typedef double T;
    typedef VECTOR<T,2> TV;

    RECONSTRUCTION_PARTICLES<TV> particles;
    
    RANGE<TV> box(TV(0,0),TV(1,1));
    particles.Initialize_Particles(box,0.3);


    return 0;
}
