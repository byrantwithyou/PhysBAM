//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include "MPM_PARTICLE.h"
#include "MPM_CONSTITUTIVE_MODEL.h"

using namespace PhysBAM;

int main(int argc,char *argv[])
{
    static const int dimension=2;
    static const int influence_n=4;
    typedef double T;
    typedef VECTOR<T,dimension> TV;
    typedef VECTOR<int,dimension> TV_INT;
    
    TV_INT grid_counts(6,7);
    RANGE<TV> grid_box(TV(0,0),TV(5,6));
    GRID<TV> grid(grid_counts,grid_box);
    
    MPM_PARTICLE<TV,influence_n> p(1,1,TV(2.75,3.5),TV(0,0),1,grid);

    LOG::cout<<p.influence_corner<<std::endl;
    LOG::cout<<p.weights<<std::endl;

    return 0;
}
