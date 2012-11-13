//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/ELLIPSOID_SIGNED_DISTANCE.h>

namespace PhysBAM{
namespace SIGNED_DISTANCE{
//#####################################################################
// Function Calculate
//#####################################################################
template<class TV,class T_GRID,class T_ARRAY> void Calculate_Approximate(const ELLIPSOID<TV>& ellipsoid,const T_GRID& grid,T_ARRAY& phi,bool verbose)
{
    for(int i=0;i<grid.counts.x;i++)for(int j=0;j<grid.counts.y;j++)for(int ij=0;ij<grid.counts.z;ij++)
        phi(i,j,ij)=ellipsoid.Approximate_Signed_Distance(grid.X(i,j,ij));
}
//#####################################################################

}
}
