//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __REINITIALIZATION__
#define __REINITIALIZATION__

#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{
template<class T,class TV> void Reinitialize(LEVELSET<TV>& levelset,int time_steps,T time,T half_band_width,T extra_band,T cfl,int temporal_order,int spatial_order,int process_sign);
}
#endif
