//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_PARTICLES_FORWARD__
#define __MPM_PARTICLES_FORWARD__

#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class MPM_PARTICLES;

const ATTRIBUTE_ID ATTRIBUTE_ID_DENSITY(126);
const ATTRIBUTE_ID ATTRIBUTE_ID_VOLUME(127);
const ATTRIBUTE_ID ATTRIBUTE_ID_FE(128);
const ATTRIBUTE_ID ATTRIBUTE_ID_FP(129);
}
#endif
