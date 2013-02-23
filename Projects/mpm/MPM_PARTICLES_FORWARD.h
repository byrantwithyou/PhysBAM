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
const ATTRIBUTE_ID ATTRIBUTE_ID_MATERIAL_COORDINATE(127);
const ATTRIBUTE_ID ATTRIBUTE_ID_VOLUME(128);
const ATTRIBUTE_ID ATTRIBUTE_ID_FE(129);
const ATTRIBUTE_ID ATTRIBUTE_ID_FP(130);
}
#endif
