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
const ATTRIBUTE_ID ATTRIBUTE_ID_PARTICLE_DOMAIN(131);
const ATTRIBUTE_ID ATTRIBUTE_ID_MU(132);
const ATTRIBUTE_ID ATTRIBUTE_ID_LAMBDA(133);
const ATTRIBUTE_ID ATTRIBUTE_ID_MU0(134);
const ATTRIBUTE_ID ATTRIBUTE_ID_LAMBDA0(135);
const ATTRIBUTE_ID ATTRIBUTE_ID_USE_PLASTICITY_YIELD(136);
const ATTRIBUTE_ID ATTRIBUTE_ID_USE_PLASTICITY_CLAMP(137);
const ATTRIBUTE_ID ATTRIBUTE_ID_YIELD_MIN(138);
const ATTRIBUTE_ID ATTRIBUTE_ID_YIELD_MAX(139);
const ATTRIBUTE_ID ATTRIBUTE_ID_CLAMP_MIN(140);
const ATTRIBUTE_ID ATTRIBUTE_ID_CLAMP_MAX(141);
const ATTRIBUTE_ID ATTRIBUTE_ID_USE_VISCO_PLASTICITY(142);
const ATTRIBUTE_ID ATTRIBUTE_ID_VISCO_NU(143);
const ATTRIBUTE_ID ATTRIBUTE_ID_VISCO_TAU(144);
const ATTRIBUTE_ID ATTRIBUTE_ID_VISCO_KAPPA(145);
}
#endif
