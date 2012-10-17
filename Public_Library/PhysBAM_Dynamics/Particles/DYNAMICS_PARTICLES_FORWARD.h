//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DYNAMICS_PARTICLES_FORWARD__
#define __DYNAMICS_PARTICLES_FORWARD__

#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class SPH_PARTICLES;
template<class TV> class PARTICLE_LEVELSET_PARTICLES;
template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLES;
template<class TV> class VORTICITY_PARTICLES;

const ATTRIBUTE_ID ATTRIBUTE_ID_E(10);
const ATTRIBUTE_ID ATTRIBUTE_ID_RHO(11);
const ATTRIBUTE_ID ATTRIBUTE_ID_AGE(12);
const ATTRIBUTE_ID ATTRIBUTE_ID_PHI(13);
const ATTRIBUTE_ID ATTRIBUTE_ID_GRAD_PHI(14);
const ATTRIBUTE_ID ATTRIBUTE_ID_MATERIAL_VOLUME(17);
const ATTRIBUTE_ID ATTRIBUTE_ID_QUANTIZED_COLLISION_DISTANCE(18);

}
#endif
