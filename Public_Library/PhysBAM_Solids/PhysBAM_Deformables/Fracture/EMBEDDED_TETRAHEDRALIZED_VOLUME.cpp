//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TETRAHEDRALIZED_VOLUME
//#####################################################################
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> EMBEDDED_TETRAHEDRALIZED_VOLUME<T>::
EMBEDDED_TETRAHEDRALIZED_VOLUME(TETRAHEDRALIZED_VOLUME<T>& simplicial_object_input)
    :EMBEDDED_OBJECT<TV,3>(simplicial_object_input)
{
}
//#####################################################################
template class EMBEDDED_TETRAHEDRALIZED_VOLUME<float>;
template class EMBEDDED_TETRAHEDRALIZED_VOLUME<double>;
}
