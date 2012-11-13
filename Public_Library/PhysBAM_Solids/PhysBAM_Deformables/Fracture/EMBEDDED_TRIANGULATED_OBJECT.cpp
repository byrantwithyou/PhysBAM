//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_TRIANGULATED_OBJECT
//#####################################################################
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TRIANGULATED_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EMBEDDED_TRIANGULATED_OBJECT<TV>::
EMBEDDED_TRIANGULATED_OBJECT(T_TRIANGULATED_OBJECT& simplicial_object_input)
    :EMBEDDED_OBJECT<TV,2>(simplicial_object_input)
{
}
//#####################################################################
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,2> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,3> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,2> >;
template class EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,3> >;
}
