//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TOPOLOGY_BASED_GEOMETRY_FORWARD 
//#####################################################################
#ifndef __TOPOLOGY_BASED_GEOMETRY_FORWARD__
#define __TOPOLOGY_BASED_GEOMETRY_FORWARD__

namespace PhysBAM{

template<class T> class POINT_SIMPLICES_1D;
template<class T> class SEGMENTED_CURVE_2D;
template<class TV> class SEGMENTED_CURVE;
template<class T> class TRIANGULATED_AREA;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class HEXAHEDRALIZED_VOLUME;
template<class TV> class IMPLICIT_OBJECT;
template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class TV,class TRANSFORM> class IMPLICIT_OBJECT_TRANSFORMED;

template<class TV> class ANALYTIC_IMPLICIT_OBJECT;
template<class TV> class IMPLICIT_OBJECT_COMBINED;
template<class TV> class IMPLICIT_OBJECT_COMBINED_EULERIAN;
template<class TV> class MULTIBODY_LEVELSET_IMPLICIT_OBJECT;

}
#endif
