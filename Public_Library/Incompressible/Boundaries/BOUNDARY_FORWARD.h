//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header BOUNDARY_FORWARD
//#####################################################################
#ifndef __BOUNDARY_FORWARD__
#define __BOUNDARY_FORWARD__

namespace PhysBAM{

template<class TV,class T2> class BOUNDARY_LINEAR_EXTRAPOLATION;
template<class TV,class T2> class BOUNDARY_MAC_GRID_PERIODIC;
template<class TV,class T2> class BOUNDARY_PERIODIC;
template<class TV,class T2> class BOUNDARY_REFLECTION_UNIFORM;
template<class TV,class T2> class BOUNDARY_REFLECTION_WATER;
template<class TV,class T2> class BOUNDARY_CONSTANT_EXTRAPOLATION_COLLIDABLE;
template<class TV,class T2> class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_PERIODIC;
template<class TV> class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP;
template<class TV> class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP;
template<class TV> class BOUNDARY_PHI_WATER;
template<class TV> class BOUNDARY_SOLID_WALL_SLIP_OUTFLOW;
template<class T> class BOUNDARY_EULER_EQUATIONS_CYLINDRICAL;
template<class T> class BOUNDARY_EULER_EQUATIONS_SPHERICAL;

}
#endif
