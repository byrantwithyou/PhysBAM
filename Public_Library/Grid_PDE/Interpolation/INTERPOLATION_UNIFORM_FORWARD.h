//#####################################################################
// Copyright 2009, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERPOLATION_UNIFORM_FORWARD__
#define __INTERPOLATION_UNIFORM_FORWARD__

namespace PhysBAM{
template<class TV> class FACE_LOOKUP_UNIFORM;
template<class TV> class FACE_LOOKUP_BINARY_UNIFORM;

template<class TV> class LINEAR_INTERPOLATION_MAC_HELPER;
template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class CATMULL_ROM_SPLINE_INTERPOLATION;
template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class CUBIC_MN_INTERPOLATION_UNIFORM;
template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class CUBIC_SPLINE_INTERPOLATION_UNIFORM;

template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class INTERPOLATION_UNIFORM;
template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class LINEAR_INTERPOLATION_UNIFORM;

template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM;

template<class TV,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class AVERAGING_UNIFORM;

}
#endif
