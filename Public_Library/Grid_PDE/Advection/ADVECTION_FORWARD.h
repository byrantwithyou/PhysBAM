//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Avi Robinson-Mosher, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Header ADVECTION_FORWARD
//#####################################################################
#ifndef __ADVECTION_FORWARD__
#define __ADVECTION_FORWARD__

namespace PhysBAM{

template<class TV> class FACE_LOOKUP_UNIFORM;

template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> > class ADVECTION;

}
#endif
