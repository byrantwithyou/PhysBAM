//#####################################################################
// Copyright 2017, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LEVELSET_IMPLICIT_UTILITIES__
#define __LEVELSET_IMPLICIT_UTILITIES__
#include <Core/Vectors/VECTOR.h>
#include <string>
namespace PhysBAM{
template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class SEGMENTED_CURVE_2D;

template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >*
Initialize_Implicit_Surface(TRIANGULATED_SURFACE<T>& surface,int max_res);
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,2> >*
Initialize_Implicit_Surface(SEGMENTED_CURVE_2D<T>& surface,int max_res);
template<class T> LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* 
Levelset_From_Tri_File(const std::string& filename,int max_resolution=200);
}
#endif
