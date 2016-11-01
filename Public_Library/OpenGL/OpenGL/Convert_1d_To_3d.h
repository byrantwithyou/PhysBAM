//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Convert_1d_To_3d
//#####################################################################
#ifndef __Convert_1d_To_3d__
#define __Convert_1d_To_3d__

#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> inline VECTOR<T,3> Convert_1d_To_3d(const VECTOR<T,1>& v)
{return VECTOR<T,3>(v.x,0,0);}

template<class T> inline RANGE<VECTOR<T,3> > Convert_1d_To_3d(const RANGE<VECTOR<T,1> >& box)
{return RANGE<VECTOR<T,3> >(VECTOR<T,3>(box.min_corner),VECTOR<T,3>(box.max_corner));}

template<class T> inline ROTATION<VECTOR<T,3> > Convert_1d_To_3d(const ROTATION<VECTOR<T,1> >& c)
{return ROTATION<VECTOR<T,3> >();}

template<class T> inline FRAME<VECTOR<T,3> > Convert_1d_To_3d(const FRAME<VECTOR<T,1> >& f)
{return FRAME<VECTOR<T,3> >(Convert_1d_To_3d(f.t),Convert_1d_To_3d(f.r));}

}
#endif
