//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Convert_2d_To_3d
//#####################################################################
#ifndef __Convert_2d_To_3d__
#define __Convert_2d_To_3d__

#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/ROTATION.h>
#include <Core/Vectors/COMPLEX.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> inline VECTOR<T,3> Convert_2d_To_3d(const VECTOR<T,2>& v)
{return VECTOR<T,3>(v.x,v.y,0);}

template<class T> inline RANGE<VECTOR<T,3> > Convert_2d_To_3d(const RANGE<VECTOR<T,2> >& box)
{return RANGE<VECTOR<T,3> >(box.min_corner.Append(0),box.max_corner.Append(0));}

template<class T> inline ROTATION<VECTOR<T,3> > Convert_2d_To_3d(const ROTATION<VECTOR<T,2> >& c)
{COMPLEX<T> sqrt_c=c.Complex().Sqrt();return ROTATION<VECTOR<T,3> >::From_Components(sqrt_c.re,0,0,sqrt_c.im);}

template<class T> inline FRAME<VECTOR<T,3> > Convert_2d_To_3d(const FRAME<VECTOR<T,2> >& f)
{return FRAME<VECTOR<T,3> >(Convert_2d_To_3d(f.t),Convert_2d_To_3d(f.r));}

}
#endif
