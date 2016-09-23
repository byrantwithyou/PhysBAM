//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUT_CELL
//#####################################################################
#ifndef __CUT_CELL__
#define __CUT_CELL__
#include <Core/Arrays/ARRAY.h>
#include <Core/Vectors/VECTOR_1D.h>
#include <Core/Vectors/VECTOR_2D.h>
#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Basic_Geometry/POLYGON.h>
namespace PhysBAM {
template<class T,int d>
class CUT_CELL
{
  public:
    int dominant_element;
    ARRAY<POLYGON<VECTOR<T,d> > > geometry;
    ARRAY<ARRAY<VECTOR<int,d> > > visibility;
};
}
#endif
