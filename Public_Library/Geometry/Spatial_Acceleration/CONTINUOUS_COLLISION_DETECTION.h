//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Joseph Teran, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CONTINUOUS_COLLISION_DETECTION__
#define __CONTINUOUS_COLLISION_DETECTION__
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class STRUCTURE;
template<int d>
struct CCD_PAIR
{
    int s0,s1;
    VECTOR<int,d+1> f; // pf: 0 1 1 1, ee: 0 0 1 1
};
template<class TV>
void Continuous_Collision_Detection(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
    ARRAY_VIEW<STRUCTURE<TV>*> structures,ARRAY<CCD_PAIR<TV::m> >* point_face,
    ARRAY<CCD_PAIR<TV::m> >* edge_edge,typename TV::SCALAR thickness=0);
}
#endif

