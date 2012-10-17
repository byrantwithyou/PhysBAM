//#####################################################################
// Copyright 2006-2009, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_POLICY_UNIFORM 
//#####################################################################
#ifndef __LEVELSET_POLICY_UNIFORM__
#define __LEVELSET_POLICY_UNIFORM__

#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID> class LEVELSET_UNIFORM;
template<class T_GRID> class LEVELSET_1D;
template<class T_GRID> class LEVELSET_2D;
template<class T_GRID> class LEVELSET_3D;

template<class T> struct LEVELSET_POLICY<GRID<VECTOR<T,1> > >
{
    typedef LEVELSET_1D<T> LEVELSET;
};

template<class T> struct LEVELSET_POLICY<GRID<VECTOR<T,2> > >
{
    typedef LEVELSET_2D<GRID<VECTOR<T,2> > > LEVELSET;
};

template<class T> struct LEVELSET_POLICY<GRID<VECTOR<T,3> > >
{
    typedef LEVELSET_3D<GRID<VECTOR<T,3> > > LEVELSET;
};

}
#endif
