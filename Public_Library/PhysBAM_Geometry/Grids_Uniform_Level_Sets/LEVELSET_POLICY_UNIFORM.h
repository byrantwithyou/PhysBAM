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
template<class TV> class LEVELSET_IMPLICIT_OBJECT;

template<class T_GRID> class LEVELSET_UNIFORM;
template<class T_GRID> class LEVELSET_1D;
template<class T_GRID> class LEVELSET_2D;
template<class T_GRID> class LEVELSET_3D;
template<class T_GRID> class FAST_LEVELSET;
template<class T_GRID> class LEVELSET_MULTIPLE;
template<class T_GRID> class PARTICLE_LEVELSET_UNIFORM;
template<class T_GRID> class PARTICLE_LEVELSET_EVOLUTION_UNIFORM;
template<class T_GRID,class T2> class EXTRAPOLATION_UNIFORM;

template<class T> struct LEVELSET_POLICY<GRID<VECTOR<T,1> > >
{
    typedef LEVELSET_1D<T> LEVELSET;
    typedef FAST_LEVELSET<GRID<VECTOR<T,1> > > FAST_LEVELSET_T;
    typedef PhysBAM::LEVELSET_IMPLICIT_OBJECT<VECTOR<T,1> > LEVELSET_IMPLICIT_OBJECT;
    typedef EXTRAPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T> EXTRAPOLATION_SCALAR;
    typedef EXTRAPOLATION_UNIFORM<GRID<VECTOR<T,1> >,VECTOR<T,1> > EXTRAPOLATION_VECTOR;
    typedef PARTICLE_LEVELSET_UNIFORM<GRID<VECTOR<T,1> > > PARTICLE_LEVELSET;
    typedef PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<T,1> > > PARTICLE_LEVELSET_EVOLUTION;
};

template<class T> struct LEVELSET_POLICY<GRID<VECTOR<T,2> > >
{
    typedef LEVELSET_2D<GRID<VECTOR<T,2> > > LEVELSET;
    typedef FAST_LEVELSET<GRID<VECTOR<T,2> > > FAST_LEVELSET_T;
    typedef PhysBAM::LEVELSET_IMPLICIT_OBJECT<VECTOR<T,2> > LEVELSET_IMPLICIT_OBJECT;
    typedef EXTRAPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T> EXTRAPOLATION_SCALAR;
    typedef EXTRAPOLATION_UNIFORM<GRID<VECTOR<T,2> >,VECTOR<T,2> > EXTRAPOLATION_VECTOR;
    typedef PARTICLE_LEVELSET_UNIFORM<GRID<VECTOR<T,2> > > PARTICLE_LEVELSET;
    typedef PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<T,2> > > PARTICLE_LEVELSET_EVOLUTION;
};

template<class T> struct LEVELSET_POLICY<GRID<VECTOR<T,3> > >
{
    typedef LEVELSET_3D<GRID<VECTOR<T,3> > > LEVELSET;
    typedef FAST_LEVELSET<GRID<VECTOR<T,3> > > FAST_LEVELSET_T;
    typedef PhysBAM::LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> > LEVELSET_IMPLICIT_OBJECT;
    typedef EXTRAPOLATION_UNIFORM<GRID<VECTOR<T,3> >,T> EXTRAPOLATION_SCALAR;
    typedef EXTRAPOLATION_UNIFORM<GRID<VECTOR<T,3> >,VECTOR<T,3> > EXTRAPOLATION_VECTOR;
    typedef PARTICLE_LEVELSET_UNIFORM<GRID<VECTOR<T,3> > > PARTICLE_LEVELSET;
    typedef PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<T,3> > > PARTICLE_LEVELSET_EVOLUTION;
};

}
#endif
