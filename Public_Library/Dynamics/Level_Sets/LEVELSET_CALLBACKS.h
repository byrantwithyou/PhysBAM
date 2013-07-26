//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_CALLBACKS
//##################################################################### 
#ifndef __LEVELSET_CALLBACKS__
#define __LEVELSET_CALLBACKS__

#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV> class LEVELSET_MULTIPLE;
template<class TV> class LEVELSET;
template<class TV> class PARTICLE_LEVELSET_PARTICLES;

template<class TV>
class LEVELSET_CALLBACKS
{    
    typedef typename TV::SCALAR T;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
public:

    LEVELSET_CALLBACKS()
    {}

    virtual ~LEVELSET_CALLBACKS();

//#####################################################################
    virtual void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET<TV>& levelset,T_FACE_ARRAYS_SCALAR& face_velocity,const T time=0) const;
    virtual void Get_Levelset_Velocity(const GRID<TV>& grid,LEVELSET_MULTIPLE<TV>& levelset,T_FACE_ARRAYS_SCALAR& face_velocity,const T time=0) const;
    virtual void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,
        const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    virtual bool Adjust_Particle_For_Objects(TV& X,TV& V,const T r, const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
//#####################################################################
};
}
#endif
