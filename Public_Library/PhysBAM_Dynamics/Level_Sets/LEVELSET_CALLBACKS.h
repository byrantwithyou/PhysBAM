//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_CALLBACKS
//##################################################################### 
#ifndef __LEVELSET_CALLBACKS__
#define __LEVELSET_CALLBACKS__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;
template<class T_GRID> class LEVELSET_MULTIPLE;
template<class TV> class LEVELSET;

template<class T_GRID>
class LEVELSET_CALLBACKS
{    
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
public:

    LEVELSET_CALLBACKS()
    {}

    virtual ~LEVELSET_CALLBACKS()
    {}

//#####################################################################
    virtual void Get_Levelset_Velocity(const T_GRID& grid,LEVELSET<TV>& levelset,T_FACE_ARRAYS_SCALAR& face_velocity,const T time=0) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Get_Levelset_Velocity(const T_GRID& grid,LEVELSET_MULTIPLE<T_GRID>& levelset,T_FACE_ARRAYS_SCALAR& face_velocity,const T time=0) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
        const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual bool Adjust_Particle_For_Objects(TV& X,TV& V,const T r, const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
        {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return true;} // return false if particle should be deleted
//#####################################################################
};
}
#endif
