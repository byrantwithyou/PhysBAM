//#####################################################################
// Copyright 2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPH_CALLBACKS
//##################################################################### 
#ifndef __SPH_CALLBACKS__
#define __SPH_CALLBACKS__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID;

template<class TV>
class SPH_CALLBACKS
{    
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:

    SPH_CALLBACKS()
    {}

    virtual ~SPH_CALLBACKS()
    {}

//#####################################################################
    virtual void Adjust_SPH_Particle_For_Domain_Boundaries(SPH_PARTICLES<TV>& particles,const int index,TV& V,const T dt,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual bool Adjust_SPH_Particle_For_Objects(SPH_PARTICLES<TV>& particles,const int index,TV& V,const T dt,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return true;} // return false if particle should be deleted
    virtual void Do_Something_With_Density(const GRID<TV> &grid,const ARRAY<T,TV_INT> &cell_weight)const{}
    virtual T Target_Density_Factor(const TV& location,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 1;}
//#####################################################################
};
}
#endif
