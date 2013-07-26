//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FLUIDS_FORCES
//#####################################################################
#ifndef __INCOMPRESSIBLE_FLUIDS_FORCES__
#define __INCOMPRESSIBLE_FLUIDS_FORCES__

#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class INCOMPRESSIBLE_FLUIDS_FORCES:public NONCOPYABLE
{
    
    typedef typename TV::SCALAR T;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
public:

    INCOMPRESSIBLE_FLUIDS_FORCES()
    {}

    virtual ~INCOMPRESSIBLE_FLUIDS_FORCES()
    {}

//#####################################################################
    virtual void Add_Explicit_Forces(const GRID<TV>& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Add_Implicit_Forces_Before_Projection(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Add_Implicit_Forces_Projection(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Initialize_Grids(const GRID<TV>& grid)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    
    virtual T CFL(const GRID<TV>& grid,const T_FACE_ARRAYS_SCALAR& face_velocities)
    {return 0;}
//#####################################################################
};
}
#endif
