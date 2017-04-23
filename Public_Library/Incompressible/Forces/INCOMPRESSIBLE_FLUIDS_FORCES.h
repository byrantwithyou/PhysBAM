//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FLUIDS_FORCES
//#####################################################################
#ifndef __INCOMPRESSIBLE_FLUIDS_FORCES__
#define __INCOMPRESSIBLE_FLUIDS_FORCES__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class INCOMPRESSIBLE_FLUIDS_FORCES
{
    typedef typename TV::SCALAR T;
public:

    INCOMPRESSIBLE_FLUIDS_FORCES()
    {}

    INCOMPRESSIBLE_FLUIDS_FORCES(const INCOMPRESSIBLE_FLUIDS_FORCES&) = delete;
    void operator=(const INCOMPRESSIBLE_FLUIDS_FORCES&) = delete;

    virtual ~INCOMPRESSIBLE_FLUIDS_FORCES()
    {}

//#####################################################################
    virtual void Add_Explicit_Forces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Add_Implicit_Forces_Before_Projection(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Add_Implicit_Forces_Projection(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Initialize_Grids(const GRID<TV>& grid)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    
    virtual T CFL(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
    {return 0;}
//#####################################################################
};
}
#endif
