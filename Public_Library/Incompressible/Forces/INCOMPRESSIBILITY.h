//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBILITY
//#####################################################################
#ifndef __INCOMPRESSIBILITY__
#define __INCOMPRESSIBILITY__

#include <Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{
template<class TV> class PROJECTION_UNIFORM;

template<class TV>
class INCOMPRESSIBILITY:public INCOMPRESSIBLE_FLUIDS_FORCES<TV>
{
    
    typedef typename TV::SCALAR T;
    PROJECTION_UNIFORM<TV>& projection;
public:

    INCOMPRESSIBILITY(PROJECTION_UNIFORM<TV>& projection_input);
    virtual ~INCOMPRESSIBILITY();

//#####################################################################
    void Add_Explicit_Forces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time) override {}
    void Add_Implicit_Forces_Projection(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time) override;
    void Initialize_Grids(const GRID<TV>& grid) override;
//#####################################################################
};
}
#endif
