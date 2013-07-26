//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_CONFINEMENT
//#####################################################################
#ifndef __VORTICITY_CONFINEMENT__
#define __VORTICITY_CONFINEMENT__

#include <Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV>
class VORTICITY_CONFINEMENT:public INCOMPRESSIBLE_FLUIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_body_list;
    T_FACE_ARRAYS_BOOL* valid_mask;
    bool use_variable_vorticity_confinement;
    T vorticity_confinement;
    T_ARRAYS_SCALAR variable_vorticity_confinement;

    VORTICITY_CONFINEMENT(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_body_list=0,T_FACE_ARRAYS_BOOL* valid_mask=0,const bool use_variable_vorticity_confinement=false,const T vorticity_confinement=.3);
    virtual ~VORTICITY_CONFINEMENT();

    void Set_Vorticity_Confinement(const T vorticity_confinement_input=.3)
    {vorticity_confinement=vorticity_confinement_input;}

//#####################################################################
    void Apply_Vorticity_Confinement_Force(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities,ARRAY<TV,TV_INT>& F);
    virtual void Compute_Vorticity_Confinement_Force(const GRID<TV>& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T_FACE_ARRAYS_BOOL* valid_mask,ARRAY<TV,TV_INT>& F);
    void Add_Explicit_Forces(const GRID<TV>& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Add_Implicit_Forces_Projection(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Initialize_Grids(const GRID<TV>& grid) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
