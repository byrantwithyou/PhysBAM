//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EVOLUTION_CALLBACKS
//#####################################################################
#ifndef __SOLIDS_EVOLUTION_CALLBACKS__
#define __SOLIDS_EVOLUTION_CALLBACKS__

#include <Core/Arrays/ARRAY.h>
namespace PhysBAM{

class RIGID_BODY_COLLISION_MANAGER;

template<class TV>
class SOLIDS_EVOLUTION_CALLBACKS
{
    typedef typename TV::SCALAR T;
public:
    SOLIDS_EVOLUTION_CALLBACKS()
    {}

    virtual ~SOLIDS_EVOLUTION_CALLBACKS();

//#####################################################################
    virtual void Update_Solids_Parameters(const T time);
    virtual void Self_Collisions_Begin_Callback(const T time,const int substep);
    virtual void Preprocess_Solids_Substep(const T time,const int substep);
    virtual void Postprocess_Solids_Substep(const T time,const int substep);
    virtual void Align_Deformable_Bodies_With_Rigid_Bodies();
    virtual void Apply_Constraints(const T dt,const T time);
    virtual T Constraints_CFL();
    virtual void Limit_Solids_Dt(T& dt,const T time);
    virtual void Pre_Advance_Cluster_Fracture(const T& dt, const T& time);
    virtual void Post_Advance_Cluster_Fracture(const T& dt, const T& time);
    virtual void Filter_Velocities(const T dt,const T time,const bool velocity_update);
    virtual bool Get_Solid_Source_Velocities(ARRAY<int>& deformable_simplices,ARRAY<T>& deformable_simplex_forces,ARRAY<PAIR<int,int> >& rigid_simplices,
        ARRAY<T>& rigid_simplex_forces,TV& orientation,const T time);
//#####################################################################
};
}
#endif
