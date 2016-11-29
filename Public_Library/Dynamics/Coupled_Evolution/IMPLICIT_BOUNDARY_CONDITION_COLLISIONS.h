//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS
//#####################################################################
#ifndef __IMPLICIT_BOUNDARY_CONDITION_COLLISIONS__
#define __IMPLICIT_BOUNDARY_CONDITION_COLLISIONS__
#include <Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class TV>
class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX T_SIMPLEX;

    COLLISION_BODY_COLLECTION<TV>& collision_geometry_collection;
    bool use_implicit_geometry;

public:
    IMPLICIT_BOUNDARY_CONDITION_COLLISIONS(COLLISION_BODY_COLLECTION<TV>& collision_geometry_collection_input,
        const bool use_implicit_geometry_input);

    virtual ~IMPLICIT_BOUNDARY_CONDITION_COLLISIONS();

//#####################################################################
    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::m> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T time) override;
//#####################################################################
};
}
#endif
