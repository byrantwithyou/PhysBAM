//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SOLID_FLUID_INTERFACE__
#define __SOLID_FLUID_INTERFACE__
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV> class FLUID_SOLVER;
template<class TV> class SOLID_SOLVER;
template<class TV> class FLUID_BC;
template<class TV> class SOLID_BC;
template<class TV> class FLUID_BOUNDARY_VECTOR;
template<class TV> class SOLID_BOUNDARY_VECTOR;

template<class TV>
class SOLID_FLUID_INTERFACE
{
public:
    typedef typename TV::SCALAR T;

    SOLID_FLUID_INTERFACE()=default;
    virtual ~SOLID_FLUID_INTERFACE()=default;

    virtual void Compute_BC(FLUID_SOLVER<TV>* fluid_solver,SOLID_BC<TV>* solid_bc,T time,T dt) const=0;
    virtual void Compute_BC(SOLID_SOLVER<TV>* solid_solver,FLUID_BC<TV>* fluid_bc,T time,T dt) const=0;

    virtual void Interpolate_Velocity(FLUID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v)=0; // Times W
    virtual void Distribute_Force(SOLID_BOUNDARY_VECTOR<TV>* v, const FLUID_BOUNDARY_VECTOR<TV>* u)=0; // Times W^T
    virtual void Get_Boundary(SOLID_BOUNDARY_VECTOR<TV>* v)=0;
};
}
#endif
