//#####################################################################
// Copyright 2019, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Vectors/TWIST.h>
#include "FLUID_BOUNDARY_VECTOR_PB.h"
#include "SOLID_BOUNDARY_VECTOR_PB.h"
#include "SOLID_FLUID_INTERFACE_PB.h"
namespace PhysBAM{

//#####################################################################
// Function Compute_BC
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_BC(const FLUID_SOLVER<TV>* fluid_solver,SOLID_BC<TV>* solid_bc,T time,T dt) const
{
    LOG::printf("DO THIS: %s\n",__PRETTY_FUNCTION__);
}

//#####################################################################
// Function Compute_BC
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_BC(const SOLID_SOLVER<TV>* solid_solver,FLUID_BC<TV>* fluid_bc,T time,T dt) const
{
    LOG::printf("DO THIS: %s\n",__PRETTY_FUNCTION__);
}

//#####################################################################
// Function Interpolate_Velocity
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Interpolate_Velocity(FLUID_BOUNDARY_VECTOR<TV>* u, const SOLID_BOUNDARY_VECTOR<TV>* v)
{
    auto* v_pb=dynamic_cast<const SOLID_BOUNDARY_VECTOR_PB<TV>*>(v);
    auto* u_pb=dynamic_cast<FLUID_BOUNDARY_VECTOR_PB<TV>*>(u);
    u_pb->V.Remove_All();
    for(const auto& e:deformable_weights) u_pb->V.Get_Or_Insert(e.face)+=e.w*v_pb->V.Get(e.p)(e.face.axis);
    for(const auto& e:rigid_weights)
    {
        TWIST<TV> t=v_pb->twist.Get(e.p);
        u_pb->V.Get_Or_Insert(e.face)+=e.w*t.linear(e.face.axis)+e.aw.Dot(t.angular);
    }
}

//#####################################################################
// Function Distribute_Force
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Distribute_Force(SOLID_BOUNDARY_VECTOR<TV>* v, const FLUID_BOUNDARY_VECTOR<TV>* u)
{
    auto* v_pb=dynamic_cast<SOLID_BOUNDARY_VECTOR_PB<TV>*>(v);
    auto* u_pb=dynamic_cast<const FLUID_BOUNDARY_VECTOR_PB<TV>*>(u);
    v_pb->V.Remove_All();
    v_pb->twist.Remove_All();
    for(const auto& e:deformable_weights) v_pb->V.Get_Or_Insert(e.p)(e.face.axis)+=e.w*u_pb->V.Get(e.face);
    for(const auto& e:rigid_weights)
    {
        TWIST<TV>& t=v_pb->twist.Get_Or_Insert(e.p);
        T a=u_pb->V.Get(e.face);
        t.linear(e.face.axis)+=e.w*a;
        t.angular+=e.aw*a;
    }
}

//#####################################################################
// Function Get_Boundary
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Get_Boundary(SOLID_BOUNDARY_VECTOR<TV>* v)
{
    auto* v_pb=dynamic_cast<SOLID_BOUNDARY_VECTOR_PB<TV>*>(v);
    v_pb->V.Remove_All();
    v_pb->twist.Remove_All();
    for(const auto& e:deformable_weights) v_pb->V.Set(e.p,TV());
    for(const auto& e:rigid_weights) v_pb->twist.Set(e.p,TWIST<TV>());
}

//#####################################################################
// Function Compute_Coupling_Weights
//#####################################################################
template<class TV> void SOLID_FLUID_INTERFACE_PB<TV>::
Compute_Coupling_Weights(const SOLID_SOLVER<TV>* solid_solver,const FLUID_SOLVER<TV>* fluid_solver)
{
    


    
#if 0
    struct DEFORMABLE_ENTRY
    {
        int p;
        FACE_INDEX<TV::m> face;
        TV w;
    };
    ARRAY<DEFORMABLE_ENTRY> deformable_weights;
    
    struct RIGID_ENTRY
    {
        int p;
        FACE_INDEX<TV::m> face;
        TWIST<TV> w;
    };
    ARRAY<RIGID_ENTRY> rigid_weights;
#endif

    LOG::printf("DO THIS: %s\n",__PRETTY_FUNCTION__);
}

template class SOLID_FLUID_INTERFACE_PB<VECTOR<float,2> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<float,3> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<double,2> >;
template class SOLID_FLUID_INTERFACE_PB<VECTOR<double,3> >;
}
