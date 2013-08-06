//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include "MINIMIZATION_OBJECTIVE.h"
#include "SIMULATION.h"
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SIMULATION<TV>::
SIMULATION()
{
    nm.max_iterations=100000;
    nm.max_krylov_iterations=2000;
    nm.krylov_tolerance=1;
    nm.fail_on_krylov_not_converged=false;
    nm.tolerance=1e-5;
    nm.angle_tolerance=1e-2;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SIMULATION<TV>::
~SIMULATION()
{
}
//#####################################################################
// Function Advance_One_Time_Step_Position
//#####################################################################
template<class TV> void SIMULATION<TV>::
Advance_One_Time_Step_Position(const T dt)
{
    LOG::SCOPE scope("Advance_One_Time_Step_Position");
    MINIMIZATION_OBJECTIVE<TV> obj(solid_body_collection,dt,time);
    KRYLOV_VECTOR_BASE<T>* x0 = obj.x0.Clone_Default();
    *x0=obj.x0;
    obj.Test(*x0,obj);

    bool converged=nm.Newtons_Method(obj,obj,*x0);
    PHYSBAM_ASSERT(converged);
    delete x0;
}
template class SIMULATION<VECTOR<float,2> >;
template class SIMULATION<VECTOR<float,3> >;
template class SIMULATION<VECTOR<double,2> >;
template class SIMULATION<VECTOR<double,3> >;
