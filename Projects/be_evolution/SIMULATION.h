//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMULATION
//#####################################################################
#ifndef __SIMULATION__
#define __SIMULATION__
#include <Tools/Nonlinear_Equations/NEWTONS_METHOD.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
using namespace PhysBAM;
template<class TV>
class SIMULATION
{
    typedef typename TV::SCALAR T;
public:
    SOLID_BODY_COLLECTION<TV> solid_body_collection;
    T time;
    NEWTONS_METHOD<T> nm;

    SIMULATION();
    ~SIMULATION();

    void Advance_One_Time_Step_Position(const T dt);
};
#endif
