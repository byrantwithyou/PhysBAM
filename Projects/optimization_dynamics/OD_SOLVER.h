//#####################################################################
// Copyright 2015
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __OD_SOLVER__
#define __OD_SOLVER__
#include <Tools/Vectors/VECTOR.h>
#include "OD_EXAMPLE.h"
namespace PhysBAM{

template<class TV>
class OD_SOLVER
{
    typedef typename TV::SCALAR T;
public:

    OD_EXAMPLE<TV>& example;
    int max_iterations;

    OD_SOLVER(OD_EXAMPLE<TV>& example):example(example),max_iterations(100){}
    virtual ~OD_SOLVER(){};
    
    virtual void Solve(T dt){};
//#####################################################################
};
}

#endif
