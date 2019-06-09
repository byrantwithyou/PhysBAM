//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_FEM__
#define __ANALYTIC_FEM__
#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Analytic_Tests/ANALYTIC_VECTOR.h>
#include "BLOCK_MATRIX.h"
#include "BLOCK_VECTOR.h"
#include "MATRIX_CONSTRUCTION_FEM.h"

namespace PhysBAM{

template<class TV>
struct ANALYTIC_FEM
{
    typedef typename TV::SCALAR T;
    ANALYTIC_VECTOR<TV>* analytic_velocity=0;
    ANALYTIC_SCALAR<TV>* analytic_pressure=0;

    MATRIX_CONSTRUCTION_FEM<TV>& mc;

    ANALYTIC_FEM(MATRIX_CONSTRUCTION_FEM<TV>& mc);
    ~ANALYTIC_FEM();
    
    TV Traction(const TV& N,const TV& X) const;
    TV Force(const TV& X) const;
    bool Check_Analytic_Solution(bool dump) const;
    void Compute_RHS();

    void Set_Velocity(const char* s)
    {
        analytic_velocity=new ANALYTIC_VECTOR_PROGRAM<TV>(s);
    }

    void Set_Pressure(const char* s)
    {
        analytic_pressure=new ANALYTIC_SCALAR_PROGRAM<TV>(s);
    }
};

template<class T,int d>
bool Check_Solution(const MATRIX_CONSTRUCTION_FEM<VECTOR<T,d> >& mc,
    std::function<VECTOR<T,d>(const VECTOR<T,d>&)> fv,
    std::function<T(const VECTOR<T,d>&)> fp,bool dump);
}
#endif
