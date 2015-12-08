//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MPM_PLASTICITY_MODEL__
#define __MPM_PLASTICITY_MODEL__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV>
class MPM_PLASTICITY_MODEL:public NONCOPYABLE
{
public:
    typedef typename TV::SCALAR T;

    virtual ~MPM_PLASTICITY_MODEL() {}

    virtual void Set_Lame_Constants_And_F_Elastic(T mu,T lambda,const DIAGONAL_MATRIX<T,TV::m>& Fe)=0;
    virtual T Yield_Function() const=0;
    virtual T Yield_Function_Final() const=0;
    virtual bool Project_Stress(int max_iterations,T tolerance)=0;
    virtual TV Get_Updated_Sigma() const=0;
    virtual void Set_Plastic_Deformation_Lambda(T plastic_def){PHYSBAM_NOT_IMPLEMENTED();}
    virtual T Get_Updated_Plastic_Deformation_Lambda() const {PHYSBAM_NOT_IMPLEMENTED();}
};
}
#endif
