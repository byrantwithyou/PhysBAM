//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/PARAMETRIC_LINE.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class T> PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)>::
~PARAMETRIC_LINE()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> void PARAMETRIC_LINE<T,T(KRYLOV_VECTOR_BASE<T>&)>::
Compute(const T t,T* ddg,T* dg,T* g) const PHYSBAM_OVERRIDE
{
    tmp.Copy(t,dx,x);
    PHYSBAM_ASSERT((!ddg || system) && (!dg || tmp2));
    f.Compute(tmp,ddg?system:0,dg?tmp2:0,g);
    if(dg) *dg=tmp2->Dot(dx);
    if(ddg){
        system->Multiply(dx,tmp);
        *ddg=tmp.Dot(dx);}
}
template class PARAMETRIC_LINE<double,double (KRYLOV_VECTOR_BASE<double>&)>;
template class PARAMETRIC_LINE<float,float (KRYLOV_VECTOR_BASE<float>&)>;
}
