//#####################################################################
// Copyright 2008, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> KRYLOV_SOLVER<T>::
KRYLOV_SOLVER()
    :print_diagnostics(true),print_residuals(false),nullspace_tolerance((T)1e-5),iterations_used(0),residual_magnitude_squared(0),nullspace_measure(0),restart_iterations(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> KRYLOV_SOLVER<T>::
~KRYLOV_SOLVER()
{
}
//#####################################################################
// Function Ensure_Size
//#####################################################################
template<class T> void KRYLOV_SOLVER<T>::
Ensure_Size(ARRAY<KRYLOV_VECTOR_BASE<T>*>& av,const KRYLOV_VECTOR_BASE<T>& v,int size)
{
    for(int i=0,m=min(size,av.m);i<m;i++)
        av(i)->Resize(v);

    if(size<=av.m) return;
    int old=av.m;
    av.Resize(size);
    for(int i=old;i<av.m;i++)
        av(i)=v.Clone_Default();
}
namespace PhysBAM{
template class KRYLOV_SOLVER<float>;
template class KRYLOV_SOLVER<double>;
}
