//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> KRYLOV_VECTOR_BASE<T>::
KRYLOV_VECTOR_BASE()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> KRYLOV_VECTOR_BASE<T>::
~KRYLOV_VECTOR_BASE()
{
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class T> const T& KRYLOV_VECTOR_BASE<T>::
Raw_Get(int i) const
{
    return const_cast<KRYLOV_VECTOR_BASE<T>*>(this)->Raw_Get(i);
}
namespace PhysBAM{
//#####################################################################
// Operator <<
//#####################################################################
template<class T> std::ostream& operator<<(std::ostream& output,const KRYLOV_VECTOR_BASE<T>& x)
{
    output<<"(";
    for(int i=0,n=x.Raw_Size();i<n;i++){
        output<<x.Raw_Get(i);
        if(i<n-1) output<<' ';}
    output<<")";
    return output;
}
//#####################################################################
// Function Random_Fill_Uniform
//#####################################################################
template<class T> void
Random_Fill_Uniform(RANDOM_NUMBERS<T>& rand,KRYLOV_VECTOR_BASE<T>& v,T a,T b)
{
    ARRAY<T> u(v.Raw_Size());
    rand.Fill_Uniform(u,a,b);
    v.Set(u);
}
template class KRYLOV_VECTOR_BASE<float>;
template class KRYLOV_VECTOR_BASE<double>;
template std::ostream& operator<< <float>(std::ostream&,const KRYLOV_VECTOR_BASE<float>&);
template std::ostream& operator<< <double>(std::ostream&,const KRYLOV_VECTOR_BASE<double>&);
}
