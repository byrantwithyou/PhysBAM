//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
using namespace PhysBAM;
template<class T,class SIZE> 
struct STATIC_POLYNOMIAL
{

    STATIC_TENSOR<T,SIZE> terms;

    STATIC_POLYNOMIAL();
    ~STATIC_POLYNOMIAL();

    template<class SIZE2>
    STATIC_POLYNOMIAL<T,typename STATIC_SIZES::SUM<SIZE,SIZE2>::TYPE> operator* (const STATIC_POLYNOMIAL<T,SIZE2>& p)
    {
        STATIC_POLYNOMIAL<T,typename STATIC_SIZES::SUM<SIZE,SIZE2>::TYPE> r;
        RANGE<typename SIZE::T_VECTOR> size(typename SIZE::T_VECTOR(),SIZE::As_Vector());
        RANGE<typename SIZE2::T_VECTOR> size2(typename SIZE2::T_VECTOR(),SIZE2::As_Vector());
        for(UNIFORM_ARRAY_ITERATOR<SIZE::rank> it(size);it.Valid();it.Next())
            for(UNIFORM_ARRAY_ITERATOR<SIZE2::rank> it2(size2);it2.Valid();it2.Next())
                r.terms(it.size+it2.size)+=terms(it.size)*p.terms(it2.size);
        return r;
    }
template class STATIC_POLYNOMIAL<VECTOR<float,1> >;
#ifndef COMPILATE_WITHOUT_DOUBLE_SUPPORT
template class STATIC_POLYNOMIAL<VECTOR<double,1> >;
#endif
