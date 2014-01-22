//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ADAPTIVE_QUOTIENT__
#define __ADAPTIVE_QUOTIENT__
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_BASE.h>
#include <Tools/Adaptive_Arithmetic/ADAPTIVE_FORWARD.h>
#include <Tools/Adaptive_Arithmetic/EXACT_ARITHMETIC_POLICY.h>
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Utilities/STATIC_ASSERT.h>
#include <cmath>
#include <limits>

namespace PhysBAM{
namespace ADAPTIVE_DETAIL{

template<class T_ADAPTIVE0,class T_ADAPTIVE1>
struct IS_ADAPTIVE<ADAPTIVE_QUOTIENT<T_ADAPTIVE0,T_ADAPTIVE1> > {static const bool value=true;};

template<class T_ADAPTIVE0,class T_ADAPTIVE1>
class ADAPTIVE_QUOTIENT:public ADAPTIVE_BASE<ADAPTIVE_QUOTIENT<T_ADAPTIVE0,T_ADAPTIVE1>,typename EXACT_ARITHMETIC_POLICY<typename T_ADAPTIVE0::EXACT_TYPE,typename T_ADAPTIVE1::EXACT_TYPE>::QUOTIENT_TYPE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
public:
    typedef ADAPTIVE_BASE<ADAPTIVE_QUOTIENT,typename EXACT_ARITHMETIC_POLICY<typename T_ADAPTIVE0::EXACT_TYPE,typename T_ADAPTIVE1::EXACT_TYPE>::QUOTIENT_TYPE> BASE;

    typedef typename BASE::EXACT_TYPE EXACT_TYPE;
    typedef typename BASE::FP_TYPE FP_TYPE;

private:
    const T_ADAPTIVE0 a1;
    const T_ADAPTIVE1 a2;

    ADAPTIVE_QUOTIENT();
    ADAPTIVE_QUOTIENT& operator=(const ADAPTIVE_QUOTIENT&);

public:
    ADAPTIVE_QUOTIENT(const T_ADAPTIVE0& a1_,const T_ADAPTIVE1& a2_)
        :a1(a1_),a2(a2_)
    {}

    PAIR<FP_TYPE,FP_TYPE>
    Estimate_And_Error_Implementation() const
    {const FP_TYPE epsilon=(FP_TYPE).5*std::numeric_limits<FP_TYPE>::epsilon();
    FP_TYPE a2_estimate,a2_error;a2.Estimate_And_Error().Get(a2_estimate,a2_error);
    FP_TYPE abs_a2_estimate=fabs(a2_estimate);
    FP_TYPE lower_abs_a2_estimate=abs_a2_estimate-a2_error;
    FP_TYPE estimate,error;
    if(lower_abs_a2_estimate<=0){estimate=0;error=std::numeric_limits<FP_TYPE>::max();}
    else{
        FP_TYPE a1_estimate,a1_error;a1.Estimate_And_Error().Get(a1_estimate,a1_error);
        estimate=a1_estimate/a2_estimate;
        FP_TYPE part_of_the_error =(fabs(a1_estimate)*a2_error+abs_a2_estimate*a1_error)/(abs_a2_estimate*lower_abs_a2_estimate);
        part_of_the_error+=6*epsilon*part_of_the_error;
        error=epsilon*fabs(estimate)+part_of_the_error;
        error+=epsilon*error;}
    return PAIR<FP_TYPE,FP_TYPE>(estimate,error);}

    EXACT_TYPE Exact_Implementation() const
    {return (a1.Exact()/a2.Exact());}

    int Quick_Sign_Implementation() const
    {int sign=a1.Quick_Sign()*a2.Quick_Sign();
    return (sign<=1&&sign>=-1?sign:ADAPTIVE_SIGN_UNKNOWN);}

//#####################################################################
};
}
using ADAPTIVE_DETAIL::ADAPTIVE_QUOTIENT;
}
#endif
