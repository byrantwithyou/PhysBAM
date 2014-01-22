//#####################################################################
// Copyright 2008,Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.x
//#####################################################################
#ifndef __ADAPTIVE_PRODUCT__
#define __ADAPTIVE_PRODUCT__
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
struct IS_ADAPTIVE<ADAPTIVE_PRODUCT<T_ADAPTIVE0,T_ADAPTIVE1> > {static const bool value=true;};

template<class T_ADAPTIVE0,class T_ADAPTIVE1>
class ADAPTIVE_PRODUCT:public ADAPTIVE_BASE<ADAPTIVE_PRODUCT<T_ADAPTIVE0,T_ADAPTIVE1>,typename EXACT_ARITHMETIC_POLICY<typename T_ADAPTIVE0::EXACT_TYPE,typename T_ADAPTIVE1::EXACT_TYPE>::PRODUCT_TYPE>
{
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE0>::value));
    STATIC_ASSERT((IS_ADAPTIVE<T_ADAPTIVE1>::value));
public:
    typedef ADAPTIVE_BASE<ADAPTIVE_PRODUCT,typename EXACT_ARITHMETIC_POLICY<typename T_ADAPTIVE0::EXACT_TYPE,typename T_ADAPTIVE1::EXACT_TYPE>::PRODUCT_TYPE> BASE;
    typedef typename BASE::EXACT_TYPE EXACT_TYPE;
    typedef typename BASE::FP_TYPE FP_TYPE;

private:
    const T_ADAPTIVE0 a1;
    const T_ADAPTIVE1 a2;

    ADAPTIVE_PRODUCT();
    ADAPTIVE_PRODUCT& operator=(const ADAPTIVE_PRODUCT&);

public:
    ADAPTIVE_PRODUCT(const T_ADAPTIVE0& a1_input,const T_ADAPTIVE1& a2_input)
        :a1(a1_input),a2(a2_input)
    {}

    inline PAIR<FP_TYPE,FP_TYPE>
    Estimate_And_Error_Implementation() const
    {const FP_TYPE epsilon=std::numeric_limits<FP_TYPE>::epsilon()/2;
    FP_TYPE a1_estimate,a1_error;
    a1.Estimate_And_Error().Get(a1_estimate,a1_error);
    FP_TYPE a2_estimate,a2_error;
    a2.Estimate_And_Error().Get(a2_estimate,a2_error);
    FP_TYPE estimate=a1_estimate*a2_estimate;
    FP_TYPE error=epsilon*fabs(estimate)+fabs(a1_estimate)*a2_error+a1_error*fabs(a2_estimate)+a1_error*a2_error;
    error+=5*epsilon*error;
    return PAIR<FP_TYPE,FP_TYPE>(estimate,error);}

    EXACT_TYPE Exact_Implementation() const
    {return (a1.Exact()*a2.Exact());}

    int Quick_Sign_Implementation() const
    {int sign=a1.Quick_Sign()*a2.Quick_Sign();
    return sign<=1 && sign>=-1 ? sign:ADAPTIVE_SIGN_UNKNOWN;}

//#####################################################################
};
}
}
#endif
