//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include "helper.h"

using namespace PhysBAM;

static MATRIX<T,2> m2(1.2,1.3,2,4.1);
static MATRIX<T,3> m3(1.2,1.3,2,4.1,.2,.3,.5,.2,.7);
static MATRIX<T,2,3> mx(1.2,1.3,2,4.1,.2,.3);
static SYMMETRIC_MATRIX<T,3> s3(1.2,1.3,2,4.1,.2,.3);
static SYMMETRIC_MATRIX<T,2> s2(1.2,1.3,2);

void Test_MM()
{
    TEST2(a*m2+m2);
    TEST2(a*m3+m3);
    TEST2(a*mx+mx);
    TEST2(a*s2+s2);
    TEST2(a*s3+s3);
    TEST2(a*s2+m2);
    TEST2(a*s3+m3);
    TEST2(a*m2+s2);
    TEST2(a*m3+s3);
    TEST2(a*m2-m2);
    TEST2(a*m3-m3);
    TEST2(a*mx-mx);
    TEST2(a*s2-s2);
    TEST2(a*s3-s3);
    TEST2(a*s2-m2);
    TEST2(a*s3-m3);
    TEST2(a*m2-s2);
    TEST2(a*m3-s3);
    TEST2(m2+a*m2);
    TEST2(m3+a*m3);
    TEST2(mx+a*mx);
    TEST2(s2+a*s2);
    TEST2(s3+a*s3);
    TEST2(s2+a*m2);
    TEST2(s3+a*m3);
    TEST2(m2+a*s2);
    TEST2(m3+a*s3);
    TEST2(m2-a*m2);
    TEST2(m3-a*m3);
    TEST2(mx-a*mx);
    TEST2(s2-a*s2);
    TEST2(s3-a*s3);
    TEST2(s2-a*m2);
    TEST2(s3-a*m3);
    TEST2(m2-a*s2);
    TEST2(m3-a*s3);
    TEST2(a*m2+b*m2);
    TEST2(a*m3+b*m3);
    TEST2(a*mx+b*mx);
    TEST2(a*s2+b*s2);
    TEST2(a*s3+b*s3);
    TEST2(a*s2+b*m2);
    TEST2(a*s3+b*m3);
    TEST2(a*m2+b*s2);
    TEST2(a*m3+b*s3);
    TEST2(a*m2-b*m2);
    TEST2(a*m3-b*m3);
    TEST2(a*mx-b*mx);
    TEST2(a*s2-b*s2);
    TEST2(a*s3-b*s3);
    TEST2(a*s2-b*m2);
    TEST2(a*s3-b*m3);
    TEST2(a*m2-b*s2);
    TEST2(a*m3-b*s3);

    TEST2((a*m2).Double_Contract(b*m2));
    TEST2((a*m3).Double_Contract(b*m3));
    TEST2((a*mx).Double_Contract(b*mx));
    TEST2((a*s2).Double_Contract(b*s2));
    TEST2((a*s3).Double_Contract(b*s3));
    TEST2((a*m2).Double_Contract(b*s2));
    TEST2((a*m3).Double_Contract(b*s3));
    TEST2((a*s2).Double_Contract(b*m2));
    TEST2((a*s3).Double_Contract(b*m3));

    TEST2((a*m2).Double_Contract(m2));
    TEST2((a*m3).Double_Contract(m3));
    TEST2((a*mx).Double_Contract(mx));
    TEST2((a*s2).Double_Contract(s2));
    TEST2((a*s3).Double_Contract(s3));
    TEST2((a*m2).Double_Contract(s2));
    TEST2((a*m3).Double_Contract(s3));
    TEST2((a*s2).Double_Contract(m2));
    TEST2((a*s3).Double_Contract(m3));
}
