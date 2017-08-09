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

void Test_M()
{
    TEST2(m2*a);
    TEST2(m3*a);
    TEST2(mx*a);
    TEST2(s2*a);
    TEST2(s3*a);
    TEST2(a*m2);
    TEST2(a*m3);
    TEST2(a*mx);
    TEST2(a*s2);
    TEST2(a*s3);
    TEST2(a*m2*b);
    TEST2(a*m3*b);
    TEST2(a*mx*b);
    TEST2(a*s2*b);
    TEST2(a*s3*b);
    TEST2(b*(a*m2));
    TEST2(b*(a*m3));
    TEST2(b*(a*mx));
    TEST2(b*(a*s2));
    TEST2(b*(a*s3));
    TEST2(1.2*(a*m2));
    TEST2(1.2*(a*m3));
    TEST2(1.2*(a*mx));
    TEST2(1.2*(a*s2));
    TEST2(1.2*(a*s3));
    TEST2((a*m2)*1.2);
    TEST2((a*m3)*1.2);
    TEST2((a*mx)*1.2);
    TEST2((a*s2)*1.2);
    TEST2((a*s3)*1.2);
    TEST2((a*m2).Normal_Equations_Matrix());
    TEST2((a*m3).Normal_Equations_Matrix());
    TEST2((a*mx).Normal_Equations_Matrix());
    TEST2((a*s2).Normal_Equations_Matrix());
    TEST2((a*s3).Normal_Equations_Matrix());
    TEST2((a*m2).Outer_Product_Matrix());
    TEST2((a*m3).Outer_Product_Matrix());
    TEST2((a*mx).Outer_Product_Matrix());
    TEST2((a*s2).Outer_Product_Matrix());
    TEST2((a*s3).Outer_Product_Matrix());
    TEST2((a*m2).Trace());
    TEST2((a*m3).Trace());
    TEST2((a*s2).Trace());
    TEST2((a*s3).Trace());
}

