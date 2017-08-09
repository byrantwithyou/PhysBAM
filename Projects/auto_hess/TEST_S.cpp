//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include "helper.h"

using namespace PhysBAM;

void Test_S()
{
    TEST(a);
    TEST(b);
    TEST(a+1.3);
    TEST(a-1.3);
    TEST(a*1.3);
    TEST(a/1.3);
    TEST(1.3+a);
    TEST(1.3-a);
    TEST(1.3*a);
    TEST(1.3/a);
    TEST(b+a);
    TEST(b-a);
    TEST(b*a);
    TEST(b/a);

    TEST(sin(b));
    TEST(cos(b));
    TEST(tan(b));
    TEST(log(b));
    TEST(exp(b));
    TEST(abs(b));
    TEST(sqrt(b));
    TEST(sqr(b));

    TEST(max(b,a));
    TEST(max(1.2,a));
    TEST(max(b,1.2));
    TEST(min(b,a));
    TEST(min(1.2,a));
    TEST(min(b,1.2));
    TEST(hypot(b,a));
    TEST(hypot(1.2,a));
    TEST(hypot(b,1.2));

    TEST(hypot(1.2,1.4,a*b));
    TEST(hypot(1.2,a,1.7));
    TEST(hypot(1.2,a,a*b));
    TEST(hypot(b,1.4,1.7));
    TEST(hypot(b,1.4,a*b));
    TEST(hypot(b,a,1.7));
    TEST(hypot(b,a,a*b));
    
    TEST(atan2(b,a));
    TEST(atan2(1.2,a));
    TEST(atan2(b,1.2));
}

