//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include "helper.h"

using namespace PhysBAM;

static TV o(1.2,2.3,1.4);
static VECTOR<T,2> p(1.2,2.3);
static MATRIX<T,2> m2(1.2,1.3,2,4.1);
static MATRIX<T,3> m3(1.2,1.3,2,4.1,.2,.3,.5,.2,.7);
static MATRIX<T,2,3> mx(1.2,1.3,2,4.1,.2,.3);
static SYMMETRIC_MATRIX<T,3> s3(1.2,1.3,2,4.1,.2,.3);
static SYMMETRIC_MATRIX<T,2> s2(1.2,1.3,2);

void Test_MV()
{
    TEST(m2*z);
    TEST(m3*v);
    TEST(mx*v);
    TEST(s2*z);
    TEST(s3*v);
}

