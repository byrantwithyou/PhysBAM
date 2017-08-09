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

void Test_MM_mult()
{
    TEST2((a*m2)*m2);
    TEST2((a*m2)*mx);
    TEST2((a*m3)*m3);
    TEST2((a*mx)*m3);
    TEST2((a*m2)*m2);
    TEST2((a*m3)*m3);
    TEST2((a*mx)*m3);
    TEST2((a*s2)*m2);
    TEST2((a*s3)*m3);
    TEST2((a*m2)*s2);
    TEST2((a*m3)*s3);
    TEST2((a*mx)*s3);
    TEST2((a*s2)*s2);
    TEST2((a*s3)*s3);
    TEST2((a*m2)*mx);
    TEST2((a*s2)*mx);
    TEST2((a*m2)*(b*m2));
    TEST2((a*m3)*(b*m3));
    TEST2((a*mx)*(b*m3));
    TEST2((a*s2)*(b*m2));
    TEST2((a*s3)*(b*m3));
    TEST2((a*m2)*(b*s2));
    TEST2((a*m3)*(b*s3));
    TEST2((a*mx)*(b*s3));
    TEST2((a*s2)*(b*s2));
    TEST2((a*s3)*(b*s3));
    TEST2((a*m2)*(b*mx));
    TEST2((a*s2)*(b*mx));
    TEST2(m2*(b*m2));
    TEST2(m3*(b*m3));
    TEST2(mx*(b*m3));
    TEST2(s2*(b*m2));
    TEST2(s3*(b*m3));
    TEST2(m2*(b*s2));
    TEST2(m3*(b*s3));
    TEST2(mx*(b*s3));
    TEST2(s2*(b*s2));
    TEST2(s3*(b*s3));
    TEST2(m2*(b*mx));
    TEST2(s2*(b*mx));

    TEST2((a*m2).Transpose_Times(m2));
    TEST2((a*m2).Transpose_Times(mx));
    TEST2((a*m3).Transpose_Times(m3));
    TEST2((a*mx).Transpose_Times(m2));
    TEST2((a*m2).Transpose_Times(m2));
    TEST2((a*m3).Transpose_Times(m3));
    TEST2((a*mx).Transpose_Times(m2));
    TEST2((a*s2).Transpose_Times(m2));
    TEST2((a*s3).Transpose_Times(m3));
    TEST2((a*m2).Transpose_Times(s2));
    TEST2((a*m3).Transpose_Times(s3));
    TEST2((a*mx).Transpose_Times(s2));
    TEST2((a*s2).Transpose_Times(s2));
    TEST2((a*s3).Transpose_Times(s3));
    TEST2((a*m2).Transpose_Times(mx));
    TEST2((a*s2).Transpose_Times(mx));
    TEST2((a*m2).Transpose_Times(b*m2));
    TEST2((a*m3).Transpose_Times(b*m3));
    TEST2((a*mx).Transpose_Times(b*m2));
    TEST2((a*s2).Transpose_Times(b*m2));
    TEST2((a*s3).Transpose_Times(b*m3));
    TEST2((a*m2).Transpose_Times(b*s2));
    TEST2((a*m3).Transpose_Times(b*s3));
    TEST2((a*mx).Transpose_Times(b*s2));
    TEST2((a*s2).Transpose_Times(b*s2));
    TEST2((a*s3).Transpose_Times(b*s3));
    TEST2((a*m2).Transpose_Times(b*mx));
    TEST2((a*s2).Transpose_Times(b*mx));

    TEST2((a*m2).Times_Transpose(m2));
    TEST2((a*m3).Times_Transpose(mx));
    TEST2((a*m3).Times_Transpose(m3));
    TEST2((a*mx).Times_Transpose(m3));
    TEST2((a*m2).Times_Transpose(m2));
    TEST2((a*m3).Times_Transpose(m3));
    TEST2((a*mx).Times_Transpose(m3));
    TEST2((a*s2).Times_Transpose(m2));
    TEST2((a*s3).Times_Transpose(m3));
    TEST2((a*m2).Times_Transpose(s2));
    TEST2((a*m3).Times_Transpose(s3));
    TEST2((a*mx).Times_Transpose(s3));
    TEST2((a*s2).Times_Transpose(s2));
    TEST2((a*s3).Times_Transpose(s3));
    TEST2((a*m3).Times_Transpose(mx));
    TEST2((a*s3).Times_Transpose(mx));
    TEST2((a*m2).Times_Transpose(b*m2));
    TEST2((a*m3).Times_Transpose(b*m3));
    TEST2((a*mx).Times_Transpose(b*m3));
    TEST2((a*s2).Times_Transpose(b*m2));
    TEST2((a*s3).Times_Transpose(b*m3));
    TEST2((a*m2).Times_Transpose(b*s2));
    TEST2((a*m3).Times_Transpose(b*s3));
    TEST2((a*mx).Times_Transpose(b*s3));
    TEST2((a*s2).Times_Transpose(b*s2));
    TEST2((a*s3).Times_Transpose(b*s3));
    TEST2((a*m3).Times_Transpose(b*mx));
    TEST2((a*s3).Times_Transpose(b*mx));
}
