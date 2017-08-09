//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include "helper.h"

using namespace PhysBAM;

template<class T,int d> T max(VECTOR<T,d> x){return x.Max();}
template<class T,int d> T min(VECTOR<T,d> x){return x.Min();}

static TV o(1.2,2.3,1.4);
static VECTOR<T,2> p(1.2,2.3);

void Test_V()
{
    TEST(v);
    TEST(w);
    TEST(cos(u));
    TEST(-v);
    TEST(+v);
    TEST(v*1.3);
    TEST(1.3*v);
    TEST(v+w);
    TEST(v-w);
    TEST(w/v);
    TEST(v*w);
    TEST(log(v+exp(u*2)));
    TEST(v.Dot(w));
    TEST(Dot_Product(v,w));
    TEST(v.Cross(w));
    TEST(Cross_Product(v,w));
    TEST(v.Magnitude());
    TEST(w*v.Magnitude());
    TEST(w/v.Magnitude());
    TEST(o+w);
    TEST(o-w);
    TEST(w/o);
    TEST(o*w);
    TEST(log(o+exp(u*2)));
    TEST(Dot_Product(o,w));
    TEST(Cross_Product(o,w));
    TEST(o.Magnitude_Squared());
    TEST(o.Magnitude());
    TEST(v+o);
    TEST(v-o);
    TEST(o/v);
    TEST(v*o);
    TEST(log(v+exp(u*2)));
    TEST(v.Dot(o));
    TEST(Dot_Product(v,o));
    TEST(v.Cross(o));
    TEST(Cross_Product(v,o));
    TEST(v.Magnitude());
    TEST(w(0));
    TEST(w(1));
    TEST(w(2));
    TEST(z);
    TEST(-z);
    TEST(+z);
    TEST(z*1.3);
    TEST(1.3*z);
    TEST(z.Dot(z*a));
    TEST(Dot_Product(z,z*a));
    TEST(z.Magnitude());
    TEST(w*z.Magnitude());
    TEST(w/z.Magnitude());
    TEST(z+p);
    TEST(z-p);
    TEST(p/z);
    TEST(z*p);
    TEST(log(z+exp(z*a)));
    TEST(z.Dot(p));
    TEST(Dot_Product(z,p));
    TEST(z.Magnitude());
    TEST(max(z));
    TEST(min(z));
    TEST(max(v));
    TEST(min(v));
}

