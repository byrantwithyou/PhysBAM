#include "helper.h"

void test_24()
{
    TEST((+(+v)));
    TEST((+v));
    TEST(v);
    TEST(((+(+v))/(abs(((u/(abs(a)+1)).Magnitude()))+1)));
    TEST((((+(+v))/(abs(((u/(abs(a)+1)).Magnitude()))+1)).Dot((a*(u/(abs(a)+1))))));
    TEST(((+v).Cross((+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))))));
    TEST((((+v).Cross((+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))))).Magnitude_Squared()));
    TEST(v.Cross(u));
    TEST((v.Cross((+(u/(abs(a)+1))))));
    TEST(((v.Cross((+(u/(abs(a)+1)))))-((-((+v)-v))+((-((+v)-v))/(abs(tan(0.88))+1)))));
    TEST(v.Cross(v));
    TEST((-(+((+v).Cross(((+((+v)-v))-(-((+v)-v))))))));
    TEST(((+v).Cross(((+((+v)-v))-(-((+v)-v))))));
    TEST((+((+v).Cross(((+((+v)-v))-(-((+v)-v)))))));
    TEST((((+v).Cross(((+((+v)-v))-(-((+v)-v))))).Magnitude_Squared()));
    TEST(((+((+v).Cross(((+((+v)-v))-(-((+v)-v))))))-(-(-(-((-((+v)-v))-((+v)-v)))))));
    TEST(((+(-(+v))).Dot((hypot(a,0.76)*u))));
    TEST((v.Dot((TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))));
    TEST(v.Dot(u));
    TEST(v.Dot(v));
}
