#include "helper.h"

void test_27()
{
    TEST((-((-((+v)-v))-(+v))));
    TEST(((-((+v)-v))-(+v)));
    TEST((+((-((+v)-v))-(+v))));
    TEST((-(((-((+v)-v))-(+v))*((-0.62-a)-((a*(u/(abs(a)+1))).Dot(v))))));
    TEST((((-((+v)-v))-(+v))*((-0.62-a)-((a*(u/(abs(a)+1))).Dot(v)))));
    TEST(((+((-((+v)-v))-(+v)))+(-(((a*(u/(abs(a)+1))).Dot(v))*(-((+v)-v))))));
    TEST((((-((+v)-v))-(+v)).Cross((-((-((+v)-v))-((+v)-v))))));
    TEST((-(-(-((-((+v)-v))-((+v)-v))))));
    TEST((-(-((-((+v)-v))-((+v)-v)))));
    TEST((-((-((+v)-v))-((+v)-v))));
    TEST(((-((+v)-v))-((+v)-v)));
    TEST(((+((+v)-v))-(-((+v)-v))));
    TEST(((-((-((+v)-v))-((+v)-v)))/(abs(sin(atan2(log(abs(b)+1),hypot(a,0.76))))+1)));
    TEST(((-((+v)-v))+((-((+v)-v))/(abs(tan(0.88))+1))));
    TEST(((-(-(-((-((+v)-v))-((+v)-v))))).Dot((-(+(-((+v)-v)))))));
    TEST((((-((+v)-v))-((+v)-v))*tan(-0.12)));
    TEST(((((-((+v)-v))-((+v)-v))*tan(-0.12))-(+(-(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))))));
    TEST((-((((-((+v)-v))-((+v)-v))*tan(-0.12))+(-((+v)-v)))));
    TEST(((((-((+v)-v))-((+v)-v))*tan(-0.12))+(-((+v)-v))));
    TEST((((((-((+v)-v))-((+v)-v))*tan(-0.12))+(-((+v)-v)))/(abs(tan(tan(0.88)))+1)));
}
