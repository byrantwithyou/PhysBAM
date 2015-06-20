#include "helper.h"

void test_25()
{
    TEST(((+v)*exp((-tan(0.88)))));
    TEST((-((-(+v)).Magnitude())));
    TEST(((-(+v)).Magnitude()));
    TEST(((+(+v)).Magnitude()));
    TEST((v.Magnitude()));
    TEST(((+v).Magnitude_Squared()));
    TEST((v.Magnitude_Squared()));
    TEST(((v.Magnitude_Squared())-cube(hypot(a,0.76))));
    TEST((v*min(((+v).Magnitude_Squared()),cube(log(abs(b)+1)))));
    TEST(((v*min(((+v).Magnitude_Squared()),cube(log(abs(b)+1)))).Dot((+(-((+v)-v))))));
    TEST(((+v)+(sin(atan2(log(abs(b)+1),hypot(a,0.76)))*((-0.62-a)*v))));
    TEST(((-(+v))*((+(TV(0.11,0.65,-0.89)/(abs(cube(b))+1))).Magnitude())));
    TEST((v+(u/(abs(a)+1))));
    TEST(((+(+v))+(u+(-((+v)-v)))));
    TEST((-((+v)-v)));
    TEST((-(+(-((+v)-v)))));
    TEST(((+v)-v));
    TEST((+(-((+v)-v))));
    TEST((+((+v)-v)));
    TEST((((+v)-v)-(-(((-0.62-a)*v)/(abs(((-((+v)-v)).Magnitude()))+1)))));
}
