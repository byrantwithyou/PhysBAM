#include "helper.h"

void test_26()
{
    TEST((-((-((+v)-v))/(abs(tan(0.88))+1))));
    TEST(((-((+v)-v))/(abs(tan(0.88))+1)));
    TEST((+(-((-((+v)-v))/(abs(tan(0.88))+1)))));
    TEST(((-((-((+v)-v))/(abs(tan(0.88))+1))).Magnitude()));
    TEST((((-((-((+v)-v))/(abs(tan(0.88))+1))).Magnitude())*(-(-(a*(u/(abs(a)+1)))))));
    TEST((-((-((-((+v)-v))/(abs(tan(0.88))+1))).Magnitude_Squared())));
    TEST(((-((-((+v)-v))/(abs(tan(0.88))+1))).Magnitude_Squared()));
    TEST(((+(-((+v)-v))).Cross((-(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))))));
    TEST(((-(-(+v)))-((+v).Cross(((+((+v)-v))-(-((+v)-v)))))));
    TEST(((+((+v)-v)).Cross(((-((+v)-v))+((-((+v)-v))/(abs(tan(0.88))+1))))));
    TEST(((-((+v)-v))*cube(b)));
    TEST((-(((+v)-v).Dot((cos(tan(0.88))*(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))))));
    TEST((((+v)-v).Dot((cos(tan(0.88))*(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))))));
    TEST((+(((+v)-v).Dot((cos(tan(0.88))*(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))))));
    TEST(((-((+v)-v)).Dot((TV(0.11,0.65,-0.89)/(abs(cube(b))+1)))));
    TEST(((-((+v)-v)).Magnitude()));
    TEST(((-((+v)-v)).Magnitude_Squared()));
    TEST(((+(-((+v)-v))).Magnitude_Squared()));
    TEST(((+((+v)-v)).Magnitude_Squared()));
    TEST(((-(+(-((+v)-v))))+(((+(+(u/(abs(a)+1))))+((+v)-v))-(+((+v)-v)))));
}
