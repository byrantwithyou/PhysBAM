#include "helper.h"

void test_21()
{
    TEST(((-(+(TV(0.11,0.65,-0.89)/(abs(cube(b))+1)))).Magnitude()));
    TEST(((+(TV(0.11,0.65,-0.89)/(abs(cube(b))+1))).Magnitude()));
    TEST(((+(+(TV(0.11,0.65,-0.89)/(abs(cube(b))+1)))).Magnitude_Squared()));
    TEST(((TV(0.11,0.65,-0.89)/(abs(cube(b))+1)).Magnitude_Squared()));
    TEST((((+(+(TV(0.11,0.65,-0.89)/(abs(cube(b))+1)))).Magnitude_Squared())*(-(-(+v)))));
    TEST(((TV(0.11,0.65,-0.89)/(abs(cube(b))+1))+((-((+v)-v))-((+v)-v))));
    TEST((((TV(0.11,0.65,-0.89)/(abs(cube(b))+1))+((-((+v)-v))-((+v)-v)))+(+((+v)-v))));
    TEST((-(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))));
    TEST((+(-(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))));
    TEST((+(+(+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))))));
    TEST((+(+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))));
    TEST((+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))));
    TEST((TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))));
    TEST(((TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))+(a*(u/(abs(a)+1)))));
    TEST(((+(-(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))).Magnitude()));
    TEST(((+(+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))).Magnitude()));
    TEST(((+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))).Magnitude()));
    TEST((+((+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))).Magnitude())));
    TEST(((+(+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))).Magnitude_Squared()));
    TEST(((+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))).Magnitude_Squared()));
}
