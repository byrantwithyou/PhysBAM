#include "helper.h"

void test_22()
{
    TEST(((TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))+(+(TV(0.11,0.65,-0.89)/(abs(cube(b))+1)))));
    TEST(((TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))*((+(u/(abs(a)+1))).Magnitude_Squared())));
    TEST(((+(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))*(u.Dot((u+(-((+v)-v)))))));
    TEST(((+(-(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1)))))-(u+(-((+v)-v)))));
    TEST(-u);
    TEST(+u);
    TEST(u);
    TEST(u*a);
    TEST(u/(a+3));
    TEST((+(+(u/(abs(a)+1)))));
    TEST((+(u/(abs(a)+1))));
    TEST((u/(abs(a)+1)));
    TEST(((u/(abs(a)+1))*atan2(((a*(u/(abs(a)+1))).Dot(v)),hypot(a,0.76))));
    TEST(((u/(abs(a)+1))-(cos(tan(0.88))*(TV(-0.6,-0.1,-0.2)-(u/(abs(a)+1))))));
    TEST(((u/(abs(a)+1)).Cross((u/(abs(a)+1)))));
    TEST((((u/(abs(a)+1)).Cross((u/(abs(a)+1)))).Magnitude_Squared()));
    TEST(((u/(abs(a)+1)).Dot((TV(0.11,0.65,-0.89)/(abs(cube(b))+1)))));
    TEST(((u/(abs(a)+1)).Magnitude()));
    TEST(((+(u/(abs(a)+1))).Magnitude_Squared()));
    TEST(((u/(abs(a)+1)).Magnitude_Squared()));
}
