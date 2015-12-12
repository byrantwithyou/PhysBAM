#include "helper.h"

TV o=TV::All_Ones_Vector();
void test_29()
{
    TEST(sin(v));
    TEST(cos(v));
    TEST(exp(v.Dot(o)*o));
    TEST(cos(v.Dot(o)*o));
    TEST(sin(v.Dot(o)*o));
    TEST(sin(v+o));
    TEST(exp(v.Dot(o)*o));
    TEST(exp(o*a));
    TEST(exp(v+v));
}
