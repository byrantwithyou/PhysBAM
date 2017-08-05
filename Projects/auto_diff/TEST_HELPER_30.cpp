#include "helper.h"

extern TV o;
void test_30()
{
    TEST(v(0));
    TEST(v(1));
    TEST(cos(v(0))+sin(v(1))*v(2));
    TEST(v*v);
    TEST(o*v);
    TEST(v*o);
    TEST(u*v);
    TEST(sin(u)*exp(v));
    TEST(u/v);
    TEST(o/v);
    TEST(u/o);
    TEST(sin(u)/exp(v));
    TEST(u+Make_Vector<T>(a,b,u.Dot(v)));
}
