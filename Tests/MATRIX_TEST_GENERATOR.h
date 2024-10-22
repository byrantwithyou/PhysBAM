#ifndef __MATRIX_TEST_GENERATOR__
#define __MATRIX_TEST_GENERATOR__
#include "MATRIX_TESTS.h"
namespace PhysBAM{
template<class T2,class M0,class M1,class M2,class TM0,class TM1,class TM2,int a,int b,int c>
struct TEST_GENERATOR1
{
    bool v;
    TEST_GENERATOR1(const MATRIX_TESTS<T2>& m,int n){v=regular(m,n)&square(m,n,VECTOR<int,a==b && (a<=3)>());}

    bool regular(const MATRIX_TESTS<T2>& m,int n)
    {TEST_GENERATOR1<T2,M1,M2,MATRIX<T2,a,b>,TM1,TM2,MATRIX<T2,b,a>,b,c,a> tg1(m,n);
    TEST_GENERATOR1<T2,M1,M2,MATRIX_MXN<T2>,TM1,TM2,MATRIX_MXN<T2>,b,c,a> tg2(m,n);
    return tg1.v && tg2.v;}

    bool square(const MATRIX_TESTS<T2>& m,int n,VECTOR<int,0> w){return true;}

    bool square(const MATRIX_TESTS<T2>& m,int n,VECTOR<int,1> w)
    {TEST_GENERATOR1<T2,M1,M2,SYMMETRIC_MATRIX<T2,a>,TM1,TM2,SYMMETRIC_MATRIX<T2,a>,b,c,a> tg1(m,n);
    TEST_GENERATOR1<T2,M1,M2,DIAGONAL_MATRIX<T2,a>,TM1,TM2,DIAGONAL_MATRIX<T2,a>,b,c,a> tg2(m,n);
    TEST_GENERATOR1<T2,M1,M2,UPPER_TRIANGULAR_MATRIX<T2,a>,TM1,TM2,UPPER_TRIANGULAR_MATRIX<T2,a>,b,c,a> tg3(m,n);
    return tg1.v && tg2.v && tg3.v;}
};

template<class T2,class M1,class M2,class TM1,class TM2,int a,int b,int c>
struct TEST_GENERATOR1<T2,bool,M1,M2,bool,TM1,TM2,a,b,c>
{
    bool v;
    TEST_GENERATOR1(const MATRIX_TESTS<T2>& m,int n)
    {
        v=true;
        for(int i=0;i<n;i++) v=m.Arbitrary_Test_Two_Sizes(M1(INITIAL_SIZE(b),INITIAL_SIZE(c)),M2(INITIAL_SIZE(c),INITIAL_SIZE(a)),n)&v;
        for(int i=0;i<n;i++) v=m.Arbitrary_Test_Two_Sizes_TX(TM1(INITIAL_SIZE(c),INITIAL_SIZE(b)),M2(INITIAL_SIZE(c),INITIAL_SIZE(a)),n)&v;
        for(int i=0;i<n;i++) v=m.Arbitrary_Test_Two_Sizes_XT(M1(INITIAL_SIZE(b),INITIAL_SIZE(c)),TM2(INITIAL_SIZE(a),INITIAL_SIZE(c)),n)&v;
        if(a==b) for(int i=0;i<n;i++) v=m.Arbitrary_Test_One_Size(M1(INITIAL_SIZE(b),INITIAL_SIZE(c)),n)&v;
    }
};

template<class T2,int a,int b,int x> // use x=b*b*b*b-1
struct TEST_GENERATOR
{
    bool v;
    TEST_GENERATOR(const MATRIX_TESTS<T2>& m,int n){v=TEST_GENERATOR1<T2,int,int,bool,int,int,bool,x/b/b+a,x/b%b+a,x%b+a>(m,n).v;v=(TEST_GENERATOR<T2,a,b,x-1>(m,n).v && v);}
};

template<class T2,int a,int b>
struct TEST_GENERATOR<T2,a,b,-1>
{
    bool v;
    TEST_GENERATOR(const MATRIX_TESTS<T2>& m,int n){v=true;}
};
}
#endif
