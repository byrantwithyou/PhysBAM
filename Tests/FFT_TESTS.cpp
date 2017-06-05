//#####################################################################
// Copyright 2012, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFT_TESTS
//#####################################################################

#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_BASE.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE_ITERATOR.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Utilities/TEST_BASE.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Auto_Diff/AUTO_DIFF.h>
#include <Grid_Tools/Fourier_Transforms/FFT.h>
#include <Grid_Tools/Fourier_Transforms/FFT_UTILITIES.h>
namespace PhysBAM{

static bool fft_test_ok=false;

// RANDOM_NUMBERS does not natively support std::complex, so use helpers
template<class T> inline void init(T& x)
{
    static RANDOM_NUMBERS<T> random_;
    x=random_.Get_Uniform_Number(-1,1);
}
template<class T> inline void init(std::complex<T>& x)
{
    T r,i;
    init(r);
    init(i);
    x=std::complex<T>(r,i);
}

template<class T,int d> inline void Symmetry_Test(VECTOR<int,d> size)
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
    ARRAY<std::complex<T>,TV_INT> a(size),b(size);
    for(RANGE_ITERATOR<TV::m> it(a.domain);it.Valid();it.Next())
        init(a(it.index));
    Enforce_Real_Valued_Symmetry(a);
    FFT<TV> fft;
    fft.Inverse_Transform(a,b);
    T eps=std::numeric_limits<T>::epsilon()*100;
    for(RANGE_ITERATOR<TV::m> it(b.domain);it.Valid();it.Next())
        if(fabs(imag(b(it.index)))>eps){
            LOG::printf("FAILED: ifft of Enforce_Real_Valued_Symmetry is real.  size=%P\n",size);
            fft_test_ok=false;
            return;}
}

// Need smooth periodic functions to take the first derivative of,
// as well as the derivative, for 1D, 2D, and 3D.  Use autodiff.
template<class T> inline T func(VECTOR<T,1> X,VECTOR<T,1>& df)
{
    typedef VECTOR<T,1> TV;
    auto x=AUTO_DIFF<T,TV>::From_Var(X,0);
    auto w=exp(sin(x*2))+sqrt(2+sin(x));
    df=w.dx;
    return w.x;
}

template<class T> inline T func(VECTOR<T,2> X,VECTOR<T,2>& df)
{
    typedef VECTOR<T,2> TV;
    auto x=AUTO_DIFF<T,TV>::From_Var(X,0);
    auto y=AUTO_DIFF<T,TV>::From_Var(X,1);
    auto w=exp(sin(x*2))*log(4+cos(y))+sqrt(2+sin(x+2*y));
    df=w.dx;
    return w.x;
}

template<class T> inline T func(VECTOR<T,3> X,VECTOR<T,3>& df)
{
    typedef VECTOR<T,3> TV;
    auto x=AUTO_DIFF<T,TV>::From_Var(X,0);
    auto y=AUTO_DIFF<T,TV>::From_Var(X,1);
    auto z=AUTO_DIFF<T,TV>::From_Var(X,2);
    auto w=exp(sin(x*2))*log(4+cos(y))+sqrt(2+sin(x+2*y))+atan2(2+sin(z),1);
    df=w.dx;
    return w.x;
}

template<class T,class T_ARRAY> inline T
Magnitude(const ARRAY_BASE<std::complex<T>,T_ARRAY>& a)
{
    T x=0;
    for(int i=0;i<a.Size();i++)
        x+=norm(a(i));
    return sqrt(x);
}

template<class T,class T_ARRAY> inline T
Magnitude(const ARRAY_BASE<T,T_ARRAY>& a){return a.Magnitude();}

template<class T,int d> inline void Test_First_Derivatives(VECTOR<int,d> size)
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
    ARRAY<T,TV_INT> a(size);
    ARRAY<std::complex<T>,TV_INT> b(size);
    VECTOR<ARRAY<std::complex<T>,VECTOR<int,d> >,d> c;
    VECTOR<ARRAY<T,VECTOR<int,d> >,d> g;
    c.Fill(b);
    g.Fill(a);
    TV offset,scale;
    for(int i=0;i<d;i++) init(offset(i));
    for(int i=0;i<d;i++) init(scale(i));
    scale+=2;

    GRID<TV> grid(size,RANGE<TV>(offset,2*(T)pi+offset),true);

    TV df,tdf;
    for(RANGE_ITERATOR<TV::m> it(a.domain);it.Valid();it.Next())
        a(it.index)=func(grid.Center(it.index),df);
    FFT<TV> fft;
    fft.Transform(a,b);
    First_Derivatives(b,c,grid.domain.Edge_Lengths()*scale);
    for(int i=0;i<d;i++) fft.Inverse_Transform(c(i),g(i));
    T max_error=0;
    for(RANGE_ITERATOR<TV::m> it(a.domain);it.Valid();it.Next()){
        func(grid.Center(it.index),df);
        df/=decltype(df)(scale);
        TV error;
        for(int i=0;i<d;i++) error(i)+=abs(df(i)-g(i)(it.index)),tdf(i)=g(i)(it.index);
        max_error=std::max(Magnitude(error),max_error);}
    if(max_error>1e-6){
        LOG::printf("FAILED: fft first derivative test.  size=%P max_error=%P\n",size,max_error);
        fft_test_ok=false;}
}

// Need a divergence-free velocity field and a purely divergent velocity field
// for the projection test.  2D and 3D.  Use autodiff.
template<class T> inline VECTOR<T,3> div_free(VECTOR<T,3> X)
{
    typedef VECTOR<T,3> TV;
    auto x=AUTO_DIFF<T,TV>::From_Var(X,0);
    auto y=AUTO_DIFF<T,TV>::From_Var(X,1);
    auto z=AUTO_DIFF<T,TV>::From_Var(X,2);
    auto w=TV(1,2,0)*exp(sin(x*2))*log(4+cos(y))+TV(-1,2,1)*sqrt(2+sin(x+2*y))+TV(0,1,1)*atan2(2+sin(z),1);
    return w.dx.Contract_Permutation_Tensor();
}

template<class T> inline VECTOR<T,2> div_free(VECTOR<T,2> X)
{
    typedef VECTOR<T,2> TV;
    auto x=AUTO_DIFF<T,TV>::From_Var(X,0);
    auto y=AUTO_DIFF<T,TV>::From_Var(X,1);
    auto w=exp(sin(x*2))*log(4+cos(y))+sqrt(2+sin(x+2*y));
    return w.dx.Orthogonal_Vector();
}

template<class T> inline VECTOR<T,3> pure_div(VECTOR<T,3> X)
{
    typedef VECTOR<T,3> TV;
    auto x=AUTO_DIFF<T,TV>::From_Var(X,0);
    auto y=AUTO_DIFF<T,TV>::From_Var(X,1);
    auto z=AUTO_DIFF<T,TV>::From_Var(X,2);
    auto w=exp(sin(x*3))*log(2+cos(y))+sqrt(2.2+sin(x*3+2*y))+atan2(2.2+sin(z),1.4);
    return w.dx;
}

template<class T> inline VECTOR<T,2> pure_div(VECTOR<T,2> X)
{
    typedef VECTOR<T,2> TV;
    auto x=AUTO_DIFF<T,TV>::From_Var(X,0);
    auto y=AUTO_DIFF<T,TV>::From_Var(X,1);
    auto w=exp(sin(x*3))*log(2+cos(y))+sqrt(2.2+sin(x*3+2*y));
    return w.dx;
}

template<class T,int d> inline void Test_Make_Div_Free(VECTOR<int,d> size)
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;
    VECTOR<ARRAY<T,VECTOR<int,d> >,d> a,c;
    VECTOR<ARRAY<std::complex<T>,VECTOR<int,d> >,d> b;
    a.Fill(ARRAY<T,TV_INT>(size));
    b.Fill(ARRAY<std::complex<T>,TV_INT>(size));
    c=a;
    TV offset;
    for(int i=0;i<d;i++) init(offset(i));

    RANGE<TV> range(offset,2*(T)pi+offset);
    TV scale=TV()+1;
    scale.x=2;
    range.Scale_About_Center(scale);
    GRID<TV> grid(size,range,true);

    for(RANGE_ITERATOR<TV::m> it(a(0).domain);it.Valid();it.Next()){
        TV u=pure_div(grid.Center(it.index))+div_free(grid.Center(it.index));
        for(int i=0;i<d;i++) a(i)(it.index)=u(i);}
    FFT<TV> fft;
    for(int i=0;i<d;i++) fft.Transform(a(i),b(i));
    Make_Divergence_Free(b,grid.domain.Edge_Lengths());
    for(int i=0;i<d;i++) fft.Inverse_Transform(b(i),c(i));
    T max_error=0;
    for(RANGE_ITERATOR<TV::m> it(a(0).domain);it.Valid();it.Next())
    {
        TV sol=div_free(grid.Center(it.index)),num;
        for(int i=0;i<d;i++) num(i)=c(i)(it.index);
        T error=Magnitude(num-sol);
        max_error=std::max(error,max_error);
    }
    if(max_error>2e-3){
        LOG::printf("FAILED: fft div free test.  size=%P max_error=%P\n",size,max_error);
        fft_test_ok=false;}
}

template<class T> inline void Test_Make_Div_Free(VECTOR<int,1> size){}

template<class T,int d,class I> inline void Test()
{
    typedef VECTOR<T,d> TV;
    typedef VECTOR<int,d> TV_INT;

    for(RANGE_ITERATOR<d> it(TV_INT()+6);it.Valid();it.Next())
    {
        FFT<TV> fft;
        ARRAY<I,TV_INT> a(it.index),ta,c(it.index);
        ARRAY<std::complex<T>,TV_INT> b(it.index);
        for(RANGE_ITERATOR<TV::m> it(a.domain);it.Valid();it.Next())
            init(a(it.index));
        T aa=Magnitude(a.array);
        ta=a;
        fft.Transform(ta,b);
        T bb=Magnitude(b.array);
        if(it.index.Product()) bb/=sqrt((T)it.index.Product()); // normalize fft
        fft.Inverse_Transform(b,c);
        T err=Magnitude(a.array-c.array);
        T eps=std::numeric_limits<T>::epsilon()*100;
        if(err>eps){
            LOG::printf("FAILED: Transform and Inverse_Transform are inverses.  size=%P error=%P\n",it.index,err);
            fft_test_ok=false;}
        if(fabs(aa-bb)>eps){
            LOG::printf("FAILED: Fourier transform should perserve L2 norm.  size=%P error=%P\n",it.index,err);
            fft_test_ok=false;}
        Symmetry_Test<T>(it.index);
    }
    for(RANGE_ITERATOR<d> it(TV_INT()+2);it.Valid();it.Next())
    {
        Test_First_Derivatives<T>(it.index+60);
        Test_Make_Div_Free<T>(it.index+60);
    }
}

template<class T,class I> inline void Scalar_Test()
{
    for(int size=0;size<6;size++)
    {
        FFT<T> fft;
        ARRAY<I> a(size),ta,c(size);
        ARRAY<std::complex<T> > b(size);
        for(int i=0;i<a.m;i++)
            init(a(i));
        T aa=Magnitude(a);
        ta=a;
        fft.Transform(ta,b);
        T bb=Magnitude(b);
        fft.Inverse_Transform(b,c);
        T err=Magnitude(a-c);
        if(size) bb/=sqrt((T)size);
        T eps=std::numeric_limits<T>::epsilon()*100;
        if(err>eps){
            LOG::printf("FAILED: Transform and Inverse_Transform are inverses.  size=%P error=%P\n",size,err);
            fft_test_ok=false;}
        if(fabs(aa-bb)>eps){
            LOG::printf("FAILED: Fourier transform should perserve L2 norm.  size=%P error=%P\n",size,err);
            fft_test_ok=false;}
    }
}

class FFT_TESTS:public TEST_BASE
{

public:
    FFT_TESTS()
        :TEST_BASE("fft")
    {}

    virtual ~FFT_TESTS(){}

    TEST_RESULT Run_Test(int n)
    {
        fft_test_ok=true;
        Test<double,1,double>();
        Test<double,1,std::complex<double> >();
        Test<double,2,double>();
        Test<double,2,std::complex<double> >();
        Test<double,3,double>();
        Test<double,3,std::complex<double> >();
        Scalar_Test<double,double>();
        Scalar_Test<double,std::complex<double> >();

        return fft_test_ok?success:failure;
    }

    int Number_Of_Tests() const
    {return 1;}

//#####################################################################
};
static FFT_TESTS fft_tests;
}
