#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
namespace PhysBAM{
//#####################################################################
// Function Test
//#####################################################################
template<class R,class T1> void NONLINEAR_FUNCTION<R(T1)>::
Test(const T1 x,bool test_second_diff)
{
    T1 eps=(T1)1e-6;
    RANDOM_NUMBERS<T1> random;
    T1 dx=random.Get_Uniform_Number(eps/2,eps);

    R e0=0,e1=0,g0=0,g1=0,h0=0,h1=0;
    Compute(x,test_second_diff?&h0:0,&g0,&e0);
    Compute(x+dx,test_second_diff?&h1:0,&g1,&e1);
    R test0a=(g1+g0)/2;
    R test0b=(e1-e0)/dx;
    R test0=(test0b-test0a)/max(abs(test0a),(R)1e-20);
    R test1a=(h1+h0)/2;
    R test1b=(g1-g0)/dx;
    R test1=(test1b-test1a)/max(abs(test1a),(R)1e-20);

    LOG::cout<<"first diff test "<<test0a<<"    "<<test0b<<"    "<<test0<<std::endl;
    if(test_second_diff) LOG::cout<<"second diff test "<<test1a<<"    "<<test1b<<"    "<<test1<<std::endl;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>::
~NONLINEAR_FUNCTION()
{
}
//#####################################################################
// operator()
//#####################################################################
template<class T> T NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>::
operator()(const KRYLOV_VECTOR_BASE<T>& x) const
{
    T E=0;
    Compute(x,0,0,&E);
    return E;
}
//#####################################################################
// Function Test
//#####################################################################
template<class T> void NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>::
Test(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>& h) const
{
    T eps=(T)1e-6;
    RANDOM_NUMBERS<T> random;
    KRYLOV_VECTOR_BASE<T> *dx=x.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *g0=x.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *g1=x.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *a=x.Clone_Default();
    KRYLOV_VECTOR_BASE<T> *b=x.Clone_Default();

    for(int i=0,n=dx->Raw_Size();i<n;i++)
        dx->Raw_Get(i)=random.Get_Uniform_Number(-eps,eps);
    h.Project(*dx);

    T e0=0,e1=0;
    Compute(x,&h,g0,&e0);
    h.Multiply(*dx,*a);
    h.Project(*a);

    // new x
    b->Copy(1,x,*dx);
    Compute(*b,&h,g1,&e1);
    h.Multiply(*dx,*b);
    h.Project(*b);

    T test0a=(h.Inner_Product(*g1,*dx)+h.Inner_Product(*g0,*dx))/(2*eps);
    T test0b=(e1-e0)/eps;
    T test0=(test0b-test0a)/max(abs(test0a),(T)1e-20);

    a->Copy(1,*b,*a);
    a->Copy((T).5,*a);
    T test1a=sqrt(h.Inner_Product(*a,*a))/eps;

    b->Copy(-1,*g0,*g1);
    T test1b=sqrt(h.Inner_Product(*b,*b))/eps;

    a->Copy(-1,*b,*a);
    T test1=sqrt(h.Inner_Product(*a,*a))/(eps*max(abs(test1a),(T)1e-20));

    LOG::cout<<"energy diff test "<<test0a<<"    "<<test0b<<"    "<<test0<<std::endl;
    LOG::cout<<"force diff test "<<test1a<<"    "<<test1b<<"    "<<test1<<std::endl;

    delete dx;
    delete g0;
    delete g1;
    delete a;
    delete b;
}
//#####################################################################
// Function Make_Feasible
//#####################################################################
template<class T> void NONLINEAR_FUNCTION<T(KRYLOV_VECTOR_BASE<T>&)>::
Make_Feasible(KRYLOV_VECTOR_BASE<T>& x) const
{
}
template class NONLINEAR_FUNCTION<double(KRYLOV_VECTOR_BASE<double>&)>;
template class NONLINEAR_FUNCTION<float(KRYLOV_VECTOR_BASE<float>&)>;
template class NONLINEAR_FUNCTION<double(double)>;
template class NONLINEAR_FUNCTION<float(float)>;
}
