#include <Core/Arrays/ARRAY.h>
#include <Tools/Interpolation/CUBIC_MN_INTERPOLATION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T2> CUBIC_MN_INTERPOLATION<T,T2>::
CUBIC_MN_INTERPOLATION()
{
    Set_Parameters();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T2> CUBIC_MN_INTERPOLATION<T,T2>::
~CUBIC_MN_INTERPOLATION()
{
}
//#####################################################################
// Function Set_Parameters
//#####################################################################
template<class T,class T2> void CUBIC_MN_INTERPOLATION<T,T2>::
Set_Parameters(const T b,const T c)
{
    m0=(12-9*b-6*c)/6;
    m1=(-18+12*b+6*c)/6;
    m3=(6-2*b)/6;
    n0=(-b-6*c)/6;
    n1=(6*b+30*c)/6;
    n2=(-12*b-48*c)/6;
    n3=(8*b+24*c)/6;
}
//#####################################################################
// Cubic_MN
//#####################################################################
template<class T,class T2> T2 CUBIC_MN_INTERPOLATION<T,T2>::
Cubic_MN(const T2& u_0,const T2& u_1,const T2& u_2,const T2& u_3,const T alpha,const bool monotonic) const
{
    if(monotonic) return Cubic_MN_Monotonic(u_0,u_1,u_2,u_3,alpha);
    T temp=abs(1+alpha);
    T2 result=u_0*(((n0*temp+n1)*temp+n2)*temp+n3);
    temp=abs(alpha);result+=u_1*((m0*temp+m1)*temp*temp+m3);
    temp=abs(1-alpha);result+=u_2*((m0*temp+m1)*temp*temp+m3);
    temp=abs(2-alpha);result+=u_3*(((n0*temp+n1)*temp+n2)*temp+n3);
    return result;
}
//#####################################################################
// Cubic_MN_Monotonic
//#####################################################################
template<class T,class T2> T2 CUBIC_MN_INTERPOLATION<T,T2>::
Cubic_MN_Monotonic(const T2& u_0,const T2& u_1,const T2& u_2,const T2& u_3,const T alpha) const
{
    T x2=alpha*alpha;
    T x3=x2*alpha;
    T2 dq=u_2-u_1;
    T2 di=(T).5*(u_2-u_0); if(di*dq<T2()*(T)0) di*=(T)0;
    T2 dip1=(T).5*(u_3-u_1); if(dip1*dq<T2()*(T)0) dip1*=(T)0;
    return u_1+di*alpha+((T)3.*dq-(T)2.*di-dip1)*x2+(-(T)2.*dq+di+dip1)*x3;
}
//#####################################################################
// Function Cubic_MN_Weights
//#####################################################################
template<class T,class T2> ARRAY<T> CUBIC_MN_INTERPOLATION<T,T2>::
Cubic_MN_Weights(const T alpha) const
{
    T alpha2=alpha*alpha;
    T alpha3=alpha2*alpha;
    ARRAY<T> weights;
    weights.Append((T)(-alpha3/2.+alpha2-.5*alpha)); //f_k-1
    weights.Append((T)(3.*alpha3/2.-5*alpha2/2.+1.)); //f_k 
    weights.Append((T)(-3.*alpha3/2.+2*alpha2+alpha/2.)); //f_k+1
    weights.Append((T)(alpha3/2.-alpha2/2.));  //f_k+2
    for(int i=0;i<weights.m;i++) weights(i)=max((T)0.,min((T)1.,weights(i)));
    return weights;
}
namespace PhysBAM{
template class CUBIC_MN_INTERPOLATION<float,float>;
template class CUBIC_MN_INTERPOLATION<float,VECTOR<float,3> >;
template class CUBIC_MN_INTERPOLATION<float,VECTOR<float,2> >;
template class CUBIC_MN_INTERPOLATION<double,double>;
template class CUBIC_MN_INTERPOLATION<double,VECTOR<double,3> >;
template class CUBIC_MN_INTERPOLATION<double,VECTOR<double,2> >;
}
