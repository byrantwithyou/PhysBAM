#include <Tools/Auto_Diff/AUTO_HESS.h>
#include <Tools/Math_Tools/sqr.h>
#include <Tools/Matrices/MATRIX.h>
#include <Geometry/Basic_Geometry/SMOOTH_GEAR.h>
#include <algorithm>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> void Add_Debug_Particle(const TV& X, const VECTOR<typename TV::SCALAR,3>& color);}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SMOOTH_GEAR<VECTOR<T,2> >::
SMOOTH_GEAR(T R,T S,int N)
    :r(R),s(S),n(N)
{
    Compute_Centers();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SMOOTH_GEAR<VECTOR<T,2> >::
SMOOTH_GEAR(const TV& dimensions,int N)
    :r(dimensions.x),s(dimensions.y),n(N)
{
    Compute_Centers();
}
//#####################################################################
// Function Compute_Centers
//#####################################################################
template<class T> void SMOOTH_GEAR<VECTOR<T,2> >::
Compute_Centers()
{
    T u=cos((T)pi/n);
    T r2=r*r,s2=s*s,pl=r2+s2,rm=r2-s2,q=rm/u;
    PHYSBAM_ASSERT(pl>q);
    T y=sqrt(pl*pl-q*q);
    ci=sqrt(pl-y);
    co=q/ci;
    Ci=TV(ci,0);
    Co=co*TV(cos(pi/n),sin(pi/n));
    den=2*pi/n;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,2> > SMOOTH_GEAR<VECTOR<T,2> >::
Bounding_Box() const
{
    TV ext(ci+s,ci+s);
    return RANGE<TV>(-ext,ext);
}
//#####################################################################
// Function Compute_Helper
//#####################################################################
template<class T> void SMOOTH_GEAR<VECTOR<T,2> >::
Compute_Helper(const TV& X,HELPER& h) const
{
    h.a=::std::atan2(X.y,X.x);
    h.d=X.Magnitude();
    h.k=floor(h.a/den);
    h.b=h.a-h.k*den;
    h.flip=false;
    if(2*h.b>den){h.flip=true;h.b=den-h.b;}
    h.Y=TV(cos(h.b),sin(h.b))*h.d;
    h.ui=Co.y*(h.Y.x-ci)>h.Y.y*(Co.x-ci);
    h.P=h.ui?Ci:Co;
    h.dY=h.Y-h.P;
    h.m=h.dY.Magnitude();
    h.sd=h.ui?h.m-s:s-h.m;
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T SMOOTH_GEAR<VECTOR<T,2> >::
Signed_Distance(const TV& X) const
{
    HELPER h;
    Compute_Helper(X,h);
    return h.sd;
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,2> SMOOTH_GEAR<VECTOR<T,2> >::
Surface(const TV& X,const HELPER& h) const
{
    T a=(h.k+h.flip)*den;
    T cs=cos(a),sn=sin(a);
    MATRIX<T,2> R(cs,sn,-sn,cs);
    TV e=h.P+s/h.m*h.dY;
    if(h.flip) e.y=-e.y;
    return R*e;
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,2> SMOOTH_GEAR<VECTOR<T,2> >::
Surface(const TV& X) const
{
    HELPER h;
    Compute_Helper(X,h);
    return Surface(X,h);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,2> SMOOTH_GEAR<VECTOR<T,2> >::
Normal(const TV& X,const HELPER& h) const
{
    T a=(h.k+h.flip)*den;
    T cs=cos(a),sn=sin(a);
    MATRIX<T,2> R(cs,sn,-sn,cs);
    TV d=h.dY/h.m;
    if(!h.ui) d=-d;
    if(h.flip) d.y=-d.y;
    return R*d;
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,2> SMOOTH_GEAR<VECTOR<T,2> >::
Normal(const TV& X) const
{
    HELPER h;
    Compute_Helper(X,h);
    return Normal(X,h);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,2> SMOOTH_GEAR<VECTOR<T,2> >::
Normal(const TV& X,const int aggregate) const
{
    return Normal(X);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,2> SMOOTH_GEAR<VECTOR<T,2> >::
Hessian(const TV& X) const
{
    AUTO_HESS<T,TV> X0=AUTO_HESS<T,TV>::From_Var(X,0),X1=AUTO_HESS<T,TV>::From_Var(X,1);
    AUTO_HESS<T,TV> a=atan2(X1,X0);
    AUTO_HESS<T,TV> d=hypot(X0,X1);
    T k=floor(a.x/den);
    AUTO_HESS<T,TV> b=a-k*den;
    if(2*b.x>den) b=den-b;
    AUTO_HESS<TV,TV> Y=AUTO_HESS<TV,TV>(cos(b),sin(b))*d;
    bool ui=Co.y*(Y(0).x-ci)>Y(1).x*(Co.x-ci);
    TV P=ui?Ci:Co;
    AUTO_HESS<TV,TV> dY=Y-P;
    AUTO_HESS<T,TV> m=dY.Magnitude();
    AUTO_HESS<T,TV> sd=ui?m-s:s-m;
    return sd.ddx;
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,1> SMOOTH_GEAR<VECTOR<T,2> >::
Principal_Curvatures(const TV& X) const
{
    HELPER h;
    Compute_Helper(X,h);
    return VECTOR<T,1>(h.ui?1/s:-1/s);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,2> >::
Lazy_Inside(const TV& X) const
{
    return Signed_Distance(X)<0;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,2> >::
Lazy_Outside(const TV& X) const
{
    return !Lazy_Inside(X);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,2> >::
Inside(const TV& X,const T thickness_over_two) const
{
    return Signed_Distance(X)<=-thickness_over_two;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,2> >::
Outside(const TV& X,const T thickness_over_two) const
{
    return Signed_Distance(X)>=thickness_over_two;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,2> >::
Boundary(const TV& X,const T thickness_over_two) const
{
    T sd=Signed_Distance(X);
    return (sd<thickness_over_two && sd>-thickness_over_two);
}
//#####################################################################
// Function Name
//#####################################################################
template<class T> std::string SMOOTH_GEAR<VECTOR<T,2> >::
Name()
{
    return "SMOOTH_GEAR<VECTOR<T,2> >";
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SMOOTH_GEAR<VECTOR<T,3> >::
SMOOTH_GEAR(T R,T S,int N,T W)
    :g(R,S,N),w(W)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SMOOTH_GEAR<VECTOR<T,3> >::
SMOOTH_GEAR(const TV& dimensions,int N)
    :g(dimensions.x,dimensions.y,N),w(dimensions.z)
{
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > SMOOTH_GEAR<VECTOR<T,3> >::
Bounding_Box() const
{
    TV ext(g.ci+g.s,g.ci+g.s,w);
    return RANGE<TV>(-ext,ext);
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class T> T SMOOTH_GEAR<VECTOR<T,3> >::
Signed_Distance(const TV& X) const
{
    T ds=g.Signed_Distance(X.Remove_Index(2));
    T dw=fabs(X.z)-w;
    if(dw<0) return std::max(ds,dw);
    if(ds<0) return dw;
    return sqrt(dw*dw+ds*ds);
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> SMOOTH_GEAR<VECTOR<T,3> >::
Surface(const TV& X) const
{
    typename GEAR::HELPER h;
    g.Compute_Helper(X.Remove_Index(2),h);
    VECTOR<T,2> S2=g.Surface(X.Remove_Index(2),h);
    T dw=fabs(X.z)-w,ds=h.sd;
    if(ds>dw && dw<0) return S2.Append(X.z);
    if(ds<0){
        TV Y(X);
        Y.z=X.z>0?w:-w;
        return Y;}
    return S2.Append(X.z>0?w:-w);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> SMOOTH_GEAR<VECTOR<T,3> >::
Normal(const TV& X) const
{
    typename GEAR::HELPER h;
    g.Compute_Helper(X.Remove_Index(2),h);
    T dw=fabs(X.z)-w,ds=h.sd;
    if(ds>dw && dw<0) return g.Normal(X.Remove_Index(2),h).Append(0);
    if(ds<0) return TV(0,0,X.z>0?1:-1);
    return (X-g.Surface(X.Remove_Index(2),h).Append(X.z>0?w:-w)).Normalized();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> SMOOTH_GEAR<VECTOR<T,3> >::
Normal(const TV& X,const int aggregate) const
{
    return Normal(X);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class T> SYMMETRIC_MATRIX<T,3> SMOOTH_GEAR<VECTOR<T,3> >::
Hessian(const TV& X) const
{
    AUTO_HESS<T,TV> X0=AUTO_HESS<T,TV>::From_Var(X,0),X1=AUTO_HESS<T,TV>::From_Var(X,1),X2=AUTO_HESS<T,TV>::From_Var(X,2);
    AUTO_HESS<T,TV> a=atan2(X1,X0);
    AUTO_HESS<T,TV> d=hypot(X0,X1);
    T k=floor(a.x/g.den);
    AUTO_HESS<T,TV> b=a-k*g.den;
    if(2*b.x>g.den) b=g.den-b;
    AUTO_HESS<TV,TV> Y=AUTO_HESS<TV,TV>(cos(b)*d,sin(b)*d,AUTO_HESS<T,TV>::From_Const(0));
    bool ui=g.Co.y*(Y(0).x-g.ci)>Y(1).x*(g.Co.x-g.ci);
    TV P(ui?g.Ci:g.Co);
    AUTO_HESS<TV,TV> dY=Y-P;
    AUTO_HESS<T,TV> m=dY.Magnitude();
    AUTO_HESS<T,TV> sd=ui?m-g.s:g.s-m;
    AUTO_HESS<T,TV> wd=abs(X2)-w;
    if(wd.x<0) return max(sd,wd).ddx;
    if(sd.x<0) return wd.ddx;
    return hypot(wd,sd).ddx;
}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> SMOOTH_GEAR<VECTOR<T,3> >::
Principal_Curvatures(const TV& X) const
{
    typename GEAR::HELPER h;
    g.Compute_Helper(X.Remove_Index(2),h);
    T dw=fabs(X.z)-w,ds=h.sd;
    if(ds>dw && dw<0) return h.ui?VECTOR<T,2>(0,1/g.s):VECTOR<T,2>(-1/g.s,0);
    return VECTOR<T,2>(0,0);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,3> >::
Lazy_Inside(const TV& X) const
{
    return Signed_Distance(X)<0;
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,3> >::
Lazy_Outside(const TV& X) const
{
    return !Lazy_Inside(X);
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,3> >::
Inside(const TV& X,const T thickness_over_two) const
{
    return Signed_Distance(X)<=-thickness_over_two;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,3> >::
Outside(const TV& X,const T thickness_over_two) const
{
    return Signed_Distance(X)>=thickness_over_two;
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool SMOOTH_GEAR<VECTOR<T,3> >::
Boundary(const TV& X,const T thickness_over_two) const
{
    T sd=Signed_Distance(X);
    return (sd<thickness_over_two && sd>-thickness_over_two);
}
//#####################################################################
// Function Name
//#####################################################################
template<class T> std::string SMOOTH_GEAR<VECTOR<T,3> >::
Name()
{
    return "SMOOTH_GEAR<VECTOR<T,3> >";
}
namespace PhysBAM{
template class SMOOTH_GEAR<VECTOR<float,2> >;
template class SMOOTH_GEAR<VECTOR<float,3> >;
template class SMOOTH_GEAR<VECTOR<double,2> >;
template class SMOOTH_GEAR<VECTOR<double,3> >;
}
