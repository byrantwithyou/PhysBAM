#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_WRAPPER<T,TV>::
KRYLOV_VECTOR_WRAPPER()
    :deep_copy(false)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_WRAPPER<T,TV>::
KRYLOV_VECTOR_WRAPPER(TV vector)
    :v(vector),deep_copy(false)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T,class TV> template<class VECTOR,class INDICES> KRYLOV_VECTOR_WRAPPER<T,TV>::
KRYLOV_VECTOR_WRAPPER(VECTOR& vector,const INDICES& index)
    :v(vector,index),deep_copy(false)
{
}
//#####################################################################
// Function Destroy_Helper
//#####################################################################
template<class T,class TV> void Destroy_Helper(KRYLOV_VECTOR_WRAPPER<T,TV>& v){}
template<class T,class TV> void Destroy_Helper(KRYLOV_VECTOR_WRAPPER<T,TV&>& v){delete &v.v;}
template<class T,class T2,class I> void Destroy_Helper(KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<T2>,I> >& v){delete &v.v.array;}
template<class T,class T2,class I>
void Destroy_Helper(KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY_VIEW<T2>,I> >& v){delete [] v.v.array.Get_Array_Pointer();}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_WRAPPER<T,TV>::
~KRYLOV_VECTOR_WRAPPER()
{
    if(deep_copy) Destroy_Helper(*this);
}
//#####################################################################
// Operator +=
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_BASE<T>& KRYLOV_VECTOR_WRAPPER<T,TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v+=dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv).v;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_BASE<T>& KRYLOV_VECTOR_WRAPPER<T,TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v-=dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv).v;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_BASE<T>& KRYLOV_VECTOR_WRAPPER<T,TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T,class TV> void KRYLOV_VECTOR_WRAPPER<T,TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    v.Copy(c,dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv).v,v);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class T,class TV> void KRYLOV_VECTOR_WRAPPER<T,TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    v.Copy(c,dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv1).v,dynamic_cast<const KRYLOV_VECTOR_WRAPPER&>(bv2).v);
}
namespace{
inline int Raw_Size_Helper(const float& p){return 1;}
inline int Raw_Size_Helper(const double& p){return 1;}
template<class T,class T_ARRAY> inline int Raw_Size_Helper(const ARRAY_BASE<T,T_ARRAY>& p){if(!p.Size()) return 0;return p.Size()*Raw_Size_Helper(p(0));}
inline float& Raw_Get_Helper(int i,float*,float& p){assert(i==0);return p;}
inline double& Raw_Get_Helper(int i,double*,double& p){assert(i==0);return p;}
template<class S,class T,class T_ARRAY> inline S& Raw_Get_Helper(int i,S* a,ARRAY_BASE<T,T_ARRAY>& p){int s=Raw_Size_Helper(p(0));return Raw_Get_Helper(i%s,a,p(i/s));}
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class T,class TV> int KRYLOV_VECTOR_WRAPPER<T,TV>::
Raw_Size() const
{
    return Raw_Size_Helper(v);
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class T,class TV> T& KRYLOV_VECTOR_WRAPPER<T,TV>::
Raw_Get(int i)
{
    return Raw_Get_Helper(i,(T*)0,v);
}
//#####################################################################
// Function Clone_Default_Helper
//#####################################################################
template<class T,class TV>
static KRYLOV_VECTOR_BASE<T>* Clone_Default_Helper(const KRYLOV_VECTOR_WRAPPER<T,TV>& v)
{
    KRYLOV_VECTOR_WRAPPER<T,TV>* c=new KRYLOV_VECTOR_WRAPPER<T,TV>;
    c->v.Resize(v.v.Size());
    c->deep_copy=true;
    return c;
}
template<class T,class T2>
static KRYLOV_VECTOR_BASE<T>* Clone_Default_Helper(const KRYLOV_VECTOR_WRAPPER<T,ARRAY<ARRAY<T2> > >& v)
{
    typedef ARRAY<ARRAY<T2> > TV;
    KRYLOV_VECTOR_WRAPPER<T,TV>* c=new KRYLOV_VECTOR_WRAPPER<T,TV>;
    c->v.Resize(v.v.Size());
    for(int i=0;i<c->v.m;i++)
        c->v(i).Resize(v.v(i).Size());
    c->deep_copy=true;
    return c;
}
template<class T,class TV>
static KRYLOV_VECTOR_BASE<T>* Clone_Default_Helper(const KRYLOV_VECTOR_WRAPPER<T,TV&>& v)
{
    KRYLOV_VECTOR_WRAPPER<T,TV&>* c=new KRYLOV_VECTOR_WRAPPER<T,TV&>(*new TV(v.v.Size()));
    c->deep_copy=true;
    return c;
}
template<class T,class T2,class I>
static KRYLOV_VECTOR_BASE<T>* Clone_Default_Helper(const KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<T2>,I> >& v)
{
    typedef INDIRECT_ARRAY<ARRAY<T2>,I> TV;
    KRYLOV_VECTOR_WRAPPER<T,TV>* c=new KRYLOV_VECTOR_WRAPPER<T,TV>(*new ARRAY<T2>(v.v.array.m),v.v.indices);
    c->deep_copy=true;
    return c;
}
template<class T,class T2,class I>
static KRYLOV_VECTOR_BASE<T>* Clone_Default_Helper(const KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY_VIEW<T2>,I> >& v)
{
    typedef INDIRECT_ARRAY<ARRAY_VIEW<T2>,I> TV;
    ARRAY_VIEW<T2> view(v.v.array.Size(),new T2[v.v.array.Size()]);
    for(int i=0;i<view.Size();i++) view(i)=T2();
    KRYLOV_VECTOR_WRAPPER<T,TV>* c=new KRYLOV_VECTOR_WRAPPER<T,TV>(view,v.v.indices);
    c->deep_copy=true;
    return c;
}
//#####################################################################
// Function Clone_Default
//#####################################################################
template<class T,class TV> KRYLOV_VECTOR_BASE<T>* KRYLOV_VECTOR_WRAPPER<T,TV>::
Clone_Default() const
{
    return Clone_Default_Helper(*this);
}
//#####################################################################
// Function Resize_Helper
//#####################################################################
template<class T,class TV>
static void Resize_Helper(KRYLOV_VECTOR_WRAPPER<T,TV>& v,const KRYLOV_VECTOR_WRAPPER<T,TV>& w)
{
    v.v.Resize(w.v.Size());
}
template<class T,class T2>
static void Resize_Helper(KRYLOV_VECTOR_WRAPPER<T,ARRAY<ARRAY<T2> > >& v,const KRYLOV_VECTOR_WRAPPER<T,ARRAY<ARRAY<T2> > >& w)
{
    v.v.Resize(w.v.Size());
    for(int i=0;i<w.v.m;i++)
        v.v(i).Resize(w.v(i).Size());
}
template<class T,class TV>
static void Resize_Helper(KRYLOV_VECTOR_WRAPPER<T,TV&>& v,const KRYLOV_VECTOR_WRAPPER<T,TV&>& w)
{
    v.v.Resize(w.v.Size());
}
template<class T,class T2,class I>
static void Resize_Helper(KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<T2>,I> >& v,
    const KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<T2>,I> >& w)
{
    v.v.array.Resize(w.v.array.Size());
}
template<class T,class T2,class I>
static void Resize_Helper(KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY_VIEW<T2>,I> >& v,
    const KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY_VIEW<T2>,I> >& w)
{
    if(v.v.array.Size()==w.v.array.Size()) return;
    if(v.v.array.Size()>w.v.array.Size()){
        ARRAY_VIEW<T2> view(w.v.array.Size(),v.v.array.Get_Array_Pointer());
        v.v.array.Exchange(view);
        return;}

    ARRAY_VIEW<T2> view(w.v.array.Size(),new T2[w.v.array.Size()]);
    view.Fill(T2());
    v.v.array.Exchange(view);
    delete [] view.Get_Array_Pointer();
}
//#####################################################################
// Function Resize
//#####################################################################
template<class T,class TV> void KRYLOV_VECTOR_WRAPPER<T,TV>::
Resize(const KRYLOV_VECTOR_BASE<T>& v)
{
    if(!deep_copy) return;
    Resize_Helper(*this,debug_cast<const KRYLOV_VECTOR_WRAPPER<T,TV>&>(v));
}
template class KRYLOV_VECTOR_WRAPPER<float,ARRAY<ARRAY<float> > >;
template class KRYLOV_VECTOR_WRAPPER<float,ARRAY<float> >;
template KRYLOV_VECTOR_BASE<float>& KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<float,int>,ARRAY<int>&> >::operator*=(float);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,2>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,3>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<float,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<float,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY<float,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,2>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<float,2>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<float,3>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<float,3>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<float,ARRAY<float>&>::KRYLOV_VECTOR_WRAPPER(ARRAY<float>&);
template KRYLOV_VECTOR_WRAPPER<float,ARRAY<float>&>::~KRYLOV_VECTOR_WRAPPER();
template class KRYLOV_VECTOR_WRAPPER<double,ARRAY<ARRAY<double> > >;
template class KRYLOV_VECTOR_WRAPPER<double,ARRAY<double> >;
template KRYLOV_VECTOR_BASE<double>& KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<double,int>,ARRAY<int>&> >::operator*=(double);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,2>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,3>,int>&,ARRAY<int> const&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<double,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<double,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY<double,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,2>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<double,2>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,2>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY<VECTOR<double,3>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3>,int>,ARRAY<int>&> >::KRYLOV_VECTOR_WRAPPER(ARRAY_VIEW<VECTOR<double,3>,int>&,const ARRAY<int>&);
template KRYLOV_VECTOR_WRAPPER<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3>,int>,ARRAY<int>&> >::~KRYLOV_VECTOR_WRAPPER();
template KRYLOV_VECTOR_WRAPPER<double,ARRAY<double>&>::KRYLOV_VECTOR_WRAPPER(ARRAY<double>&);
template KRYLOV_VECTOR_WRAPPER<double,ARRAY<double>&>::~KRYLOV_VECTOR_WRAPPER();
