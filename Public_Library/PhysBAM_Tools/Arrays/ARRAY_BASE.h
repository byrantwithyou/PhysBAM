//#####################################################################
// Copyright 2004-2009, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_BASE
//#####################################################################
#ifndef __ARRAY_BASE__
#define __ARRAY_BASE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Arrays/SORT.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/maxmag.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/minmag.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/ARITHMETIC_POLICY.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <algorithm>
#include <iostream>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

template<class T_ARRAY,class ENABLER=void> struct CANONICALIZE_CONST_ARRAY:public FIRST<T_ARRAY>{};

template<class T_ARRAY1,class T_ARRAY2> struct SAME_ARRAY_CANONICAL{static bool Same_Array(const T_ARRAY1& array1,const T_ARRAY2& array2)
{STATIC_ASSERT(!IS_SAME<T_ARRAY1,T_ARRAY2>::value);return false;}};

template<class T_ARRAY> struct SAME_ARRAY_CANONICAL<T_ARRAY,T_ARRAY>{static bool Same_Array(const T_ARRAY& array1,const T_ARRAY& array2)
{return T_ARRAY::Same_Array(array1,array2);}};

template<class TA1,class TA2> struct SAME_ARRAY:public SAME_ARRAY_CANONICAL<typename CANONICALIZE_CONST_ARRAY<TA1>::TYPE,typename CANONICALIZE_CONST_ARRAY<TA2>::TYPE>{};

template<class TV> struct ELEMENT_OF_VECTOR {private:struct UNUSABLE;public:typedef UNUSABLE TYPE;};
template<class T,int d> struct ELEMENT_OF_VECTOR<VECTOR<T,d> > {typedef T TYPE;};

template<class T> struct CAN_REFERENCE_ELEMENTS {static const int value=true;};
template<class T_ARRAY1,class T_ARRAY2> struct CAN_REFERENCE_ELEMENTS<ARRAY_SUM<T_ARRAY1,T_ARRAY2> > {static const int value=false;};
template<class T_ARRAY1,class T_ARRAY2> struct CAN_REFERENCE_ELEMENTS<ARRAY_DIFFERENCE<T_ARRAY1,T_ARRAY2> > {static const int value=false;};
template<class T_ARRAY1,class T_ARRAY2> struct CAN_REFERENCE_ELEMENTS<ARRAY_PRODUCT<T_ARRAY1,T_ARRAY2> > {static const int value=false;};
template<class T_ARRAY1,class T_ARRAY2> struct CAN_REFERENCE_ELEMENTS<ARRAY_QUOTIENT<T_ARRAY1,T_ARRAY2> > {static const int value=false;};
template<class T_ARRAY1> struct CAN_REFERENCE_ELEMENTS<ARRAY_NEGATION<T_ARRAY1> > {static const int value=false;};
template<class T_ARRAY1,class T> struct CAN_REFERENCE_ELEMENTS<ARRAY_PLUS_SCALAR<T_ARRAY1,T> > {static const int value=false;};
template<class T_ARRAY1,class T> struct CAN_REFERENCE_ELEMENTS<ARRAY_LEFT_MULTIPLE<T_ARRAY1,T> > {static const int value=false;};

template<class T_ARRAY1,class enabler=void> struct EQUIVALENT_ARRAY;
template<class T> struct EQUIVALENT_ARRAY<T,typename ENABLE_IF<!IS_ARRAY<T>::value>::TYPE> {typedef T TYPE;};
template<class T> struct EQUIVALENT_ARRAY<T,typename ENABLE_IF<IS_ARRAY<T>::value>::TYPE> {typedef ARRAY<typename EQUIVALENT_ARRAY<typename T::ELEMENT>::TYPE> TYPE;};

template<class T,class T_ARRAY,class ID>
class ARRAY_BASE
{
    struct UNUSABLE{};
public:
    typedef T ELEMENT;
    typedef ID INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    typedef const T* const_iterator; // for stl
    typedef int difference_type; // for stl

    typedef typename IF<CAN_REFERENCE_ELEMENTS<T_ARRAY>::value,T&,T>::TYPE T_REF_IF_POSSIBLE;
    typedef typename IF<CAN_REFERENCE_ELEMENTS<T_ARRAY>::value,const T&,const T>::TYPE CONST_T_REF_IF_POSSIBLE;
    typedef T& RESULT_TYPE;
    typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;
    template<class T2> struct REBIND{typedef ARRAY_BASE<T2,typename T_ARRAY::template REBIND<T2>::TYPE,ID> TYPE;};

protected:
    ARRAY_BASE(){}
    ARRAY_BASE(const ARRAY_BASE&){}
    ~ARRAY_BASE(){}
public:

    T_ARRAY& Derived()
    {return static_cast<T_ARRAY&>(*this);}

    const T_ARRAY& Derived() const
    {return static_cast<const T_ARRAY&>(*this);}

protected:
    template<class T_ARRAY2> struct IS_ARRAY_BASE {static const bool value=false;};
    template<class T2,class T_ARRAY2> struct IS_ARRAY_BASE<ARRAY_BASE<T2,T_ARRAY2,ID> > {static const bool value=true;};
public:

    template<class T_ARRAY1> static bool Same_Array(const T_ARRAY1& array1,const T_ARRAY1& array2)
    {return &array1==&array2;}

    template<class T_ARRAY1,class T_ARRAY2> static bool
    Same_Array(const T_ARRAY1& array1,const T_ARRAY2& array2)
    {return SAME_ARRAY<T_ARRAY1,T_ARRAY2>::Same_Array(array1,array2);}

protected:
    T_ARRAY& operator=(const ARRAY_BASE& source)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY& source_=source.Derived();assert(m==source_.Size());
    if(!T_ARRAY::Same_Array(self,source_)) for(ID i(0);i<m;i++) self(i)=source_(i);
    return self;}

    template<class T_ARRAY1>
    T_ARRAY& operator=(const T_ARRAY1& source)
    {STATIC_ASSERT(CAN_ASSIGN<T,typename T_ARRAY1::ELEMENT>::value);
    T_ARRAY& self=Derived();ID m=self.Size();assert(m==source.Size());
    if(!T_ARRAY::Same_Array(self,source)) for(ID i(0);i<m;i++) self(i)=source(i);
    return self;}
public:

    template<class T_INDICES>
    INDIRECT_ARRAY<T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices)
    {return INDIRECT_ARRAY<T_ARRAY,T_INDICES&>(Derived(),indices);}

    template<class T_INDICES>
    INDIRECT_ARRAY<const T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices) const
    {return INDIRECT_ARRAY<const T_ARRAY,T_INDICES&>(Derived(),indices);}

    INDIRECT_ARRAY<T_ARRAY,IDENTITY_ARRAY<> > Prefix(const ID prefix_size)
    {assert(prefix_size<=Derived().Size());return INDIRECT_ARRAY<T_ARRAY,IDENTITY_ARRAY<> >(Derived(),IDENTITY_ARRAY<>(prefix_size));}

    INDIRECT_ARRAY<const T_ARRAY,IDENTITY_ARRAY<> > Prefix(const ID prefix_size) const
    {assert(prefix_size<=Derived().Size());return INDIRECT_ARRAY<const T_ARRAY,IDENTITY_ARRAY<> >(Derived(),IDENTITY_ARRAY<>(prefix_size));}

private:
    typedef typename IF<IS_CLASS<T>::value,T,UNUSABLE>::TYPE T_IF_CLASS;
public:

    template<int d>
    static void Static_Assert_Not_Small(const VECTOR<T,d>*)
    {STATIC_ASSERT(d>3);}

    template<class T_ARRAY1>
    static void Static_Assert_Not_Small(const T_ARRAY1*)
    {}

    void Static_Assert_Not_Small() const
    {Static_Assert_Not_Small((T_ARRAY*)0);}

    template<class T2,class T3,int p,int q>
    static void Assert_Same_Size_Helper(const VECTOR<T2,p>&,const VECTOR<T3,q>&)
    {STATIC_ASSERT(p==q);}

    template<class T2,class T3,class T_ARRAY2,class T_ARRAY3>
    static void Assert_Same_Size_Helper(const ARRAY_BASE<T2,T_ARRAY2>& u,const ARRAY_BASE<T3,T_ARRAY3>& v)
    {assert(u.Size()==v.Size());}

    template<class T2,class T3,class T_ARRAY2,class T_ARRAY3>
    static void Assert_Same_Size(const ARRAY_BASE<T2,T_ARRAY2>& u,const ARRAY_BASE<T3,T_ARRAY3>& v)
    {Assert_Same_Size_Helper(u.Derived(),v.Derived());}

    template<class T_FIELD,T_FIELD T_IF_CLASS::* field>
    PROJECTED_ARRAY<T_ARRAY,FIELD_PROJECTOR<T_IF_CLASS,T_FIELD,field> > Project()
    {return PROJECTED_ARRAY<T_ARRAY,FIELD_PROJECTOR<ELEMENT,T_FIELD,field> >(Derived());}

    template<class T_FIELD,T_FIELD T_IF_CLASS::* field>
    PROJECTED_ARRAY<const T_ARRAY,FIELD_PROJECTOR<T_IF_CLASS,T_FIELD,field> > Project() const
    {return PROJECTED_ARRAY<const T_ARRAY,FIELD_PROJECTOR<ELEMENT,T_FIELD,field> >(Derived());}

    PROJECTED_ARRAY<T_ARRAY,INDEX_PROJECTOR> Project(const ID index)
    {return PROJECTED_ARRAY<T_ARRAY,INDEX_PROJECTOR>(Derived(),INDEX_PROJECTOR(index));}

    PROJECTED_ARRAY<const T_ARRAY,INDEX_PROJECTOR> Project(const ID index) const
    {return PROJECTED_ARRAY<const T_ARRAY,INDEX_PROJECTOR>(Derived(),INDEX_PROJECTOR(index));}

    T_REF_IF_POSSIBLE operator()(const ID i)
    {return Derived()(i);}

    CONST_T_REF_IF_POSSIBLE operator()(const ID i) const
    {return Derived()(i);}

    bool Valid_Index(const ID i) const
    {return (unsigned)Value(i)<Value(Size());}

    ARRAY_VIEW<typename ELEMENT_OF_VECTOR<T>::TYPE> Flattened() // valid only for contiguous arrays of VECTOR<T,d>
    {T_ARRAY& self=Derived();return ARRAY_VIEW<typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    ARRAY_VIEW<const typename ELEMENT_OF_VECTOR<T>::TYPE> Flattened() const // valid only for contiguous arrays of VECTOR<T,d>
    {const T_ARRAY& self=Derived();return ARRAY_VIEW<const typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    template<class ID2>
    ARRAY_VIEW<const T,ID2> Array_View(const ID first,const ID2 length) const
    {const T_ARRAY& self=Derived();assert((unsigned)(Value(first)+Value(length))<=(unsigned)Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<const T,ID2>(length,self.Get_Array_Pointer()+Value(first));}

    template<class ID2>
    ARRAY_VIEW<T,ID2> Array_View(const ID first,const ID2 length)
    {T_ARRAY& self=Derived();assert((unsigned)(Value(first)+Value(length))<=(unsigned)Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<T,ID2>(length,(T*)self.Get_Array_Pointer()+Value(first));}

    ARRAY_VIEW<const T,ID> Array_View(INTERVAL<ID> I) const
    {return Array_View(I.min_corner,I.max_corner-I.min_corner);}

    ARRAY_VIEW<T,ID> Array_View(INTERVAL<ID> I)
    {return Array_View(I.min_corner,I.max_corner-I.min_corner);}

    template<class T_ARRAY1>
    bool operator==(const T_ARRAY1& v) const
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    const T_ARRAY& self=Derived();ID m=self.Size();
    if(m!=v.Size()) return false;for(ID i(0);i<m;i++) if(self(i)!=v(i)) return false;return true;}

    template<class T_ARRAY1>
    bool operator!=(const T_ARRAY1& v) const
    {return !(*this==v);}

    template<class T_ARRAY1>
    T_ARRAY& operator+=(const ARRAY_BASE<T,T_ARRAY1,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY1& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++) self(i)+=v_(i);return self;}

    T_ARRAY& operator+=(const T& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)+=a;return self;}

    template<class T_ARRAY1>
    T_ARRAY& operator-=(const ARRAY_BASE<T,T_ARRAY1,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY1& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++) self(i)-=v_(i);return self;}

    T_ARRAY& operator-=(const T& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)-=a;return self;}

    template<class T2,class T_ARRAY_T2>
    T_ARRAY& operator*=(const ARRAY_BASE<T2,T_ARRAY_T2,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY_T2& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++) self(i)*=v_(i);return self;}

    T_ARRAY& operator*=(const SCALAR& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)*=a;return self;}

    template<class T2,class T_ARRAY_T2>
    T_ARRAY& operator/=(const ARRAY_BASE<T2,T_ARRAY_T2,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY_T2& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++){assert(v_(i));self(i)/=v_(i);}return self;}

    T_ARRAY& operator/=(const SCALAR& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)/=a;return self;}

    T& Last()
    {T_ARRAY& self=Derived();return self(self.Size()-1);}

    const T& Last() const
    {const T_ARRAY& self=Derived();return self(self.Size()-1);}

    ID Size() const
    {return Derived().Size();}

    template<class T_ARRAY2> SCALAR
    Inner_Product(const ARRAY_BASE<SCALAR,T_ARRAY2,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());return (m*(*this*a2)).Sum();}

    template<class T2,class T_ARRAY2> typename DISABLE_IF<IS_SCALAR<T2>::value,SCALAR>::TYPE
    Inner_Product(const ARRAY_BASE<T2,T_ARRAY2,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());typename T_ARRAY2::SCALAR result(0);ID size=Size();for(ID i(0);i<size;i++) result+=m(i).Inner_Product((*this)(i),a2(i));return result;}

    template<class T_ARRAY2> double
    Inner_Product_Double_Precision(const ARRAY_BASE<SCALAR,T_ARRAY2,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());double d=0;for(ID i(0);i<Size();i++) d+=m(i)*(*this)(i).Dot(a2(i));return d;}

    template<class T2,class T_ARRAY2> typename DISABLE_IF<IS_SCALAR<T2>::value,double>::TYPE
    Inner_Product_Double_Precision(const ARRAY_BASE<T2,T_ARRAY2,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());double result(0);ID size=Size();for(ID i(0);i<size;i++) result+=m(i).Inner_Product((*this)(i),a2(i));return result;}

    T Max() const
    {const T_ARRAY& self=Derived();T result=self(ID(0));ID m=self.Size();for(ID i(1);i<m;i++) result=PhysBAM::max(result,self(i));return result;}

    T Max_Abs() const
    {const T_ARRAY& self=Derived();T result=T();;ID m=self.Size();for(ID i(0);i<m;i++) result=PhysBAM::max(result,abs(self(i)));return result;}

    T Max_Mag() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result=PhysBAM::maxmag(result,self(i));return result;}

    ID Arg_Max() const
    {const T_ARRAY& self=Derived();ID result(0),m=self.Size();for(ID i(1);i<m;i++) if(self(i)>self(result)) result=i;return result;}

    T Min() const
    {const T_ARRAY& self=Derived();T result=self(ID(0));ID m=self.Size();for(ID i(1);i<m;i++) result=PhysBAM::min(result,self(i));return result;}

    T Min_Mag() const
    {const T_ARRAY& self=Derived();T result=self(ID(0));ID m=self.Size();for(ID i(1);i<m;i++) result=PhysBAM::minmag(result,self(i));return result;}

    ID Arg_Min() const
    {const T_ARRAY& self=Derived();ID result(0),m=self.Size();for(ID i(1);i<m;i++) if(self(i)<self(result)) result=i;return result;}

    T Componentwise_Max_Abs() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result=T::Componentwise_Max(result,abs(self(i)));return result;}

    T Sum() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result+=self(i);return result;}

    double Sum_Double_Precision() const
    {const T_ARRAY& self=Derived();double result=0;ID m=self.Size();for(ID i(0);i<m;i++) result+=self(i);return result;}

    T Sum_Abs() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result+=abs(self(i));return result;}
    
    T Average() const
    {const T_ARRAY& self=Derived();return self.Size()?Sum()/typename ARRAY_BASE<T,T_ARRAY,ID>::SCALAR(self.Size()):T();}

    template<class T_ARRAY2> typename SCALAR_POLICY<T>::TYPE
    Dot(const ARRAY_BASE<T,T_ARRAY2,ID>& a) const
    {assert(Size()==a.Size());
    typename SCALAR_POLICY<T>::TYPE result(0);ID m=Size();for(ID i(0);i<m;i++) result+=PhysBAM::Dot_Product((*this)(i),a(i));return result;}

    template<class T_ARRAY2> double
    Dot_Double_Precision(const ARRAY_BASE<T,T_ARRAY2,ID>& a) const
    {assert(Size()==a.Size());
    double result(0);ID m=Size();for(ID i(0);i<m;i++) result+=PhysBAM::Dot_Product((*this)(i),a(i));return result;}

    template<class T_ARRAY2> static typename SCALAR_POLICY<T>::TYPE
    Dot_Product(ARRAY_BASE& a1,const ARRAY_BASE<T,T_ARRAY2,ID>& a2)
    {return a1.Dot(a2);}

    template<class T_ARRAY2> static double
    Dot_Product_Double_Precision(const ARRAY_BASE& a1,const ARRAY_BASE<T,T_ARRAY2,ID>& a2)
    {return a1.Dot_Double_Precision(a2);}

    typename SCALAR_POLICY<T>::TYPE Magnitude_Squared() const
    {const T_ARRAY& self=Derived();typename SCALAR_POLICY<T>::TYPE result(0);ID m=self.Size();for(ID i(0);i<m;i++) result+=PhysBAM::Magnitude_Squared(self(i));return result;}

    typename SCALAR_POLICY<T>::TYPE Magnitude() const
    {return sqrt(Magnitude_Squared());}

    T Normalize()
    {Static_Assert_Not_Small();T magnitude=Magnitude();if(magnitude) Derived()*=1/magnitude;else{Fill(0);(*this)(0)=(T)1;}return magnitude;}

    typename SCALAR_POLICY<T>::TYPE Maximum_Magnitude() const
    {return Maximum_Magnitude((T*)0);}

    ID Arg_Maximum_Magnitude() const
    {return Arg_Maximum_Magnitude((ID*)0);}

private:
    template<class U>
    typename ENABLE_IF<IS_SCALAR<T>::value,U>::TYPE Maximum_Magnitude(U*) const
    {T result=(T)0;for(int i=0;i<Size();i++) result=PhysBAM::max(result,abs((*this)(i)));return result;}

    template<class U>
    typename U::SCALAR Maximum_Magnitude(U*) const
    {typename T::SCALAR result(0);for(int i=0;i<Size();i++) result=PhysBAM::max(result,PhysBAM::Magnitude_Squared((*this)(i)));return sqrt(result);}

    template<class U>
    typename ENABLE_IF<IS_SCALAR<T>::value,U>::TYPE Arg_Maximum_Magnitude(U*) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    T maximum=-1;ID argmax=ID();
    for(ID i(0);i<m;i++){T current=abs(self(i));if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}

    template<class U>
    typename DISABLE_IF<IS_SCALAR<T>::value,U>::TYPE Arg_Maximum_Magnitude(U*) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    typename T::SCALAR maximum=-1;ID argmax=ID();
    for(ID i(0);i<m;i++){
        typename T::SCALAR current=self(i).Magnitude_Squared();
        if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}
public:

    template<class T_ARRAY1>
    ELEMENT Weighted_Sum(const T_ARRAY1& weights) const
    {STATIC_ASSERT_SAME(typename T_ARRAY1::ELEMENT,SCALAR);assert(weights.Size()==Size());
    ELEMENT result((ELEMENT()));INDEX m=Size();for(INDEX i(0);i<m;i++) result+=weights(i)*(*this)(i);return result;}

    T_ARRAY Householder_Vector(const int k) const
    {T_ARRAY v((INITIAL_SIZE)Size());T v_dot_v=0;for(int i=k;i<Size();i++){v(i)=(*this)(i);v_dot_v+=sqr(v(i));}
    if((*this)(k)>=0) v(k)+=sqrt(v_dot_v);else v(k)-=sqrt(v_dot_v);
    return v;}

    template<class T_ARRAY2>
    ARRAY_DIFFERENCE<T_ARRAY,ARRAY_LEFT_MULTIPLE<T,T_ARRAY2> > Householder_Transform(const ARRAY_BASE<T,T_ARRAY2>& v) const
    {Assert_Same_Size(*this,v);
    T v_dot_a=0,v_dot_v=0;for(int i=0;i<Size();i++){v_dot_a+=v(i)*(*this)(i);v_dot_v+=sqr(v(i));}
    return *this-2*v_dot_a/v_dot_v*v;}

    void Givens_Rotate(const int i,const int j,const T c,const T s)
    {if(i>j){Givens_Rotate(j,i,c,-s);return;}assert(0<=i && i<j && j<Size());T u=(*this)(i),v=(*this)(j);(*this)(i)=c*u-s*v;(*this)(j)=s*u+c*v;}

    void Set_To_Orthogonal_Vector() // result isn't normalized
    {assert(Size()>=2);int m1=0,m2=1;
    if(abs((*this)(m1))>abs((*this)(m2))) exchange(m1,m2);
    for(int i=2;i<Size();i++) if(abs((*this)(i))>abs((*this)(m1))){
        m1=i;if(abs((*this)(m1))>abs((*this)(m2))) exchange(m1,m2);}
    T x1=(*this)(m1),x2=(*this)(m2);
    Fill(0);
    (*this)(m2)=x1;(*this)(m1)=-x2;}

    template<class T_ARRAY1,class T_ARRAY2>
    static T Angle_Between(const ARRAY_BASE<T,T_ARRAY1>& u,const ARRAY_BASE<T,T_ARRAY2>& v) // 0 .. pi
    {Assert_Same_Size(u,v);T u2=0,u1=u(0),v1=v(0),uv=0;for(int i=1;i<u.Size();i++){T ui=u(i);u2+=sqr(ui);uv+=ui*v(i);}u1+=sign_nonzero(u1)*sqrt(u2+sqr(u1));
    T factor=2*(uv+u1*v1)/(u2+sqr(u1)),R12=v1-factor*u1,R22=0;for(int i=1;i<u.Size();i++) R22+=sqr(v(i)-factor*u(i));return atan2(sqrt(R22),-R12*sign_nonzero(u1));}

    template<class T_ARRAY1,class T_ARRAY2>
    static typename T_ARRAY1::template REBIND<bool>::TYPE Componentwise_Greater_Equal(const ARRAY_BASE<T,T_ARRAY1>& u,const ARRAY_BASE<T,T_ARRAY2>& v)
    {Assert_Same_Size(u,v);typename T_ARRAY1::template REBIND<bool>::TYPE result(INITIAL_SIZE(u.Size()));for(int i=0;i<u.Size();i++) result(i)=u(i)>=v(i);return result;}

    template<class T_ARRAY1,class T_ARRAY2>
    static T_ARRAY1 Componentwise_And(const ARRAY_BASE<bool,T_ARRAY1>& u,const ARRAY_BASE<bool,T_ARRAY2>& v)
    {Assert_Same_Size(u,v);T_ARRAY1 result(INITIAL_SIZE(u.Size()));for(int i=0;i<u.Size();i++) result(i)=(u(i) && v(i));return result;}

    void Reverse()
    {for(ID i(0),s=Size();i<s-1-i;i++) exchange((*this)(i),(*this)(s-1-i));}

    ID Find(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(0);i<m;i++) if(self(i)==element) return i;return -1;}

    bool Find(const T& element,ID& index) const // returns the first occurence of an element in an array
    {return Find(element,0,index);}

    bool Find(const T& element,const ID start_index,ID& index) const // returns the first occurence after start_index of an element in an array
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i=start_index;i<m;i++) if(self(i)==element){index=i;return true;}return false;}

    bool Contains(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(0);i<m;i++) if(self(i)==element) return true;return false;}

    bool Contains_Only(const T& element) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    for(ID i(0);i<m;i++) if(self(i)!=element) return false;return true;}

    int Count_Matches(const T& value) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    int count=0;for(ID i(0);i<m;i++) if(self(i)==value) count++;return count;}

    int Number_True() const
    {STATIC_ASSERT_SAME(T,bool);return Count_Matches(true);}

    int Number_False() const
    {STATIC_ASSERT_SAME(T,bool);return Count_Matches(false);}

    void Get_Unique(ARRAY<T>& array) const
    {const T_ARRAY& self=Derived();HASHTABLE<T> hash(Value(self.Size())*3/2);array.Remove_All();for(int i=0;i<self.Size();i++) if(hash.Set(self(i))) array.Append(self(i));}

    void Prune_Duplicates()
    {T_ARRAY& self=Derived();HASHTABLE<T> hash(Value(self.Size())*3/2);int j=0;for(int i=0;i<self.Size();i++) if(hash.Set(self(i))) self(j++)=self(i);self.Resize(j);}

    void Coalesce()
    {Sort();T_ARRAY& self=Derived();int j=-1;if(self.Size()>0) j=0;for(int i=1;i<self.Size();i++){if(!(self(j)<self(i))) self(j).Merge(self(i));else self(++j)=self(i);}self.Resize(j+1);}

    template<class T_COMPARE>
    void Sort(const T_COMPARE comparison)
    {std::sort(begin(),end(),comparison);}

    template<class T_COMPARE>
    void Stable_Sort(const T_COMPARE& comparison)
    {std::stable_sort(begin(),end(),comparison);}

    void Sort()
    {Sort(std::less<typename T_ARRAY::ELEMENT>());}

    void Stable_Sort()
    {Stable_Sort(std::less<typename T_ARRAY::ELEMENT>());}

    void Fill(T value)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)=value;}

    template<class T_ARRAY1>
    void Copy(const T_ARRAY1& old_copy)
    {*this=old_copy;}

    template<class T_ARRAY1>
    void Copy(const SCALAR constant,const T_ARRAY1& array)
    {*this=constant*array;}

    template<class T_ARRAY1,class T_ARRAY2>
    void Copy(const SCALAR c1,const T_ARRAY1& v1,const T_ARRAY2& v2)
    {*this=c1*v1+v2;}

    template<class T_ARRAY1,class T_ARRAY2>
    void Copy(const SCALAR c1,const T_ARRAY1& v1,const SCALAR c2,const T_ARRAY2& v2)
    {*this=c1*v1+c2*v2;}

    template<class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    void Copy(const SCALAR c1,const T_ARRAY1& v1,const SCALAR c2,const T_ARRAY2& v2,const SCALAR c3,const T_ARRAY3& v3)
    {*this=c1*v1+c2*v2+c3*v3;}

    static void Get(T_ARRAY& new_copy,const T_ARRAY& old_copy)
    {if(&old_copy!=&new_copy) new_copy=old_copy.Prefix(new_copy.Size());}

    static void Put(const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {if(&old_copy!=&new_copy) new_copy.Prefix(old_copy.Size())=old_copy;}

    template<class T2>
    static void Put(const T2 constant,const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {new_copy.Prefix(old_copy.Size())=constant*old_copy;}

    template<class T1,class T_ARRAY1>
    void Copy_With_Offset(const ARRAY_BASE<T1,T_ARRAY1,ID>& old_copy,const ID offset)
    {STATIC_ASSERT(CAN_ASSIGN<T1,T>::value);
    ID m=old_copy.Size();assert(m+offset<=Size());
    for(ID i(0);i<m;i++) (*this)(i+offset)=old_copy(i);}

    void Clamp_Below(const T& value)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)=clamp_min(self(i),value);}

    static void Find_Common_Elements(const T_ARRAY& a,const T_ARRAY& b,T_ARRAY& result)
    {assert(&a!=&result);assert(&b!=&result);result.Remove_All();
    ID m=a.Size();for(ID i(0);i<m;i++) if(b.Contains(a(i))) result.Append(a(i));}

    template<class T_ARRAY1,class T_ARRAY2>
    static bool Equal_Dimensions(const T_ARRAY1& a,const T_ARRAY2& b)
    {return a.Size()==b.Size();}

    template<class T_ARRAY1>
    void Remove_Sorted_Indices(const T_ARRAY1& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY1::ELEMENT);
    T_ARRAY& self=Derived();ID m=self.Size(),index_m=index.Size();
    if(index_m==0) return;
    for(ID kk(0);kk<index_m-1;kk++){
        assert((unsigned)index(kk)<(unsigned)m);
        for(ID i=index(kk)+1-kk;i<=index(kk+1)-1-kk;i++) self(i)=self(i+kk);}
    for(ID i=index(index_m)+1-index_m;i<=m-index_m;i++) self(i)=self(i+index_m);
    self.Resize(m-index_m);}

    template<class T_ARRAY1>
    void Remove_Sorted_Indices_Lazy(const T_ARRAY1& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY1::ELEMENT);
    T_ARRAY& self=Derived();ID index_m=index.Size();
    if(index_m==0) return;
    ID curr=0;
    for(ID k=index_m-1;k>=ID(0);k--)if(index(k)!=curr){curr=index(k);self.Remove_Index_Lazy(curr);}
    self.Compact();}

    int Pack_Size() const
    {return sizeof(T);}

    template<class T_ARRAY1>
    void Pack(T_ARRAY1& buffer,typename T_ARRAY1::INDEX& position,const ID p)
    {T_ARRAY& self=Derived();*(T*)(&buffer(position+1))=self(p);position+=sizeof(T);}

    template<class T_ARRAY1>
    void Unpack(const T_ARRAY1& buffer,typename T_ARRAY1::INDEX& position,const ID p)
    {T_ARRAY& self=Derived();self(p)=*(const T*)(&buffer(position+1));position+=sizeof(T);}

    template<class T_ARRAY1,class T_ARRAY2>
    void Extract(ARRAY_BASE<T,T_ARRAY1,ID>& a,ARRAY_BASE<T,T_ARRAY2,ID>& b) const
    {assert(a.Size()+b.Size()==Size());
    for(ID i(0);i<a.Size();i++) a(i)=(*this)(i);
    for(ID i(0);i<b.Size();i++) b(i)=(*this)(i+a.Size());}

    template<class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    void Extract(ARRAY_BASE<T,T_ARRAY1,ID>& a,ARRAY_BASE<T,T_ARRAY2,ID>& b,ARRAY_BASE<T,T_ARRAY3,ID>& c) const
    {assert(a.Size()+b.Size()+c.Size()==Size());
    for(ID i(0);i<a.Size();i++) a(i)=(*this)(i);
    for(ID i(0);i<b.Size();i++) b(i)=(*this)(i+a.Size());
    for(ID i(0);i<c.Size();i++) c(i)=(*this)(i+a.Size()+b.Size());}

    template<class T_ARRAY1,class T_ARRAY2>
    void Combine(const ARRAY_BASE<T,T_ARRAY1,ID>& a,const ARRAY_BASE<T,T_ARRAY2,ID>& b)
    {assert(a.Size()+b.Size()==Size());
    for(ID i(0);i<a.Size();i++) (*this)(i)=a(i);
    for(ID i(0);i<b.Size();i++) (*this)(i+a.Size())=b(i);}

    template<class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    void Combine(const ARRAY_BASE<T,T_ARRAY1,ID>& a,const ARRAY_BASE<T,T_ARRAY2,ID>& b,const ARRAY_BASE<T,T_ARRAY3,ID>& c)
    {assert(a.Size()+b.Size()+c.Size()==Size());
    for(ID i(0);i<a.Size();i++) (*this)(i)=a(i);
    for(ID i(0);i<b.Size();i++) (*this)(i+a.Size())=b(i);
    for(ID i(0);i<c.Size();i++) (*this)(i+a.Size()+b.Size())=b(i);}

    T* begin() // for stl
    {return Derived().begin();}

    const T* begin() const // for stl
    {return Derived().begin();}

    T* end() // for stl
    {return Derived().end();}

    const T* end() const // for stl
    {return Derived().end();}

    ID Binary_Search(const T& value) const// lower_bound binary search
    {return ID(std::lower_bound(begin(),end(),value)-begin());}

    void Write_Raw(std::ostream& output) const
    {const T_ARRAY& a=Derived();ID m=a.Size();for(ID i(0);i<m;i++){output<<a(i);if(i<m-1) output<<" ";}}

//#####################################################################
};
template<class T,class T_ARRAY,class ID>
inline std::ostream& operator<<(std::ostream& output,const ARRAY_BASE<T,T_ARRAY,ID>& a)
{output<<"(";a.Write_Raw(output);output<<")";return output;}
//#####################################################################
template<class T_ARRAY1,class T_ARRAY2> struct CAN_ASSIGN<T_ARRAY1,T_ARRAY2,typename ENABLE_IF<IS_ARRAY<T_ARRAY1>::value && IS_ARRAY<T_ARRAY2>::value && IS_SAME<typename T_ARRAY1::ELEMENT,typename T_ARRAY2::ELEMENT>::value && !IS_SAME<T_ARRAY1,T_ARRAY2>::value>::TYPE>
{static const bool value=true;};

template<class T,class T_ARRAY,class ID>
static int Hash_Reduce_Array_Helper(const ARRAY_BASE<T,T_ARRAY,ID>& array)
{
    ID i(0);
    int hash=Value(array.Size())%2==0?missing_element_hash:HASH_REDUCE<T>::H(array(i++));
    for(;i<array.Size()-1;i+=2)
        hash=int_hash(hash,HASH_REDUCE<T>::H(array(i)),HASH_REDUCE<T>::H(array(i+1)));
    return Value(array.Size())==1?int_hash(hash):hash;
};

template<class T_ARRAY> struct HASH_REDUCE<T_ARRAY,typename ENABLE_IF<(IS_ARRAY<T_ARRAY>::value && !FIXED_SIZE_VECTOR<T_ARRAY>::value)>::TYPE>
{static int H(const T_ARRAY& key){return Hash_Reduce_Array_Helper(key);}};
}
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_DIFFERENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_LEFT_MULTIPLE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_NEGATION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_PLUS_SCALAR.h>
#include <PhysBAM_Tools/Arrays/ARRAY_PRODUCT.h>
#include <PhysBAM_Tools/Arrays/ARRAY_QUOTIENT.h>
#include <PhysBAM_Tools/Arrays/ARRAY_SUM.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#endif
