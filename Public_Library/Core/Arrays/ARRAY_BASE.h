//#####################################################################
// Copyright 2004-2009, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_BASE
//#####################################################################
#ifndef __ARRAY_BASE__
#define __ARRAY_BASE__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Arrays/SORT.h>
#include <Core/Data_Structures/ELEMENT_ID.h>
#include <Core/Data_Structures/HASH_REDUCE.h>
#include <Core/Math_Tools/exchange.h>
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/maxabs.h>
#include <Core/Math_Tools/maxmag.h>
#include <Core/Math_Tools/min.h>
#include <Core/Math_Tools/minmag.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Utilities/STATIC_ASSERT.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Vectors/ARITHMETIC_POLICY.h>
#include <Core/Vectors/Dot_Product.h>
#include <Core/Vectors/SCALAR_POLICY.h>
#include <algorithm>
#include <iostream>
namespace PhysBAM{

template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

template<class T_ARRAY,class ENABLER=void> struct CANONICALIZE_CONST_ARRAY:public FIRST<T_ARRAY>{};

template<class T_ARRAY0,class T_ARRAY1> struct SAME_ARRAY_CANONICAL{static bool Same_Array(const T_ARRAY0& array0,const T_ARRAY1& array1)
{STATIC_ASSERT(!is_same<T_ARRAY0,T_ARRAY1>::value);return false;}};

template<class T_ARRAY> struct SAME_ARRAY_CANONICAL<T_ARRAY,T_ARRAY>{static bool Same_Array(const T_ARRAY& array0,const T_ARRAY& array1)
{return T_ARRAY::Same_Array(array0,array1);}};

template<class TA0,class TA1> struct SAME_ARRAY:public SAME_ARRAY_CANONICAL<typename CANONICALIZE_CONST_ARRAY<TA0>::TYPE,typename CANONICALIZE_CONST_ARRAY<TA1>::TYPE>{};

template<class TV> struct ELEMENT_OF_VECTOR {private:struct UNUSABLE;public:typedef UNUSABLE TYPE;};
template<class T,int d> struct ELEMENT_OF_VECTOR<VECTOR<T,d> > {typedef T TYPE;};

template<class T> struct CAN_REFERENCE_ELEMENTS {static const int value=true;};
template<class OP,class ID> struct CAN_REFERENCE_ELEMENTS<ARRAY_EXPRESSION<OP,ID> > {static const int value=false;};

template<class T_ARRAY0,class enabler=void> struct EQUIVALENT_ARRAY;
template<class T> struct EQUIVALENT_ARRAY<T,typename enable_if<!IS_ARRAY<T>::value>::type> {typedef T TYPE;};
template<class T> struct EQUIVALENT_ARRAY<T,typename enable_if<IS_ARRAY<T>::value>::type> {typedef ARRAY<typename EQUIVALENT_ARRAY<typename T::ELEMENT>::TYPE> TYPE;};

template<class T,class T_ARRAY,class ID>
class ARRAY_BASE
{
    struct UNUSABLE{};
public:
    typedef T ELEMENT;
    typedef ID INDEX;
    typedef T value_type; // for stl
    typedef int difference_type; // for stl

    typedef typename conditional<CAN_REFERENCE_ELEMENTS<T_ARRAY>::value,T&,T>::type T_REF_IF_POSSIBLE;
    typedef typename conditional<CAN_REFERENCE_ELEMENTS<T_ARRAY>::value,const T&,const T>::type CONST_T_REF_IF_POSSIBLE;
    typedef T& RESULT_TYPE;
    typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;

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
    template<class T_ARRAY1> struct IS_ARRAY_BASE {static const bool value=false;};
    template<class T2,class T_ARRAY1> struct IS_ARRAY_BASE<ARRAY_BASE<T2,T_ARRAY1,ID> > {static const bool value=true;};
public:

    template<class T_ARRAY0> static bool Same_Array(const T_ARRAY0& array0,const T_ARRAY0& array1)
    {return &array0==&array1;}

    template<class T_ARRAY0,class T_ARRAY1> static bool
    Same_Array(const T_ARRAY0& array0,const T_ARRAY1& array1)
    {return SAME_ARRAY<T_ARRAY0,T_ARRAY1>::Same_Array(array0,array1);}

protected:
    T_ARRAY& operator=(const ARRAY_BASE& source)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY& source_=source.Derived();assert(m==source_.Size());
    if(!T_ARRAY::Same_Array(self,source_)) for(ID i(0);i<m;i++) self(i)=source_(i);
    return self;}

    template<class T_ARRAY0>
    T_ARRAY& operator=(const T_ARRAY0& source)
    {STATIC_ASSERT(CAN_ASSIGN<T,typename T_ARRAY0::ELEMENT>::value);
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
    typedef typename conditional<is_class<T>::value,T,UNUSABLE>::type T_IF_CLASS;
public:

    template<int d>
    static void Static_Assert_Not_Small(const VECTOR<T,d>*)
    {STATIC_ASSERT(d>3);}

    template<class T_ARRAY0>
    static void Static_Assert_Not_Small(const T_ARRAY0*)
    {}

    void Static_Assert_Not_Small() const
    {Static_Assert_Not_Small((T_ARRAY*)0);}

    template<class T2,class T3,int p,int q>
    static void Assert_Same_Size_Helper(const VECTOR<T2,p>&,const VECTOR<T3,q>&)
    {STATIC_ASSERT(p==q);}

    template<class T2,class T3,class T_ARRAY1,class T_ARRAY2>
    static void Assert_Same_Size_Helper(const ARRAY_BASE<T2,T_ARRAY1>& u,const ARRAY_BASE<T3,T_ARRAY2>& v)
    {assert(u.Size()==v.Size());}

    template<class T2,class T3,class T_ARRAY1,class T_ARRAY2>
    static void Assert_Same_Size(const ARRAY_BASE<T2,T_ARRAY1>& u,const ARRAY_BASE<T3,T_ARRAY2>& v)
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

    template<class ID1>
    ARRAY_VIEW<const T,ID1> Array_View(const ID first,const ID1 length) const
    {const T_ARRAY& self=Derived();assert((unsigned)(Value(first)+Value(length))<=(unsigned)Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<const T,ID1>(length,self.Get_Array_Pointer()+Value(first));}

    template<class ID1>
    ARRAY_VIEW<T,ID1> Array_View(const ID first,const ID1 length)
    {T_ARRAY& self=Derived();assert((unsigned)(Value(first)+Value(length))<=(unsigned)Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<T,ID1>(length,(T*)self.Get_Array_Pointer()+Value(first));}

    ARRAY_VIEW<const T,ID> Array_View(INTERVAL<ID> I) const
    {return Array_View(I.min_corner,I.max_corner-I.min_corner);}

    ARRAY_VIEW<T,ID> Array_View(INTERVAL<ID> I)
    {return Array_View(I.min_corner,I.max_corner-I.min_corner);}

    template<class T_ARRAY0>
    bool operator==(const T_ARRAY0& v) const
    {STATIC_ASSERT_SAME(T,typename T_ARRAY0::ELEMENT);
    const T_ARRAY& self=Derived();ID m=self.Size();
    if(m!=v.Size()) return false;for(ID i(0);i<m;i++) if(self(i)!=v(i)) return false;return true;}

    template<class T_ARRAY0>
    bool operator!=(const T_ARRAY0& v) const
    {return !(*this==v);}

    template<class T_ARRAY0>
    T_ARRAY& operator+=(const ARRAY_BASE<T,T_ARRAY0,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY0& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++) self(i)+=v_(i);return self;}

    T_ARRAY& operator+=(const T& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)+=a;return self;}

    const T_ARRAY& operator+() const
    {return *this;};

    template<class T_ARRAY0>
    T_ARRAY& operator-=(const ARRAY_BASE<T,T_ARRAY0,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY0& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++) self(i)-=v_(i);return self;}

    T_ARRAY& operator-=(const T& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)-=a;return self;}

    template<class T2,class T_ARRAY_T1>
    T_ARRAY& operator*=(const ARRAY_BASE<T2,T_ARRAY_T1,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY_T1& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++) self(i)*=v_(i);return self;}

    T_ARRAY& operator*=(const SCALAR& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)*=a;return self;}

    template<class T2,class T_ARRAY_T1>
    T_ARRAY& operator/=(const ARRAY_BASE<T2,T_ARRAY_T1,ID>& v)
    {T_ARRAY& self=Derived();ID m=self.Size();const T_ARRAY_T1& v_=v.Derived();assert(m==v_.Size());for(ID i(0);i<m;i++){assert(v_(i));self(i)/=v_(i);}return self;}

    T_ARRAY& operator/=(const SCALAR& a)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)/=a;return self;}

    T& Last()
    {T_ARRAY& self=Derived();return self(self.Size()-1);}

    const T& Last() const
    {const T_ARRAY& self=Derived();return self(self.Size()-1);}

    ID Size() const
    {return Derived().Size();}

    template<class T_ARRAY1> SCALAR
    Inner_Product(const ARRAY_BASE<SCALAR,T_ARRAY1,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());return (m*(*this*a2)).Sum();}

    template<class T2,class T_ARRAY1> typename enable_if<!is_scalar<T2>::value,SCALAR>::type
    Inner_Product(const ARRAY_BASE<T2,T_ARRAY1,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());typename T_ARRAY1::SCALAR result(0);ID size=Size();for(ID i(0);i<size;i++) result+=m(i).Inner_Product((*this)(i),a2(i));return result;}

    template<class T_ARRAY1> double
    Inner_Product_Double_Precision(const ARRAY_BASE<SCALAR,T_ARRAY1,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());double d=0;for(ID i(0);i<Size();i++) d+=m(i)*(*this)(i).Dot(a2(i));return d;}

    template<class T2,class T_ARRAY1> typename enable_if<!is_scalar<T2>::value,double>::type
    Inner_Product_Double_Precision(const ARRAY_BASE<T2,T_ARRAY1,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
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

    T Product() const
    {const T_ARRAY& self=Derived();T result(1);ID m=self.Size();for(ID i(0);i<m;i++) result*=self(i);return result;}

    double Sum_Double_Precision() const
    {const T_ARRAY& self=Derived();double result=0;ID m=self.Size();for(ID i(0);i<m;i++) result+=self(i);return result;}

    T Sum_Abs() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result+=abs(self(i));return result;}
    
    T Average() const
    {const T_ARRAY& self=Derived();return self.Size()?Sum()/typename ARRAY_BASE<T,T_ARRAY,ID>::SCALAR(self.Size()):T();}

    T Lp_Norm(const T& p) const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result+=pow(abs(self(i)),p);return pow(result,1/p);}

    template<class T_ARRAY1> typename SCALAR_POLICY<T>::TYPE
    Dot(const ARRAY_BASE<T,T_ARRAY1,ID>& a) const
    {assert(Size()==a.Size());
    typename SCALAR_POLICY<T>::TYPE result(0);ID m=Size();for(ID i(0);i<m;i++) result+=PhysBAM::Dot_Product((*this)(i),a(i));return result;}

    template<class T_ARRAY1> double
    Dot_Double_Precision(const ARRAY_BASE<T,T_ARRAY1,ID>& a) const
    {assert(Size()==a.Size());
    double result(0);ID m=Size();for(ID i(0);i<m;i++) result+=PhysBAM::Dot_Product((*this)(i),a(i));return result;}

    template<class T_ARRAY1> static typename SCALAR_POLICY<T>::TYPE
    Dot_Product(ARRAY_BASE& a1,const ARRAY_BASE<T,T_ARRAY1,ID>& a2)
    {return a1.Dot(a2);}

    template<class T_ARRAY1> static double
    Dot_Product_Double_Precision(const ARRAY_BASE& a1,const ARRAY_BASE<T,T_ARRAY1,ID>& a2)
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
    typename enable_if<is_scalar<T>::value,U>::type Maximum_Magnitude(U*) const
    {T result=(T)0;for(int i=0;i<Size();i++) result=PhysBAM::max(result,abs((*this)(i)));return result;}

    template<class U>
    typename U::SCALAR Maximum_Magnitude(U*) const
    {typename T::SCALAR result(0);for(int i=0;i<Size();i++) result=PhysBAM::max(result,PhysBAM::Magnitude_Squared((*this)(i)));return sqrt(result);}

    template<class U>
    typename enable_if<is_scalar<T>::value,U>::type Arg_Maximum_Magnitude(U*) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    T maximum=-1;ID argmax=ID();
    for(ID i(0);i<m;i++){T current=abs(self(i));if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}

    template<class U>
    typename enable_if<!is_scalar<T>::value,U>::type Arg_Maximum_Magnitude(U*) const
    {const T_ARRAY& self=Derived();ID m=self.Size();
    typename T::SCALAR maximum=-1;ID argmax=ID();
    for(ID i(0);i<m;i++){
        typename T::SCALAR current=self(i).Magnitude_Squared();
        if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}
public:

    template<class T_ARRAY0>
    ELEMENT Weighted_Sum(const T_ARRAY0& weights) const
    {STATIC_ASSERT_SAME(typename T_ARRAY0::ELEMENT,SCALAR);assert(weights.Size()==Size());
    ELEMENT result((ELEMENT()));INDEX m=Size();for(INDEX i(0);i<m;i++) result+=weights(i)*(*this)(i);return result;}

    T_ARRAY Householder_Vector(const int k) const
    {T_ARRAY v((INITIAL_SIZE)Size());T v_dot_v=0;for(int i=k;i<Size();i++){v(i)=(*this)(i);v_dot_v+=sqr(v(i));}
    if((*this)(k)>=0) v(k)+=sqrt(v_dot_v);else v(k)-=sqrt(v_dot_v);
    return v;}

    template<class T_ARRAY1>
    auto Householder_Transform(const ARRAY_BASE<T,T_ARRAY1>& v) const
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
    T x0=(*this)(m1),x1=(*this)(m2);
    Fill(0);
    (*this)(m2)=x0;(*this)(m1)=-x1;}

    template<class T_ARRAY0,class T_ARRAY1>
    static T Angle_Between(const ARRAY_BASE<T,T_ARRAY0>& u,const ARRAY_BASE<T,T_ARRAY1>& v) // 0 .. pi
    {Assert_Same_Size(u,v);T u2=0,u1=u(0),v1=v(0),uv=0;for(int i=1;i<u.Size();i++){T ui=u(i);u2+=sqr(ui);uv+=ui*v(i);}u1+=sign_nonzero(u1)*sqrt(u2+sqr(u1));
    T factor=2*(uv+u1*v1)/(u2+sqr(u1)),R01=v1-factor*u1,R11=0;for(int i=1;i<u.Size();i++) R11+=sqr(v(i)-factor*u(i));return atan2(sqrt(R11),-R01*sign_nonzero(u1));}

    template<class T_ARRAY1>
    auto Componentwise_Greater_Equal(const ARRAY_BASE<T,T_ARRAY1>& v) const
    {Assert_Same_Size(*this,v);return Array_Expression_Helper(*this,v,[](const T& a,const T& b){return a>=b;});}

    template<class T_ARRAY1>
    auto Componentwise_And(const ARRAY_BASE<bool,T_ARRAY1>& v) const
    {Assert_Same_Size(*this,v);STATIC_ASSERT((is_same<bool,T>::value));return Array_Expression_Helper(*this,v,[](bool a,bool b){return a && b;});}

    void Reverse()
    {for(ID i(0),s=Size();i<s-1-i;i++) exchange((*this)(i),(*this)(s-1-i));}

    template<class T_ARRAY0>
    bool All_Greater_Equal(const ARRAY_BASE<T,T_ARRAY0,ID>& v) const
    {Assert_Same_Size(*this,v);const T_ARRAY& self=Derived();const T_ARRAY0& vd=v.Derived();
    for(int i=0,n=self.Size();i<n;i++) if(self(i)<vd(i)) return false;return true;}

    template<class T_ARRAY0>
    bool All_Less_Equal(const ARRAY_BASE<T,T_ARRAY0,ID>& v) const
    {return v.All_Greater_Equal(*this);}

    template<class T_ARRAY0>
    bool All_Greater(const ARRAY_BASE<T,T_ARRAY0,ID>& v) const
    {return !v.All_Greater_Equal(*this);}

    template<class T_ARRAY0>
    bool All_Less(const ARRAY_BASE<T,T_ARRAY0,ID>& v) const
    {return !All_Greater_Equal(v);}

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

    template<class T_ARRAY0>
    void Copy(const T_ARRAY0& old_copy)
    {*this=old_copy;}

    template<class T_ARRAY0>
    void Copy(const SCALAR constant,const T_ARRAY0& array)
    {*this=constant*array;}

    template<class T_ARRAY0,class T_ARRAY1>
    void Copy(const SCALAR c1,const T_ARRAY0& v1,const T_ARRAY1& v2)
    {*this=c1*v1+v2;}

    template<class T_ARRAY0,class T_ARRAY1>
    void Copy(const SCALAR c1,const T_ARRAY0& v1,const SCALAR c2,const T_ARRAY1& v2)
    {*this=c1*v1+c2*v2;}

    template<class T_ARRAY0,class T_ARRAY1,class T_ARRAY2>
    void Copy(const SCALAR c1,const T_ARRAY0& v1,const SCALAR c2,const T_ARRAY1& v2,const SCALAR c3,const T_ARRAY2& v3)
    {*this=c1*v1+c2*v2+c3*v3;}

    static void Get(T_ARRAY& new_copy,const T_ARRAY& old_copy)
    {if(&old_copy!=&new_copy) new_copy=old_copy.Prefix(new_copy.Size());}

    static void Put(const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {if(&old_copy!=&new_copy) new_copy.Prefix(old_copy.Size())=old_copy;}

    template<class T2>
    static void Put(const T2 constant,const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {new_copy.Prefix(old_copy.Size())=constant*old_copy;}

    template<class T0,class T_ARRAY0>
    void Copy_With_Offset(const ARRAY_BASE<T0,T_ARRAY0,ID>& old_copy,const ID offset)
    {STATIC_ASSERT(CAN_ASSIGN<T0,T>::value);
    ID m=old_copy.Size();assert(m+offset<=Size());
    for(ID i(0);i<m;i++) (*this)(i+offset)=old_copy(i);}

    void Clamp_Below(const T& value)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)=clamp_min(self(i),value);}

    template<class T_ARRAY0,class T_ARRAY1>
    static bool Equal_Dimensions(const T_ARRAY0& a,const T_ARRAY1& b)
    {return a.Size()==b.Size();}

    template<class T_ARRAY0>
    void Remove_Sorted_Indices(const T_ARRAY0& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY0::ELEMENT);
    T_ARRAY& self=Derived();ID m=self.Size(),index_m=index.Size();
    if(index_m==0) return;
    for(ID kk(0);kk<index_m-1;kk++){
        assert((unsigned)index(kk)<(unsigned)m);
        for(ID i=index(kk)+1-kk;i<=index(kk+1)-1-kk;i++) self(i)=self(i+kk);}
    for(ID i=index(index_m)+1-index_m;i<=m-index_m;i++) self(i)=self(i+index_m);
    self.Resize(m-index_m);}

    template<class T_ARRAY0>
    void Remove_Sorted_Indices_Lazy(const T_ARRAY0& index)
    {STATIC_ASSERT_SAME(ID,typename T_ARRAY0::ELEMENT);
    T_ARRAY& self=Derived();ID index_m=index.Size();
    if(index_m==0) return;
    ID curr=0;
    for(ID k=index_m-1;k>=ID(0);k--)if(index(k)!=curr){curr=index(k);self.Remove_Index_Lazy(curr);}
    self.Compact();}

    int Pack_Size() const
    {return sizeof(T);}

    template<class T_ARRAY0>
    void Pack(T_ARRAY0& buffer,typename T_ARRAY0::INDEX& position,const ID p)
    {T_ARRAY& self=Derived();*(T*)(&buffer(position+1))=self(p);position+=sizeof(T);}

    template<class T_ARRAY0>
    void Unpack(const T_ARRAY0& buffer,typename T_ARRAY0::INDEX& position,const ID p)
    {T_ARRAY& self=Derived();self(p)=*(const T*)(&buffer(position+1));position+=sizeof(T);}

    template<class T_ARRAY0,class T_ARRAY1>
    void Extract(ARRAY_BASE<T,T_ARRAY0,ID>& a,ARRAY_BASE<T,T_ARRAY1,ID>& b) const
    {assert(a.Size()+b.Size()==Size());
    for(ID i(0);i<a.Size();i++) a(i)=(*this)(i);
    for(ID i(0);i<b.Size();i++) b(i)=(*this)(i+a.Size());}

    template<class T_ARRAY0,class T_ARRAY1,class T_ARRAY2>
    void Extract(ARRAY_BASE<T,T_ARRAY0,ID>& a,ARRAY_BASE<T,T_ARRAY1,ID>& b,ARRAY_BASE<T,T_ARRAY2,ID>& c) const
    {assert(a.Size()+b.Size()+c.Size()==Size());
    for(ID i(0);i<a.Size();i++) a(i)=(*this)(i);
    for(ID i(0);i<b.Size();i++) b(i)=(*this)(i+a.Size());
    for(ID i(0);i<c.Size();i++) c(i)=(*this)(i+a.Size()+b.Size());}

    template<class T_ARRAY0,class T_ARRAY1>
    void Combine(const ARRAY_BASE<T,T_ARRAY0,ID>& a,const ARRAY_BASE<T,T_ARRAY1,ID>& b)
    {assert(a.Size()+b.Size()==Size());
    for(ID i(0);i<a.Size();i++) (*this)(i)=a(i);
    for(ID i(0);i<b.Size();i++) (*this)(i+a.Size())=b(i);}

    template<class T_ARRAY0,class T_ARRAY1,class T_ARRAY2>
    void Combine(const ARRAY_BASE<T,T_ARRAY0,ID>& a,const ARRAY_BASE<T,T_ARRAY1,ID>& b,const ARRAY_BASE<T,T_ARRAY2,ID>& c)
    {assert(a.Size()+b.Size()+c.Size()==Size());
    for(ID i(0);i<a.Size();i++) (*this)(i)=a(i);
    for(ID i(0);i<b.Size();i++) (*this)(i+a.Size())=b(i);
    for(ID i(0);i<c.Size();i++) (*this)(i+a.Size()+b.Size())=b(i);}

    auto begin() // for stl
    {return Derived().begin();}

    auto begin() const // for stl
    {return Derived().begin();}

    auto end() // for stl
    {return Derived().end();}

    auto end() const // for stl
    {return Derived().end();}

    ID Binary_Search(const T& value) const // lower_bound binary search
    {return ID(std::lower_bound(begin(),end(),value)-begin());}

    template<class T_COMPARE>
    ID Binary_Search(const T& value,const T_COMPARE comparison) const // lower_bound binary search
    {return ID(std::lower_bound(begin(),end(),value,comparison)-begin());}

    void Write_Raw(std::ostream& output) const
    {const T_ARRAY& a=Derived();ID m=a.Size();for(ID i(0);i<m;i++){output<<a(i);if(i<m-1) output<<" ";}}

//#####################################################################
};

template<class T,class T_ARRAY,class ID>
inline std::ostream& operator<<(std::ostream& output,const ARRAY_BASE<T,T_ARRAY,ID>& a)
{output<<"(";a.Write_Raw(output);output<<")";return output;}
//#####################################################################
template<class T_ARRAY0,class T_ARRAY1> struct CAN_ASSIGN<T_ARRAY0,T_ARRAY1,typename enable_if<IS_ARRAY<T_ARRAY0>::value && IS_ARRAY<T_ARRAY1>::value && is_same<typename T_ARRAY0::ELEMENT,typename T_ARRAY1::ELEMENT>::value && !is_same<T_ARRAY0,T_ARRAY1>::value>::type>
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

template<class T_ARRAY> struct HASH_REDUCE<T_ARRAY,typename enable_if<(IS_ARRAY<T_ARRAY>::value && !FIXED_SIZE_VECTOR<T_ARRAY>::value)>::type>
{static int H(const T_ARRAY& key){return Hash_Reduce_Array_Helper(key);}};
}
//#####################################################################
#include <Core/Arrays/ARRAY_EXPRESSION.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Arrays/IDENTITY_ARRAY.h>
#include <Core/Arrays/INDIRECT_ARRAY.h>
#endif
