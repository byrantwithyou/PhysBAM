//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR
//#####################################################################
#ifndef __VECTOR__
#define __VECTOR__

#include <Tools/Data_Structures/ELEMENT_ID.h>
#include <Tools/Math_Tools/Inverse.h>
#include <Tools/Utilities/STATIC_ASSERT.h>
#include <Tools/Vectors/VECTOR_0D.h>
#include <Tools/Vectors/VECTOR_1D.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <cmath>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
using ::std::abs;
using ::std::floor;
using ::std::ceil;
using ::std::sqrt;
using ::std::exp;
using ::std::sin;
using ::std::cos;
using ::std::pow;
namespace PhysBAM{

template<class T_ARRAY,class T_INDICES> class INDIRECT_ARRAY;
template<class T,int d> struct IS_ARRAY<VECTOR<T,d> > {static const bool value=true;};

template<class T,int d>
class VECTOR:public ARRAY_BASE<T,VECTOR<T,d> >
{
    STATIC_ASSERT(d>3);
    struct UNUSABLE{};
public:
    typedef ARRAY_BASE<T,VECTOR<T,d> > BASE;
    template<class T2> struct REBIND{typedef VECTOR<T2,d> TYPE;};
    typedef typename IF<IS_SCALAR<T>::value,T,UNUSABLE>::TYPE SCALAR;
    typedef T ELEMENT;
    typedef UNUSABLE SPIN;
    typedef int INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    enum WORKAROUND1 {dimension=d};
    enum WORKAROUND2 {m=d};
    typedef int HAS_UNTYPED_READ_WRITE;
    using BASE::Assert_Same_Size;

    T array[d];

    VECTOR()
    {
        STATIC_ASSERT(sizeof(VECTOR)==d*sizeof(T));
        for(int i=0;i<d;i++) array[i]=T();
    }

    explicit VECTOR(INITIAL_SIZE n)
    {
        assert(n==INITIAL_SIZE(d));
        for(int i=0;i<d;i++) array[i]=T();
    }

    template<class ...Args>
    VECTOR(const T& a,const T& b,const T& c,const T& e,Args &&... args)
    {STATIC_ASSERT(d==4+sizeof...(Args));Set_Helper(a,b,c,e,args...);}

    template<class ...Args>
    void Set_Helper(const T& a,Args&&... args)
    {array[d-1-sizeof...(Args)]=a;Set_Helper(args...);}

    template<class ...Args>
    void Set_Helper()
    {}

    template<class ...Args>
    void Get_Helper(T& a,Args&... args) const
    {a=array[d-1-sizeof...(Args)];Get_Helper(args...);}

    template<class ...Args>
    void Get_Helper() const
    {}

    template<class T_VECTOR>
    explicit VECTOR(const ARRAY_BASE<T,T_VECTOR>& v)
    {
        Assert_Same_Size(*this,v);
        for(int i=0;i<d;i++) array[i]=v(i);
    }

    template<class T2,int d2>
    explicit VECTOR(const VECTOR<T2,d2>& v)
    {
        STATIC_ASSERT(d2<=d);
        for(int i=0;i<d2;i++) array[i]=(T)v[i];
        for(int i=d2;i<d;i++) array[i]=T();
    }

    VECTOR(const VECTOR& v)
        :BASE()
    {
        for(int i=0;i<d;i++) array[i]=v.array[i];
    }

    template<int n>
    VECTOR(const VECTOR<T,n>& v1,const VECTOR<T,d-n>& v2)
    {
        for(int i=0;i<n;i++) (*this)(i)=v1(i);
        for(int i=n;i<d;i++) (*this)(i)=v2(i-n);
    }

    VECTOR& operator=(const VECTOR& v)
    {
        for(int i=0;i<d;i++) array[i]=v(i);
        return *this;
    }

    template<class T_VECTOR>
    VECTOR& operator=(const ARRAY_BASE<T,T_VECTOR>& v)
    {
        Assert_Same_Size(*this,v);
        for(int i=0;i<d;i++) array[i]=v(i);
        return *this;
    }

    int Size() const
    {return m;}

    const T& operator[](const int i) const
    {assert((unsigned)i<d);return array[i];}

    T& operator[](const int i)
    {assert((unsigned)i<d);return array[i];}

    const T& operator()(const int i) const
    {assert((unsigned)i<d);return array[i];}

    T& operator()(const int i)
    {assert((unsigned)i<d);return array[i];}

    bool operator==(const VECTOR& v) const
    {for(int i=0;i<d;i++) if(array[i]!=v.array[i]) return false;return true;}

    bool operator!=(const VECTOR& v) const
    {return !((*this)==v);}

    VECTOR operator-() const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=-array[i];return r;}

    VECTOR operator+(const T& a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]+a;return r;}

    VECTOR operator-(const T& a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]-a;return r;}

    VECTOR operator*(const T& a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]*a;return r;}

    VECTOR operator/(const T& a) const
    {return *this*Inverse(a);}

    VECTOR operator+(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]+v.array[i];return r;}

    VECTOR operator-(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]-v.array[i];return r;}

    VECTOR operator*(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]*v.array[i];return r;}

    VECTOR operator/(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]/v.array[i];return r;}

    VECTOR operator*(const INT_INVERSE a) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=array[i]*a;return r;}

    VECTOR Normalized() const
    {VECTOR r(*this);r.Normalize();return r;}

    bool Elements_Equal() const
    {bool equal=true;for(int i=1;i<d;i++) equal&=(array[i]==array[0]);return equal;}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {return v1.Componentwise_Min(v2);}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {return v1.Componentwise_Max(v2);}

    VECTOR Componentwise_Min(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=min(array[i],v.array[i]);return r;}

    VECTOR Componentwise_Max(const VECTOR& v) const
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=max(array[i],v.array[i]);return r;}

    VECTOR Projected_On_Unit_Direction(const VECTOR& direction) const
    {return Dot_Product(*this,direction)*direction;}

    VECTOR Projected(const VECTOR& direction) const // un-normalized direction
    {return Dot_Product(*this,direction)/direction.Magnitude_Squared()*direction;}

    void Project_On_Unit_Direction(const VECTOR& direction)
    {*this=Dot_Product(*this,direction)*direction;}

    void Project(const VECTOR& direction) // un-normalized direction
    {*this=Dot_Product(*this,direction)/direction.Magnitude_Squared()*direction;}

    VECTOR Projected_Orthogonal_To_Unit_Direction(const VECTOR& direction) const
    {return *this-Dot_Product(*this,direction)*direction;}

    void Project_Orthogonal_To_Unit_Direction(const VECTOR& direction)
    {*this-=Dot_Product(*this,direction)*direction;}

    int Number_True() const
    {STATIC_ASSERT((IS_SAME<T,bool>::value));int count=0;for(int i=0;i<d;i++)if(array[i]) count++;return count;}

    static VECTOR Axis_Vector(const int axis)
    {VECTOR r;r(axis)=(T)1;return r;}

    static VECTOR Constant_Vector(const T& constant)
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=constant;return r;}

    static VECTOR All_Ones_Vector()
    {return Constant_Vector(1);}

    void Fill(const T& constant)
    {for(int i=0;i<d;i++) array[i]=constant;}

    template<class ...Args>
    void Get(Args&... args) const
    {STATIC_ASSERT(d==sizeof...(Args));Get_Helper(args...);}

    template<class ...Args>
    void Set(Args&&... args)
    {STATIC_ASSERT(d==sizeof...(Args));Set_Helper(args...);}

    template<class T_FUNCTION>
    static VECTOR Map(const T_FUNCTION& f,const VECTOR& v)
    {VECTOR r;for(int i=0;i<d;i++) r.array[i]=f(v.array[i]);return r;}

    int Find(const T& element) const
    {for(int i=0;i<d;i++) if(array[i]==element) return i;return -1;}

    bool Contains(const T& element) const
    {for(int i=0;i<d;i++) if(array[i]==element) return true;return false;}

    template<class T_ARRAY>
    bool Contains_All(const T_ARRAY& elements) const
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=0;i<elements.Size();i++) if(!Contains(elements(i))) return false;
    return true;}

    template<class T_ARRAY>
    bool Contains_Any(const T_ARRAY& elements) const
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=0;i<elements.Size();i++) if(Contains(elements(i))) return true;
    return false;}

    VECTOR<T,d-1> Remove_Index(const int index) const
    {assert((unsigned)index<d);VECTOR<T,d-1> r;for(int i=0;i<d-1;i++) r[i]=(*this)[i+(i>=index)];return r;}

    VECTOR<T,d+1> Insert(const T& element,const int index) const
    {VECTOR<T,d+1> r;r[index]=element;for(int i=0;i<d;i++) r[i+(i>=index)]=(*this)[i];return r;}

    VECTOR<T,d+1> Prepend(const T& element) const
    {VECTOR<T,d+1> r;r[0]=element;for(int i=0;i<d;i++) r[i+1]=(*this)[i];return r;}

    VECTOR<T,d+1> Append(const T& element) const
    {VECTOR<T,d+1> r;for(int i=0;i<d;i++) r[i]=(*this)[i];r[d]=element;return r;}

    template<int d2> VECTOR<T,d+d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {VECTOR<T,d+d2> r;
    for(int i=0;i<d;i++) r[i]=(*this)[i];
    for(int i=0;i<d2;i++) r[i+d]=elements[i];
    return r;}

    VECTOR Sorted() const
    {VECTOR r(*this);r.Sort();return r;}

    VECTOR Reversed() const
    {VECTOR r;for(int i=0;i<d;i++) r.array[d-i-1]=array[i];return r;}

    template<int d1,int d2> VECTOR<T,d2-d1+1> Slice() const
    {STATIC_ASSERT(((0<=d1) && (d2<d)));
    VECTOR<T,d2-d1+1> r;for(int i=d1;i<=d2;i++) r[i-d1]=(*this)[i];return r;}

    template<int n> void Split(VECTOR<T,n>& v1,VECTOR<T,d-n>& v2) const
    {for(int i=0;i<n;i++) v1(i)=(*this)(i);
    for(int i=n;i<d;i++) v2(i-n)=(*this)(i);}

    VECTOR Orthogonal_Vector() const
    {int a=Dominant_Axis(),b=!a;VECTOR r;r(a)=(*this)(b);r(b)=-(*this)(a);return r;}

    VECTOR Unit_Orthogonal_Vector() const
    {return Orthogonal_Vector().Normalized();}

    template<class T_VECTOR>
    void Set_Subvector(const int istart,const T_VECTOR& v)
    {for(int i=0;i<v.Size();i++) (*this)(istart+i)=v(i);}

    template<class T_VECTOR>
    void Add_Subvector(const int istart,const T_VECTOR& v)
    {for(int i=0;i<v.Size();i++) (*this)(istart+i)+=v(i);}
    
    template<class T_VECTOR>
    void Get_Subvector(const int istart,T_VECTOR& v) const
    {for(int i=0;i<v.Size();i++) v(i)=(*this)(istart+i);}

    T* Get_Array_Pointer()
    {return array;}

    const T* Get_Array_Pointer() const
    {return array;}

    T* begin() // for stl
    {return array;}

    const T* begin() const // for stl
    {return array;}

    T* end() // for stl
    {return array+d;}

    const T* end() const // for stl
    {return array+d;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary_Array<RW>(input,array,d);}

    template<class RW> void Write(std::ostream& output) const
    {Write_Binary_Array<RW>(output,array,d);}
//#####################################################################
    int Dominant_Axis() const;
//#####################################################################
};

//#####################################################################
// Miscellaneous free operators and functions
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
operator+(const T& a,const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=a+v.array[i];return r;}

template<class T,int d> inline VECTOR<T,d>
operator*(const T& a,const VECTOR<T,d>& v)
{return v*a;}

template<class T,int d> inline VECTOR<T,d>
operator-(const T& a,const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=a-v.array[i];return r;}

template<class T,int d> inline VECTOR<T,d>
operator/(const T& a,const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=a/v.array[i];return r;}

template<class T,int d> inline VECTOR<T,d>
Inverse(const VECTOR<T,d>& v)
{VECTOR<T,d> r;for(int i=0;i<d;i++) r.array[i]=1/v.array[i];return r;}

//#####################################################################
// Function Dominant_Axis
//#####################################################################
template<class T,int d> inline int VECTOR<T,d>::
Dominant_Axis() const
{
    int axis=0;
    T abs_max=abs(array[0]);
    for(int i=1;i<d;i++){T abs_i=abs(array[i]);if(abs_max<abs_i){abs_max=abs_i;axis=i;}}
    return axis;
}
//#####################################################################
// Function abs
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
abs(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=abs(v.array[i]);
    return r;
}
//#####################################################################
// Function floor
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
floor(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=floor(v.array[i]);
    return r;
}
//#####################################################################
// Function ceil
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
ceil(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=ceil(v.array[i]);
    return r;
}
//#####################################################################
// Function rint
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
rint(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=rint(v.array[i]);
    return r;
}
//#####################################################################
// Function exp
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
exp(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=exp(v.array[i]);
    return r;
}
//#####################################################################
// Function log
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
log(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=log(v.array[i]);
    return r;
}
//#####################################################################
// Function sin
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
sin(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=sin(v.array[i]);
    return r;
}
//#####################################################################
// Function cos
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
cos(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=cos(v.array[i]);
    return r;
}
//#####################################################################
// Function sqrt
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
sqrt(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=sqrt(v.array[i]);
    return r;
}
//#####################################################################
// Function clamp
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
clamp(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin,const VECTOR<T,d>& vmax)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp(v.array[i],vmin.array[i],vmax.array[i]);
    return r;
}
//#####################################################################
// Function clamp
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
clamp(const VECTOR<T,d>& v,const T& min,const T& max)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp(v.array[i],min,max);
    return r;
}
//#####################################################################
// Function clamp_min
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
clamp_min(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_min(v.array[i],vmin.array[i]);
    return r;
}
//#####################################################################
// Function clamp_min
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
clamp_min(const VECTOR<T,d>& v,const T& min)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_min(v.array[i],min);
    return r;
}
//#####################################################################
// Function clamp_max
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
clamp_max(const VECTOR<T,d>& v,const VECTOR<T,d>& vmax)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_max(v.array[i],vmax.array[i]);
    return r;
}
//#####################################################################
// Function clamp_max
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
clamp_max(const VECTOR<T,d>& v,const T& max)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_max(v.array[i],max);
    return r;
}
//#####################################################################
// Function in_bounds
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
in_bounds(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin,const VECTOR<T,d>& vmax)
{
    for(int i=0;i<d;i++) if(!in_bounds(v.array[i],vmin.array[i],vmax.array[i])) return false;
    return true;
}

//#####################################################################
// Function wrap
//#####################################################################
template<class T,int d> inline VECTOR<T,d>
wrap(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin,const VECTOR<T,d>& vmax)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r(i)=wrap(v(i),vmin(i),vmax(i));
    return r;
}
//#####################################################################

//#####################################################################
// Vector construction
//#####################################################################
template<class T> VECTOR<T,0>
Vector()
{return VECTOR<T,0>();}

template<class T,class ...Args> VECTOR<T,1+sizeof...(Args)>
Vector(const T& d1,const Args&... d2)
{return VECTOR<T,1+sizeof...(Args)>(d1,d2...);}

//#####################################################################
template<class T,int d> struct SUM<VECTOR<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct SUM<VECTOR<T,d>,T>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct SUM<T,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct DIFFERENCE<VECTOR<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct DIFFERENCE<VECTOR<T,d>,T>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct DIFFERENCE<T,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<VECTOR<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<VECTOR<T,d>,T>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct PRODUCT<T,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct QUOTIENT<VECTOR<T,d>,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct QUOTIENT<VECTOR<T,d>,T>{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct QUOTIENT<T,VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
template<class T,int d> struct NEGATION<VECTOR<T,d> >{typedef VECTOR<T,d> TYPE;};
//#####################################################################
template<class T> struct HASH_REDUCE<VECTOR<T,0> >
{static int H(const VECTOR<T,0>& key){return missing_element_hash;}};
template<class T> struct HASH_REDUCE<VECTOR<T,1> >
{static int H(const VECTOR<T,1>& key){return int_hash(HASH_REDUCE<T>::H(key.x));}};
template<class T> struct HASH_REDUCE<VECTOR<T,2> >
{static int H(const VECTOR<T,2>& key){return int_hash(HASH_REDUCE<T>::H(key.x),HASH_REDUCE<T>::H(key.y));}};
template<class T> struct HASH_REDUCE<VECTOR<T,3> >
{static int H(const VECTOR<T,3>& key){return int_hash(HASH_REDUCE<T>::H(key.x),HASH_REDUCE<T>::H(key.y),HASH_REDUCE<T>::H(key.z));}};
template<class T> struct HASH_REDUCE<VECTOR<T,4> >
{static int H(const VECTOR<T,4>& key){return int_hash(int_hash(HASH_REDUCE<T>::H(key[0]),HASH_REDUCE<T>::H(key[1])),HASH_REDUCE<T>::H(key[2]),HASH_REDUCE<T>::H(key[3]));}};
template<class T,int d> struct HASH_REDUCE<VECTOR<T,d> >
{static int H(const VECTOR<T,d>& key){return Hash_Reduce_Array_Helper(key);}};
template<class T,int d>
inline std::istream& operator>>(std::istream& input,VECTOR<T,d>& v)
{
    FILE_UTILITIES::Ignore(input,'(');
    FILE_UTILITIES::Ignore(input,'[');
    for(int i=0;i<d;i++) input>>v(i);
    FILE_UTILITIES::Ignore(input,']');
    FILE_UTILITIES::Ignore(input,')');
    return input;
}

}
#endif
