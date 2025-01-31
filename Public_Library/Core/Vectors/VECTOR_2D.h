//#####################################################################
// Copyright 2002-2008, Robert Bridson, Kevin Der, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Igor Neverov, Duc Nguyen, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_2D
//#####################################################################
#ifndef __VECTOR_2D__
#define __VECTOR_2D__

#include <Core/Math_Tools/argmax.h>
#include <Core/Math_Tools/argmin.h>
#include <Core/Math_Tools/clamp.h>
#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/exchange_sort.h>
#include <Core/Math_Tools/Inverse.h>
#include <Core/Math_Tools/max.h>
#include <Core/Math_Tools/maxabs.h>
#include <Core/Math_Tools/min.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Math_Tools/wrap.h>
#include <Core/Utilities/STATIC_ASSERT.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Core/Vectors/VECTOR_1D.h>
#include <cmath>
#include <cstdlib>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

using ::std::floor;
using ::std::rint;
using ::std::ceil;
using ::std::sin;
using ::std::cos;
using ::std::abs;
using ::std::sqrt;

template<class T>
class VECTOR<T,2>:public ARRAY_BASE<T,VECTOR<T,2> >
{
    struct UNUSABLE{};
public:
    typedef ARRAY_BASE<T,VECTOR<T,2> > BASE;
    using BASE::Assert_Same_Size;
    template<class T2> struct REBIND{typedef VECTOR<T2,2> TYPE;};
    typedef typename conditional<is_scalar<T>::value,T,UNUSABLE>::type SCALAR;
    typedef T ELEMENT;
    typedef VECTOR<T,1> SPIN;
    typedef int INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    typedef const T* const_iterator; // for stl
    enum WORKAROUND1 {dimension=2};
    enum WORKAROUND2 {m=2};
    typedef int HAS_UNTYPED_READ_WRITE;

    union
    {
        struct{T x,y;};
        T array[2];
    };

    VECTOR()
        :x(),y()
    {
        STATIC_ASSERT(sizeof(VECTOR)==2*sizeof(T));
    }

    explicit VECTOR(INITIAL_SIZE n)
        :x(),y()
    {
        assert(n==INITIAL_SIZE(2));
    }

    VECTOR(const T& x_input,const T& y_input)
        :x(x_input),y(y_input)
    {}

    VECTOR(const VECTOR& vector_input)
        :x(vector_input.x),y(vector_input.y)
    {}

    VECTOR(const VECTOR&& vector_input)
        :x(vector_input.x),y(vector_input.y)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,2>& vector_input)
        :x((T)vector_input.x),y((T)vector_input.y)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,1>& vector_input)
        :x((T)vector_input.x),y(T())
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,0>& vector_input)
        :x(T()),y(T())
    {}

    template<class T_VECTOR>
    explicit VECTOR(const ARRAY_BASE<T,T_VECTOR>& v)
        :x(v(0)),y(v(1))
    {
        Assert_Same_Size(*this,v);
    }

    ~VECTOR()
    {
        x.~T();
        y.~T();
    }

    template<class T_VECTOR>
    VECTOR& operator=(const ARRAY_BASE<T,T_VECTOR>& v)
    {
        Assert_Same_Size(*this,v);
        x=v(0);y=v(1);
        return *this;
    }

    VECTOR& operator=(const VECTOR& v)
    {
        x=v(0);y=v(1);return *this;
    }

    int Size() const
    {return 2;}

    const T& operator[](const int i) const
    {assert((unsigned)i<2);return array[i];}

    T& operator[](const int i)
    {assert((unsigned)i<2);return array[i];}

    const T& operator()(const int i) const
    {assert((unsigned)i<2);return array[i];}

    T& operator()(const int i)
    {assert((unsigned)i<2);return array[i];}

    bool operator==(const VECTOR& v) const
    {return x==v.x && y==v.y;}

    bool operator!=(const VECTOR& v) const
    {return x!=v.x || y!=v.y;}

    const VECTOR& operator+() const
    {return *this;};

    VECTOR operator-() const
    {return VECTOR(-x,-y);}

    VECTOR& operator+=(const VECTOR& v)
    {x+=v.x;y+=v.y;return *this;}

    VECTOR& operator-=(const VECTOR& v)
    {x-=v.x;y-=v.y;return *this;}

    VECTOR& operator*=(const VECTOR& v)
    {x*=v.x;y*=v.y;return *this;}

    VECTOR& operator+=(const T& a)
    {x+=a;y+=a;return *this;}

    VECTOR& operator-=(const T& a)
    {x-=a;y-=a;return *this;}

    VECTOR& operator*=(const T& a)
    {x*=a;y*=a;return *this;}

    VECTOR& operator*=(const INT_INVERSE a)
    {x*=a;y*=a;return *this;}

    VECTOR& operator/=(const T& a)
    {return *this*=Inverse(a);}

    VECTOR& operator/=(const VECTOR& v)
    {x/=v.x;y/=v.y;return *this;}

    VECTOR operator+(const VECTOR& v) const
    {return VECTOR(x+v.x,y+v.y);}

    VECTOR operator-(const VECTOR& v) const
    {return VECTOR(x-v.x,y-v.y);}

    VECTOR operator*(const VECTOR& v) const
    {return VECTOR(x*v.x,y*v.y);}

    VECTOR operator/(const VECTOR& v) const
    {return VECTOR(x/v.x,y/v.y);}

    VECTOR operator+(const T& a) const
    {return VECTOR(x+a,y+a);}

    VECTOR operator-(const T& a) const
    {return VECTOR(x-a,y-a);}

    VECTOR operator*(const T& a) const
    {return VECTOR(x*a,y*a);}

    VECTOR operator*(const INT_INVERSE a) const
    {return VECTOR(x*a,y*a);}

    VECTOR operator/(const T& a) const
    {return *this*Inverse(a);}

    FIXED_NUMBER<T,0> Dot(const ZERO_VECTOR<T,m>&) const
    {return FIXED_NUMBER<T,0>();}

    T Magnitude_Squared() const
    {return sqr(x)+sqr(y);}

    T Magnitude() const
    {return sqrt(sqr(x)+sqr(y));}

    T Lp_Norm(const T& p) const
    {return pow(pow(abs(x),p)+pow(abs(y),p),1/p);}

    T Sum_Abs() const
    {return abs(x)+abs(y);}

    T Normalize()
    {T magnitude=Magnitude();
    if(magnitude) *this*=1/magnitude;
    else *this=VECTOR(1,0);
    return magnitude;}

    VECTOR Normalized() const
    {T magnitude=Magnitude();
    if(magnitude) return *this*(1/magnitude);
    return VECTOR(1,0);}

    VECTOR Rotate_Clockwise_90() const
    {return VECTOR(y,-x);}

    VECTOR Rotate_Counterclockwise_90() const
    {return VECTOR(-y,x);}

    VECTOR Rotate_Counterclockwise_Multiple_90(const int n) const
    {VECTOR r(*this);
    if(n&2) r=-r;
    return n&1?r.Rotate_Counterclockwise_90():r;}

    VECTOR Perpendicular() const
    {return VECTOR(-y,x);}

    VECTOR Orthogonal_Vector() const
    {return VECTOR(-y,x);}

    VECTOR Unit_Orthogonal_Vector() const
    {return Orthogonal_Vector().Normalized();}

    T Min() const
    {return min(x,y);}

    T Max() const
    {return max(x,y);}

    T Max_Abs() const
    {return maxabs(x,y);}

    int Arg_Min() const
    {return argmin(x,y);}

    int Arg_Max() const
    {return argmax(x,y);}

    bool Elements_Equal() const
    {return x==y;}

    int Dominant_Axis() const
    {return (abs(x)>abs(y))?0:1;}

    T Dot(const VECTOR& v) const
    {return x*v.x+y*v.y;}
 
    static T Dot_Product(const VECTOR& v1,const VECTOR& v2)
    {return v1.Dot(v2);}

    VECTOR Componentwise_Min(const VECTOR& v) const
    {return VECTOR(min(x,v.x),min(y,v.y));}

    VECTOR Componentwise_Max(const VECTOR& v) const
    {return VECTOR(max(x,v.x),max(y,v.y));}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {return v1.Componentwise_Min(v2);}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {return v1.Componentwise_Max(v2);}

    bool All_Greater(const VECTOR& v) const
    {return x>v.x && y>v.y;}

    bool All_Less(const VECTOR& v) const
    {return x<v.x && y<v.y;}

    bool All_Greater_Equal(const VECTOR& v) const
    {return x>=v.x && y>=v.y;}

    bool All_Less_Equal(const VECTOR& v) const
    {return x<=v.x && y<=v.y;}

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

    static VECTOR<T,1> Cross_Product(const VECTOR& v1,const VECTOR& v2)
    {return v1.Cross(v2);}

    static VECTOR Cross_Product(const VECTOR& v1,const VECTOR<T,1>& v2) // v2 is out of plane
    {return v1.Cross(v2);}

    static VECTOR Cross_Product(const VECTOR<T,1>& v1,const VECTOR& v2) // v1 is out of plane
    {return v1.Cross(v2);}

    VECTOR<T,1> Cross(const VECTOR& v) const
    {return VECTOR<T,1>(x*v.y-y*v.x);}

    VECTOR Cross(const VECTOR<T,1>& v) const // v2 is out of plane
    {return VECTOR(y*v.x,-x*v.x);}

    static T Angle_Between(const VECTOR& u,const VECTOR& v)
    {T s=Cross_Product(u,v).Magnitude(),c=Dot_Product(u,v);return atan2(s,c);}

    static T Oriented_Angle_Between(const VECTOR& u,const VECTOR& v)
    {T s=Cross_Product(u,v).x,c=Dot_Product(u,v);return atan2(s,c);}

    VECTOR Givens_Transpose_Times(const VECTOR& u) const
    {return VECTOR(x*u.x+y*u.y,x*u.y-y*u.x);}

    T Sum() const
    {return x+y;}

    T Average() const
    {return Sum()/2;}

    T Product() const
    {return x*y;}

    int Number_True() const
    {STATIC_ASSERT((is_same<T,bool>::value));return x+y;}

    static VECTOR Axis_Vector(const int axis)
    {VECTOR vec;vec[axis]=(T)1;return vec;}

    VECTOR Add_Axis(int axis,T value) const
    {VECTOR vec=*this;vec(axis)+=value;return vec;}

    static VECTOR Constant_Vector(const T& constant)
    {return VECTOR(constant,constant);}

    static VECTOR All_Ones_Vector()
    {return Constant_Vector(1);}

    void Fill(const T& constant)
    {x=y=constant;}

    void Get(T& element1,T& element2) const
    {element1=x;element2=y;}

    void Set(const T& element1,const T& element2)
    {x=element1;y=element2;}

    int Find(const T& element) const
    {return x==element?0:y==element?1:-1;}

    bool Contains(const T& element) const
    {return x==element || y==element;}

    template<class T_ARRAY>
    bool Contains_All(const T_ARRAY& elements) const
    {STATIC_ASSERT((is_same<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=0;i<elements.Size();i++) if(!Contains(elements(i))) return false;
    return true;}

    template<class T_ARRAY>
    bool Contains_Any(const T_ARRAY& elements) const
    {STATIC_ASSERT((is_same<typename T_ARRAY::ELEMENT,T>::value));
    for(int i=0;i<elements.Size();i++) if(Contains(elements(i))) return true;
    return false;}

    VECTOR<T,1> Remove_Index(const int index) const
    {assert((unsigned)index<2);return VECTOR<T,1>((*this)[1-index]);}

    VECTOR<T,3> Insert(const T& element,const int index) const
    {VECTOR<T,3> r;r[index]=element;
    for(int i=0;i<2;i++) r[i+(i>=index)]=(*this)[i];
    return r;}

    VECTOR<T,3> Prepend(const T& element) const
    {return VECTOR<T,3>(element,x,y);}

    VECTOR<T,3> Append(const T& element) const
    {return VECTOR<T,3>(x,y,element);}

    template<int d2> VECTOR<T,2+d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {VECTOR<T,2+d2> r;r[0]=x;r[1]=y;
    for(int i=0;i<d2;i++) r[i+2]=elements[i];
    return r;}

    VECTOR<T,2> Sorted() const
    {VECTOR<T,2> r(*this);exchange_sort(r.x,r.y);return r;}

    void Sort()
    {exchange_sort(x,y);}

    template<class T_COMPARE>
    void Sort(const T_COMPARE comparison)
    {exchange_sort(x,y,comparison);}

    VECTOR Reversed() const
    {return VECTOR(y,x);}

    template<int d1,int d2> VECTOR<T,d2-d1+1> Slice() const
    {STATIC_ASSERT(((0<=d1) && (d2<2)));
    VECTOR<T,d2-d1+1> r;
    for(int i=d1;i<=d2;i++) r[i-d1]=(*this)[i];
    return r;}

    template<int n> void Split(VECTOR<T,n>& v1,VECTOR<T,2-n>& v2) const
    {for(int i=0;i<n;i++) v1(i)=(*this)(i);
    for(int i=n+1;i<2;i++) v2(i-n)=(*this)(i);}

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
    {return &x;}

    const T* begin() const // for stl
    {return &x;}

    T* end() // for stl
    {return &y+1;}

    const T* end() const // for stl
    {return &y+1;}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,x,y);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x,y);}
//#####################################################################
};

//#####################################################################
// Miscellaneous free operators and functions
//#####################################################################
template<class T> inline VECTOR<T,2>
operator+(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a+v.x,a+v.y);}

template<class T> inline VECTOR<T,2>
operator-(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a-v.x,a-v.y);}

template<class T> inline VECTOR<T,2>
operator*(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a*v.x,a*v.y);}

template<class T> inline VECTOR<T,2>
operator/(const T& a,const VECTOR<T,2>& v)
{return VECTOR<T,2>(a/v.x,a/v.y);}

template<class T> inline VECTOR<T,2>
abs(const VECTOR<T,2>& v)
{return VECTOR<T,2>(abs(v.x),abs(v.y));}

template<class T> inline VECTOR<T,2>
floor(const VECTOR<T,2>& v)
{return VECTOR<T,2>(floor(v.x),floor(v.y));}

template<class T> inline VECTOR<T,2>
ceil(const VECTOR<T,2>& v)
{return VECTOR<T,2>(ceil(v.x),ceil(v.y));}

template<class T> inline VECTOR<T,2>
rint(const VECTOR<T,2>& v)
{return VECTOR<T,2>(rint(v.x),rint(v.y));}

inline VECTOR<int,2>
fdiv(const VECTOR<int,2>& a,const VECTOR<int,2>& b)
{return VECTOR<int,2>(fdiv(a.x,b.x),fdiv(a.y,b.y));}

inline VECTOR<int,2>
fdiv(const VECTOR<int,2>& a,int b)
{return VECTOR<int,2>(fdiv(a.x,b),fdiv(a.y,b));}

inline VECTOR<int,2>
cdiv(const VECTOR<int,2>& a,const VECTOR<int,2>& b)
{return VECTOR<int,2>(cdiv(a.x,b.x),cdiv(a.y,b.y));}

inline VECTOR<int,2>
cdiv(const VECTOR<int,2>& a,int b)
{return VECTOR<int,2>(cdiv(a.x,b),cdiv(a.y,b));}

template<class T> inline VECTOR<T,2>
exp(const VECTOR<T,2>& v)
{return VECTOR<T,2>(exp(v.x),exp(v.y));}

template<class T> inline VECTOR<T,2>
sin(const VECTOR<T,2>& v)
{return VECTOR<T,2>(sin(v.x),sin(v.y));}

template<class T> inline VECTOR<T,2>
cos(const VECTOR<T,2>& v)
{return VECTOR<T,2>(cos(v.x),cos(v.y));}

template<class T> inline VECTOR<T,2>
sqrt(const VECTOR<T,2>& v)
{return VECTOR<T,2>(sqrt(v.x),sqrt(v.y));}

template<class T> inline VECTOR<T,2>
log(const VECTOR<T,2>& v)
{return VECTOR<T,2>(log(v.x),log(v.y));}

template<class T> inline VECTOR<T,2>
Inverse(const VECTOR<T,2>& v)
{return VECTOR<T,2>(1/v.x,1/v.y);}

//#####################################################################
// Functions clamp, clamp_min, clamp_max, in_bounds
//#####################################################################
template<class T> inline VECTOR<T,2>
clamp(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin,const VECTOR<T,2>& vmax)
{return VECTOR<T,2>(clamp(v.x,vmin.x,vmax.x),clamp(v.y,vmin.y,vmax.y));}

template<class T> inline VECTOR<T,2>
clamp(const VECTOR<T,2>& v,T min,T max)
{return VECTOR<T,2>(clamp(v.x,min,max),clamp(v.y,min,max));}

template<class T> inline VECTOR<T,2>
clamp_min(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin)
{return VECTOR<T,2>(clamp_min(v.x,vmin.x),clamp_min(v.y,vmin.y));}

template<class T> inline VECTOR<T,2>
clamp_min(const VECTOR<T,2>& v,const T& min)
{return VECTOR<T,2>(clamp_min(v.x,min),clamp_min(v.y,min));}

template<class T> inline VECTOR<T,2>
clamp_max(const VECTOR<T,2>& v,const VECTOR<T,2>& vmax)
{return VECTOR<T,2>(clamp_max(v.x,vmax.x),clamp_max(v.y,vmax.y));}

template<class T> inline VECTOR<T,2>
clamp_max(const VECTOR<T,2>& v,const T& max)
{return VECTOR<T,2>(clamp_max(v.x,max),clamp_max(v.y,max));}

template<class T> inline bool
in_bounds(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin,const VECTOR<T,2>& vmax)
{return in_bounds(v.x,vmin.x,vmax.x) && in_bounds(v.y,vmin.y,vmax.y);}

template<class T> inline VECTOR<T,2>
wrap(const VECTOR<T,2>& v,const VECTOR<T,2>& vmin,const VECTOR<T,2>& vmax)
{return VECTOR<T,2>(wrap(v.x,vmin.x,vmax.x),wrap(v.y,vmin.y,vmax.y));}
//#####################################################################
}
#endif
