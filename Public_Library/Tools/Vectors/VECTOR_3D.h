//#####################################################################
// Copyright 2002-2008, Robert Bridson, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Igor Neverov, Duc Nguyen, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_3D
//#####################################################################
#ifndef __VECTOR_3D__
#define __VECTOR_3D__
#include <Tools/Math_Tools/argmax.h>
#include <Tools/Math_Tools/argmin.h>
#include <Tools/Math_Tools/clamp.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Math_Tools/Inverse.h>
#include <Tools/Math_Tools/max.h>
#include <Tools/Math_Tools/maxabs.h>
#include <Tools/Math_Tools/min.h>
#include <Tools/Math_Tools/wrap.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <cmath>
#include <cstdlib>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

using ::std::exp;
using ::std::sin;
using ::std::cos;
using ::std::abs;
using ::std::floor;
using ::std::ceil;
using ::std::sqrt;

template<class T_ARRAY,class T_INDICES> class INDIRECT_ARRAY;

template<class T>
class VECTOR<T,3>:public ARRAY_BASE<T,VECTOR<T,3> >
{
    struct UNUSABLE{};
public:
    typedef ARRAY_BASE<T,VECTOR<T,3> > BASE;
    using BASE::Assert_Same_Size;
    typedef int HAS_UNTYPED_READ_WRITE;
    template<class T2> struct REBIND{typedef VECTOR<T2,3> TYPE;};
    typedef typename conditional<is_scalar<T>::value,T,UNUSABLE>::type SCALAR;
    typedef T ELEMENT;
    typedef VECTOR<T,3> SPIN;
    typedef int INDEX;
    typedef T value_type; // for stl
    typedef T* iterator; // for stl
    enum WORKAROUND1 {dimension=3};
    enum WORKAROUND2 {m=3};

    union
    {
        struct{T x,y,z;};
        T array[3];
    };

    VECTOR()
        :x(),y(),z()
    {
        STATIC_ASSERT(sizeof(VECTOR)==3*sizeof(T));
    }

    explicit VECTOR(INITIAL_SIZE n)
        :x(),y(),z()
    {
        assert(n==INITIAL_SIZE(3));
    }

    VECTOR(const T& x_input,const T& y_input,const T& z_input)
        :x(x_input),y(y_input),z(z_input)
    {}

    VECTOR(const VECTOR& vector_input)
        :x(vector_input.x),y(vector_input.y),z(vector_input.z)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,3>& vector_input)
        :x((T)vector_input.x),y((T)vector_input.y),z((T)vector_input.z)
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,2>& vector_input)
        :x((T)vector_input.x),y((T)vector_input.y),z()
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,1>& vector_input)
        :x((T)vector_input.x),y(),z()
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,0>& vector_input)
        :x(),y(),z()
    {}

    template<class T_VECTOR>
    explicit VECTOR(const ARRAY_BASE<T,T_VECTOR>& v)
        :x(v(0)),y(v(1)),z(v(2))
    {
        Assert_Same_Size(*this,v);
    }

    template<int n>
    VECTOR(const VECTOR<T,n>& v1,const VECTOR<T,3-n>& v2)
        :x(),y(),z()
    {
        for(int i=0;i<n;i++) (*this)(i)=v1(i);for(int i=n;i<3;i++) (*this)(i)=v2(i-n);
    }

    ~VECTOR()
    {
        x.~T();
        y.~T();
        z.~T();
    }

    template<class T_VECTOR>
    VECTOR& operator=(const ARRAY_BASE<T,T_VECTOR>& v)
    {
        Assert_Same_Size(*this,v);
        x=v(0);y=v(1);z=v(2);
        return *this;
    }

    VECTOR& operator=(const VECTOR& v)
    {
        x=v(0);y=v(1);z=v(2);
        return *this;
    }

    int Size() const
    {return 3;}

    const T& operator[](const int i) const
    {assert((unsigned)i<3);return array[i];}

    T& operator[](const int i)
    {assert((unsigned)i<3);return array[i];}

    const T& operator()(const int i) const
    {assert((unsigned)i<3);return array[i];}

    T& operator()(const int i)
    {assert((unsigned)i<3);return array[i];}

    bool operator==(const VECTOR& v) const
    {return x==v.x && y==v.y && z==v.z;}

    bool operator!=(const VECTOR& v) const
    {return x!=v.x || y!=v.y || z!=v.z;}

    VECTOR operator-() const
    {return VECTOR(-x,-y,-z);}

    VECTOR& operator+=(const VECTOR& v)
    {x+=v.x;y+=v.y;z+=v.z;return *this;}

    VECTOR& operator-=(const VECTOR& v)
    {x-=v.x;y-=v.y;z-=v.z;return *this;}

    VECTOR& operator*=(const VECTOR& v)
    {x*=v.x;y*=v.y;z*=v.z;return *this;}

    VECTOR& operator+=(const T& a)
    {x+=a;y+=a;z+=a;return *this;}

    VECTOR& operator-=(const T& a)
    {x-=a;y-=a;z-=a;return *this;}

    VECTOR& operator*=(const T& a)
    {x*=a;y*=a;z*=a;return *this;}

    VECTOR& operator*=(const INT_INVERSE a)
    {x*=a;y*=a;z*=a;return *this;}

    VECTOR& operator/=(const T& a)
    {return *this*=Inverse(a);}

    VECTOR& operator/=(const VECTOR& v)
    {x/=v.x;y/=v.y;z/=v.z;return *this;}

    VECTOR operator+(const VECTOR& v) const
    {return VECTOR(x+v.x,y+v.y,z+v.z);}

    VECTOR operator-(const VECTOR& v) const
    {return VECTOR(x-v.x,y-v.y,z-v.z);}

    VECTOR operator*(const VECTOR& v) const
    {return VECTOR(x*v.x,y*v.y,z*v.z);}

    VECTOR operator/(const VECTOR& v) const
    {return VECTOR(x/v.x,y/v.y,z/v.z);}

    VECTOR operator+(const T& a) const
    {return VECTOR(x+a,y+a,z+a);}

    VECTOR operator-(const T& a) const
    {return VECTOR(x-a,y-a,z-a);}

    VECTOR operator*(const T& a) const
    {return VECTOR(x*a,y*a,z*a);}

    VECTOR operator*(const INT_INVERSE a) const
    {return VECTOR(x*a,y*a,z*a);}

    VECTOR operator/(const T& a) const
    {return *this*Inverse(a);}

    FIXED_NUMBER<T,0> Dot(const ZERO_VECTOR<T,m>&) const
    {return FIXED_NUMBER<T,0>();}

    T Magnitude_Squared() const
    {return x*x+y*y+z*z;}

    T Magnitude() const
    {return sqrt(x*x+y*y+z*z);}

    T Lp_Norm(const T& p) const
    {return pow(pow(abs(x),p)+pow(abs(y),p)+pow(abs(z),p),1/p);}

    T Sum_Abs() const
    {return abs(x)+abs(y)+abs(z);}

    T Normalize()
    {T magnitude=Magnitude();if(magnitude) *this*=1/magnitude;else *this=VECTOR(1,0,0);return magnitude;}

    VECTOR Normalized() const // 6 mults, 2 adds, 1 div, 1 sqrt
    {T magnitude=Magnitude();if(magnitude) return *this*(1/magnitude);else return VECTOR(1,0,0);}

    VECTOR Orthogonal_Vector() const
    {T abs_x=abs(x),abs_y=abs(y),abs_z=abs(z);
    if(abs_x<abs_y) return abs_x<abs_z?VECTOR(0,z,-y):VECTOR(y,-x,0);
    else return abs_y<abs_z?VECTOR(-z,0,x):VECTOR(y,-x,0);}

    VECTOR Unit_Orthogonal_Vector() const // roughly 6 mults, 2 adds, 1 div, 1 sqrt
    {return Orthogonal_Vector().Normalized();}

    T Min() const
    {return min(x,y,z);}

    T Max() const
    {return max(x,y,z);}

    T Max_Abs() const
    {return maxabs(x,y,z);}

    int Arg_Min() const
    {return argmin(x,y,z);}

    int Arg_Max() const
    {return argmax(x,y,z);}

    bool Elements_Equal() const
    {return x==y && x==z;}

    bool All_Greater(const VECTOR& v) const
    {return x>v.x && y>v.y && z>v.z;}

    bool All_Less(const VECTOR& v) const
    {return x<v.x && y<v.y && z<v.z;}

    bool All_Greater_Equal(const VECTOR& v) const
    {return x>=v.x && y>=v.y && z>=v.z;}

    bool All_Less_Equal(const VECTOR& v) const
    {return x<=v.x && y<=v.y && z<=v.z;}

    int Dominant_Axis() const
    {T abs_x=abs(x),abs_y=abs(y),abs_z=abs(z);return (abs_x>abs_y && abs_x>abs_z)?0:((abs_y>abs_z)?1:2);}

    static T Dot_Product(const VECTOR& v1,const VECTOR& v2)
    {return v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;}

    T Dot(const VECTOR& v) const
    {return x*v.x+y*v.y+z*v.z;}

    VECTOR Componentwise_Min(const VECTOR& v) const
    {return VECTOR(min(x,v.x),min(y,v.y),min(z,v.z));}

    VECTOR Componentwise_Max(const VECTOR& v) const
    {return VECTOR(max(x,v.x),max(y,v.y),max(z,v.z));}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {return v1.Componentwise_Min(v2);}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {return v1.Componentwise_Max(v2);}

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

    static VECTOR Cross_Product(const VECTOR& v1,const VECTOR& v2) // 6 mults, 3 adds
    {return v1.Cross(v2);}

    VECTOR Cross(const VECTOR& v) const // 6 mults, 3 adds
    {return VECTOR(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);}

    static T Angle_Between(const VECTOR& u,const VECTOR& v) // 0 .. pi
    {T s=Cross_Product(u,v).Magnitude(),c=Dot_Product(u,v);return atan2(s,c);}

    static T Triple_Product(const VECTOR& u,const VECTOR& v,const VECTOR& w)
    {return Dot_Product(u,Cross_Product(v,w));}

    T Sum() const
    {return x+y+z;}

    T Average() const
    {return ((T)1/3)*Sum();}

    T Product() const
    {return x*y*z;}

    int Number_True() const
    {STATIC_ASSERT((is_same<T,bool>::value));return x+y+z;}

    static VECTOR Axis_Vector(const int axis)
    {VECTOR vec;vec[axis]=(T)1;return vec;}

    static VECTOR Constant_Vector(const T& constant)
    {return VECTOR(constant,constant,constant);}

    static VECTOR All_Ones_Vector()
    {return Constant_Vector(1);}

    void Fill(const T& constant)
    {x=y=z=constant;}

    void Get(T& element1,T& element2,T& element3) const
    {element1=x;element2=y;element3=z;}

    void Set(const T& element1,const T& element2,const T& element3)
    {x=element1;y=element2;z=element3;}

    template<class T_FUNCTION>
    static VECTOR Map(const T_FUNCTION& f,const VECTOR& v)
    {return VECTOR(f(v.x),f(v.y),f(v.z));}

    int Find(const T& element) const
    {return x==element?0:y==element?1:z==element?2:-1;}

    bool Contains(const T& element) const
    {return x==element || y==element || z==element;}

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

    VECTOR<T,2> Remove_Index(const int index) const
    {assert((unsigned)index<3);return VECTOR<T,2>(index>0?x:y,index<2?z:y);}

    VECTOR<T,4> Insert(const T& element,const int index) const
    {VECTOR<T,4> r;r[index]=element;for(int i=0;i<3;i++) r[i+(i>=index)]=(*this)[i];return r;}

    VECTOR<T,4> Prepend(const T& element) const
    {return VECTOR<T,4>(element,x,y,z);}

    VECTOR<T,4> Append(const T& element) const
    {return VECTOR<T,4>(x,y,z,element);}

    template<int d2> VECTOR<T,3+d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {VECTOR<T,3+d2> r;r[0]=x;r[1]=y;r[2]=z;for(int i=0;i<d2;i++) r[i+3]=elements[i];return r;}

    VECTOR Sorted() const
    {VECTOR r(*this);exchange_sort(r.x,r.y,r.z);return r;}

    void Sort()
    {exchange_sort(x,y,z);}

    template<class T_COMPARE>
    void Sort(const T_COMPARE comparison)
    {exchange_sort(x,y,z,comparison);}

    VECTOR Reversed() const
    {return VECTOR(z,y,x);}

    template<int d1,int d2> VECTOR<T,d2-d1+1> Slice() const
    {STATIC_ASSERT(((0<=d1) && (d2<3)));
    VECTOR<T,d2-d1+1> r;for(int i=d1;i<=d2;i++) r[i-d1]=(*this)[i];return r;}

    template<int n> void Split(VECTOR<T,n>& v1,VECTOR<T,3-n>& v2) const
    {for(int i=0;i<n;i++) v1(i)=(*this)(i);
    for(int i=n;i<3;i++) v2(i-n)=(*this)(i);}

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
    {return &z+1;}

    const T* end() const // for stl
    {return &z+1;}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,x,y,z);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,x,y,z);}
//#####################################################################
};

//#####################################################################
// Miscellaneous free operators and functions
//#####################################################################
template<class T> inline VECTOR<T,3>
operator+(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a+v.x,a+v.y,a+v.z);}

template<class T> inline VECTOR<T,3>
operator-(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a-v.x,a-v.y,a-v.z);}

template<class T> inline VECTOR<T,3>
operator*(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a*v.x,a*v.y,a*v.z);}

template<class T> inline VECTOR<T,3>
operator/(const T& a,const VECTOR<T,3>& v)
{return VECTOR<T,3>(a/v.x,a/v.y,a/v.z);}

template<class T> inline VECTOR<T,3>
abs(const VECTOR<T,3>& v)
{return VECTOR<T,3>(abs(v.x),abs(v.y),abs(v.z));}

template<class T> inline VECTOR<T,3>
floor(const VECTOR<T,3>& v)
{return VECTOR<T,3>(floor(v.x),floor(v.y),floor(v.z));}

template<class T> inline VECTOR<T,3>
ceil(const VECTOR<T,3>& v)
{return VECTOR<T,3>(ceil(v.x),ceil(v.y),ceil(v.z));}

template<class T> inline VECTOR<T,3>
rint(const VECTOR<T,3>& v)
{return VECTOR<T,3>(rint(v.x),rint(v.y),rint(v.z));}

template<class T> inline VECTOR<T,3>
exp(const VECTOR<T,3>& v)
{return VECTOR<T,3>(exp(v.x),exp(v.y),exp(v.z));}

template<class T> inline VECTOR<T,3>
sin(const VECTOR<T,3>& v)
{return VECTOR<T,3>(sin(v.x),sin(v.y),sin(v.z));}

template<class T> inline VECTOR<T,3>
cos(const VECTOR<T,3>& v)
{return VECTOR<T,3>(cos(v.x),cos(v.y),cos(v.z));}

template<class T> inline VECTOR<T,3>
sqrt(const VECTOR<T,3>& v)
{return VECTOR<T,3>(sqrt(v.x),sqrt(v.y),sqrt(v.z));}

template<class T> inline VECTOR<T,3>
log(const VECTOR<T,3>& v)
{return VECTOR<T,3>(log(v.x),log(v.y),log(v.z));}

template<class T> inline VECTOR<T,3>
Inverse(const VECTOR<T,3>& v)
{return VECTOR<T,3>(1/v.x,1/v.y,1/v.z);}

//#####################################################################
// Functions clamp, clamp_min, clamp_max, in_bounds
//#####################################################################
template<class T> inline VECTOR<T,3>
clamp(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin,const VECTOR<T,3>& vmax)
{return VECTOR<T,3>(clamp(v.x,vmin.x,vmax.x),clamp(v.y,vmin.y,vmax.y),clamp(v.z,vmin.z,vmax.z));}

template<class T> inline VECTOR<T,3>
clamp(const VECTOR<T,3>& v,T min,T max)
{return VECTOR<T,3>(clamp(v.x,min,max),clamp(v.y,min,max),clamp(v.z,min,max));}

template<class T> inline VECTOR<T,3>
clamp_min(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin)
{return VECTOR<T,3>(clamp_min(v.x,vmin.x),clamp_min(v.y,vmin.y),clamp_min(v.z,vmin.z));}

template<class T> inline VECTOR<T,3>
clamp_min(const VECTOR<T,3>& v,const T& min)
{return VECTOR<T,3>(clamp_min(v.x,min),clamp_min(v.y,min),clamp_min(v.z,min));}

template<class T> inline VECTOR<T,3>
clamp_max(const VECTOR<T,3>& v,const VECTOR<T,3>& vmax)
{return VECTOR<T,3>(clamp_max(v.x,vmax.x),clamp_max(v.y,vmax.y),clamp_max(v.z,vmax.z));}

template<class T> inline VECTOR<T,3>
clamp_max(const VECTOR<T,3>& v,const T& max)
{return VECTOR<T,3>(clamp_max(v.x,max),clamp_max(v.y,max),clamp_max(v.z,max));}

template<class T> inline bool
in_bounds(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin,const VECTOR<T,3>& vmax)
{return in_bounds(v.x,vmin.x,vmax.x) && in_bounds(v.y,vmin.y,vmax.y) && in_bounds(v.z,vmin.z,vmax.z);}

template<class T> inline VECTOR<T,3>
wrap(const VECTOR<T,3>& v,const VECTOR<T,3>& vmin,const VECTOR<T,3>& vmax)
{return VECTOR<T,3>(wrap(v.x,vmin.x,vmax.x),wrap(v.y,vmin.y,vmax.y),wrap(v.z,vmin.z,vmax.z));}

//#####################################################################
}
#endif
