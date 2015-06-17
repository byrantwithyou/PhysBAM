//#####################################################################
// Copyright 2007, Nipun Kwatra, Craig Schroeder, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_0D
//#####################################################################
#ifndef __VECTOR_0D__
#define __VECTOR_0D__

#include <Tools/Arrays/ARRAY_BASE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Math_Tools/Inverse.h>
#include <Tools/Vectors/SCALAR_POLICY.h>
#include <cmath>
#include <cstdlib>
#ifdef DIFFERENCE // Windows workaround.
#undef DIFFERENCE
#endif
namespace PhysBAM{

template<class T>
class VECTOR<T,0>:public ARRAY_BASE<T,VECTOR<T,0> >
{
    struct UNUSABLE{};
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    enum WORKAROUND1 {dimension=0};
    enum WORKAROUND2 {m=0};
    typedef typename conditional<is_scalar<T>::value,T,UNUSABLE>::type SCALAR;
    template<class T2> struct REBIND{typedef VECTOR<T2,0> TYPE;};
    typedef T ELEMENT;
    typedef UNUSABLE SPIN;
    typedef int INDEX;
    typedef ARRAY_BASE<T,VECTOR<T,0> > BASE;
    using BASE::Assert_Same_Size;

    VECTOR()
    {
    }

    explicit VECTOR(INITIAL_SIZE n)
    {
        assert(n==INITIAL_SIZE());
    }

    VECTOR(const VECTOR& vector_input)
        :BASE()
    {}

    template<class T2> explicit VECTOR(const VECTOR<T2,0>& vector_input)
    {}

    template<class T_VECTOR>
    explicit VECTOR(const ARRAY_BASE<T,T_VECTOR>& v)
    {
        Assert_Same_Size(*this,v);
    }

    VECTOR(const VECTOR& v1,const VECTOR& v2)
    {
    }

    template<class T_VECTOR>
    VECTOR& operator=(const ARRAY_BASE<T,T_VECTOR>& v)
    {
        Assert_Same_Size(*this,v);
        return *this;
    }

    VECTOR& operator=(const VECTOR& v)
    {
        return *this;
    }

    const T& operator[](const int) const
    {PHYSBAM_FATAL_ERROR();}

    T& operator[](const int)
    {PHYSBAM_FATAL_ERROR();}

    int Size() const
    {return 0;}

    const T& operator()(const int i) const
    {PHYSBAM_FATAL_ERROR();}

    T& operator()(const int i)
    {PHYSBAM_FATAL_ERROR();}

    bool operator==(const VECTOR& v) const
    {return true;}

    bool operator!=(const VECTOR& v) const
    {return false;}

    VECTOR operator-() const
    {return *this;}

    VECTOR& operator+=(const VECTOR&)
    {return *this;}

    VECTOR& operator-=(const VECTOR&)
    {return *this;}

    VECTOR& operator*=(const VECTOR&)
    {return *this;}

    VECTOR& operator-=(const T&)
    {return *this;}

    VECTOR& operator+=(const T&)
    {return *this;}

    VECTOR& operator*=(const T&)
    {return *this;}

    VECTOR& operator/=(const T&)
    {return *this;}

    VECTOR& operator/=(const VECTOR&)
    {return *this;}

    VECTOR operator+(const VECTOR&) const
    {return *this;}

    VECTOR operator+(const T&) const
    {return *this;}

    VECTOR operator-(const T&) const
    {return *this;}

    VECTOR operator-(const VECTOR&) const
    {return *this;}

    VECTOR operator*(const VECTOR&) const
    {return *this;}

    VECTOR operator/(const VECTOR&) const
    {return *this;}

    VECTOR operator*(const T&) const
    {return *this;}

    VECTOR operator/(const T&) const
    {return *this;}

    VECTOR operator*(const INT_INVERSE a) const
    {return VECTOR();}

    VECTOR& operator*=(const INT_INVERSE a)
    {return *this;}

    bool Contains(const T&) const
    {return false;}

    template<class T_ARRAY>
    bool Contains_All(const T_ARRAY& elements) const
    {STATIC_ASSERT((is_same<typename T_ARRAY::ELEMENT,T>::value));return !elements.Size();}

    template<class T_ARRAY>
    bool Contains_Any(const T_ARRAY& elements) const
    {STATIC_ASSERT((is_same<typename T_ARRAY::ELEMENT,T>::value));return false;}

    T Magnitude_Squared() const
    {return 0;}

    T Magnitude() const
    {return 0;}

    T Normalize()
    {return T();}

    VECTOR Normalized() const
    {return *this;}

    T Dot(const VECTOR&) const
    {return T();}

    static T Dot_Product(const VECTOR&,const VECTOR&)
    {return T();}

    T Sum() const
    {return 0;}

    T Average() const
    {PHYSBAM_FATAL_ERROR();}

    T Product() const
    {return 1;}

    T Min() const
    {PHYSBAM_FATAL_ERROR();}

    T Max() const
    {PHYSBAM_FATAL_ERROR();}

    T Max_Abs() const
    {return 0;}

    T Sum_Abs() const
    {return 0;}

    int Arg_Min() const
    {PHYSBAM_FATAL_ERROR();}

    int Arg_Max() const
    {PHYSBAM_FATAL_ERROR();}

    VECTOR<T,1> Append(const T& element) const
    {return VECTOR<T,1>(element);}

    VECTOR<T,1> Prepend(const T& element) const
    {return VECTOR<T,1>(element);}

    VECTOR<T,1> Insert(const T& element,const int index) const
    {VECTOR<T,1> r;r[index]=element;return r;}

    bool All_Greater(const VECTOR& v) const
    {return true;}

    bool All_Less(const VECTOR& v) const
    {return true;}

    bool All_Greater_Equal(const VECTOR& v) const
    {return true;}

    bool All_Less_Equal(const VECTOR& v) const
    {return true;}

    static VECTOR Constant_Vector(const T& constant)
    {return VECTOR();}

    static VECTOR All_Ones_Vector()
    {return VECTOR();}

    VECTOR Componentwise_Min(const VECTOR& v) const
    {return VECTOR();}

    VECTOR Componentwise_Max(const VECTOR& v) const
    {return VECTOR();}

    static VECTOR Componentwise_Min(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR();}

    static VECTOR Componentwise_Max(const VECTOR& v1,const VECTOR& v2)
    {return VECTOR();}

    VECTOR Projected_On_Unit_Direction(const VECTOR& direction) const
    {return VECTOR();}

    VECTOR Projected(const VECTOR& direction) const // un-normalized direction
    {return VECTOR();}

    void Project_On_Unit_Direction(const VECTOR& direction)
    {}

    void Project(const VECTOR& direction) // un-normalized direction
    {}

    VECTOR Projected_Orthogonal_To_Unit_Direction(const VECTOR& direction) const
    {return VECTOR();}

    void Project_Orthogonal_To_Unit_Direction(const VECTOR& direction)
    {}

    VECTOR Reversed() const
    {return *this;}

    VECTOR<T,1> Cross(const VECTOR<T,1>) const
    {return VECTOR<T,1>();}

    bool Elements_Equal() const
    {return true;}

    int Dominant_Axis() const
    {PHYSBAM_FATAL_ERROR();}

    int Number_True() const
    {STATIC_ASSERT((is_same<T,bool>::value));return 0;}

    int Find(const T& element) const
    {return -1;}

    static VECTOR Axis_Vector(const int axis)
    {PHYSBAM_FATAL_ERROR();}

    template<class T_FUNCTION>
    static VECTOR Map(const T_FUNCTION& f,const VECTOR& v)
    {return VECTOR();}

    void Split(VECTOR& v1,VECTOR& v2) const
    {}

    template<int d2> VECTOR<T,d2> Append_Elements(const VECTOR<T,d2>& elements) const
    {return elements;}

    void Fill(const T& constant)
    {}

    T Lp_Norm(const T& p) const
    {return 0;}

    VECTOR Orthogonal_Vector() const
    {PHYSBAM_FATAL_ERROR();}

    VECTOR Unit_Orthogonal_Vector() const
    {PHYSBAM_FATAL_ERROR();}

    VECTOR Sorted() const
    {return *this;}

    void Sort()
    {}

    template<class T_COMPARE>
    void Sort(const T_COMPARE comparison)
    {}

    static T Angle_Between(const VECTOR& u,const VECTOR& v)
    {return 0;}

    template<class T_VECTOR>
    void Set_Subvector(const int istart,const T_VECTOR& v)
    {assert(istart==0 && v.Size()==0);}

    template<class T_VECTOR>
    void Add_Subvector(const int istart,const T_VECTOR& v)
    {assert(istart==0 && v.Size()==0);}
    
    template<class T_VECTOR>
    void Get_Subvector(const int istart,T_VECTOR& v) const
    {assert(istart==0 && v.Size()==0);}

    T* Get_Array_Pointer()
    {return (T*)this;}

    const T* Get_Array_Pointer() const
    {return (const T*)this;}

    T* begin() // for stl
    {return (T*)this;}

    const T* begin() const // for stl
    {return (T*)this;}

    T* end() // for stl
    {return (T*)this;}

    const T* end() const // for stl
    {return (T*)this;}

    template<class RW>
    void Read(std::istream& input)
    {}

    template<class RW>
    void Write(std::ostream& output) const
    {}
};

template<class T> inline VECTOR<T,0>
operator+(const T&,const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
operator-(const T&,const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
operator*(const T&,const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
operator/(const T&,const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
Inverse(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
abs(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
floor(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
ceil(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
rint(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
exp(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
sin(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
cos(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
sqrt(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
log(const VECTOR<T,0>& v)
{return v;}

template<class T> inline VECTOR<T,0>
wrap(const VECTOR<T,0>& v,const VECTOR<T,0>& vmin,const VECTOR<T,0>& vmax)
{return VECTOR<T,0 >();}

//#####################################################################
// Functions clamp, clamp_min, clamp_max, in_bounds
//#####################################################################
template<class T> inline VECTOR<T,0>
clamp(const VECTOR<T,0>& v,const VECTOR<T,0>& vmin,const VECTOR<T,0>& vmax)
{return v;}

template<class T> inline VECTOR<T,0>
clamp(const VECTOR<T,0>& v,T min,T max)
{return v;}

template<class T> inline VECTOR<T,0>
clamp_min(const VECTOR<T,0>& v,const VECTOR<T,0>& vmin)
{return v;}

template<class T> inline VECTOR<T,0>
clamp_min(const VECTOR<T,0>& v,const T& min)
{return v;}

template<class T> inline VECTOR<T,0>
clamp_max(const VECTOR<T,0>& v,const VECTOR<T,0>& vmax)
{return v;}

template<class T> inline VECTOR<T,0>
clamp_max(const VECTOR<T,0>& v,const T& max)
{return v;}

template<class T> inline bool
in_bounds(const VECTOR<T,0>& v,const VECTOR<T,0>& vmin,const VECTOR<T,0>& vmax)
{return true;}

//#####################################################################
}
#include <Tools/Vectors/VECTOR.h>
#endif
