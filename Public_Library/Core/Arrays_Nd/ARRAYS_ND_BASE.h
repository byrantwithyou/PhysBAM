//#####################################################################
// Copyright 2007-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAYS_ND_BASE
//#####################################################################
#ifndef __ARRAYS_ND_BASE__
#define __ARRAYS_ND_BASE__

#include <Core/Arrays/ARRAY_BASE.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,class T_ARRAY,int dimension>
class ARRAY_BASE<T,T_ARRAY,VECTOR<int,dimension> >
{
public:
    typedef VECTOR<int,dimension> TV_INT;

    T_ARRAY& Derived()
    {return static_cast<T_ARRAY&>(*this);}

    const T_ARRAY& Derived() const
    {return static_cast<const T_ARRAY&>(*this);}

    template<class T_ARRAY0> static bool Same_Array(const T_ARRAY0& array0,const T_ARRAY0& array1)
    {return &array0==&array1;}

    template<class T_ARRAY0,class T_ARRAY1> static bool
    Same_Array(const T_ARRAY0& array0,const T_ARRAY1& array1)
    {return SAME_ARRAY<T_ARRAY0,T_ARRAY1>::Same_Array(array0,array1);}

    template<class T_INDICES>
    INDIRECT_ARRAY<T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices)
    {return INDIRECT_ARRAY<T_ARRAY,T_INDICES&>(Derived(),indices);}

    template<class T_INDICES>
    INDIRECT_ARRAY<const T_ARRAY,T_INDICES&> Subset(const T_INDICES& indices) const
    {return INDIRECT_ARRAY<const T_ARRAY,T_INDICES&>(Derived(),indices);}

protected:
    T_ARRAY& operator=(const ARRAY_BASE& source)
    {T_ARRAY& self=Derived();RANGE<TV_INT> domain_indices=self.Domain_Indices();const T_ARRAY& source_=source.Derived();assert(domain_indices==source_.Domain_Indices());
    if(!T_ARRAY::Same_Array(self,source_)) for(RANGE_ITERATOR<dimension> iterator(domain_indices);iterator.Valid();iterator.Next()) self(iterator.index)=source_(iterator.index);
    return self;}

    template<class T_ARRAY0>
    T_ARRAY& operator=(const T_ARRAY0& source)
    {STATIC_ASSERT(CAN_ASSIGN<T,typename T_ARRAY0::ELEMENT>::value);
    T_ARRAY& self=Derived();RANGE<TV_INT> domain_indices=self.Domain_Indices();assert(domain_indices==source.Domain_Indices());
    if(!T_ARRAY::Same_Array(self,source)) for(RANGE_ITERATOR<dimension> iterator(domain_indices);iterator.Valid();iterator.Next()) self(iterator.index)=source(iterator.index);
    return self;}

    ARRAY_BASE& operator=(ARRAY_BASE&&) = default;
};

template<class T,class TV_INT>
class ARRAYS_ND_BASE:public ARRAY_BASE<T,ARRAYS_ND_BASE<T,TV_INT>,TV_INT>
{
    struct UNUSABLE{};
public:
    enum WORKAROUND {d=TV_INT::m};
    typedef T ELEMENT;typedef T& RESULT_TYPE;typedef const T& CONST_RESULT_TYPE;
    typedef typename SCALAR_POLICY<T>::TYPE SCALAR;
    template<class T2> struct REBIND{typedef ARRAYS_ND_BASE<T2,VECTOR<int,d> > TYPE;};
    typedef ARRAY_BASE<T,ARRAYS_ND_BASE<T,TV_INT>,TV_INT> BASE;
    typedef TV_INT INDEX;

    RANGE<TV_INT> domain;
    TV_INT stride;
    int offset;
    ARRAY_VIEW<T> array; // one-dimensional data storage

protected:
    ARRAYS_ND_BASE()
        :domain(TV_INT(),TV_INT()),offset(0),array(0,0)
    {}

    ARRAYS_ND_BASE(const RANGE<TV_INT>& range)
        :domain(TV_INT(),TV_INT()),array(0,0)
    {Calculate_Acceleration_Constants(range);}

    ARRAYS_ND_BASE(const ARRAYS_ND_BASE&) = default;
    ARRAYS_ND_BASE(ARRAYS_ND_BASE&&) = default;
    
    void Calculate_Acceleration_Constants(const RANGE<TV_INT>& range)
    {domain=range;
    TV_INT counts(domain.Edge_Lengths());
    stride(TV_INT::m-1)=1;
    for(int i=TV_INT::m-1;i>0;i--) stride(i-1)=stride(i)*counts(i);
    offset=-domain.min_corner.Dot(stride);}

public:
    template<class T_ARRAY1>
    ARRAYS_ND_BASE& operator=(const T_ARRAY1& source)
    {assert(Equal_Dimensions(*this,source));array=source.array;return *this;}

    ARRAYS_ND_BASE& operator=(const ARRAYS_ND_BASE& source)
    {assert(Equal_Dimensions(*this,source));array=source.array;return *this;}

    ARRAYS_ND_BASE& operator=(ARRAYS_ND_BASE&& source)
    {assert(Equal_Dimensions(*this,source));array=std::move(source.array);return *this;}

    const RANGE<TV_INT>& Domain_Indices() const
    {return domain;}

    const TV_INT Size() const
    {return domain.Edge_Lengths();}

    bool operator==(const ARRAYS_ND_BASE& v) const
    {return Equal_Dimensions(*this,v) && array==v.array;}

    bool operator!=(const ARRAYS_ND_BASE& v) const
    {return !(*this==v);}

    ARRAYS_ND_BASE& operator+=(const ARRAYS_ND_BASE& v)
    {assert(Equal_Dimensions(*this,v));array+=v.array;return *this;}

    ARRAYS_ND_BASE& operator+=(const T& a)
    {array+=a;return *this;}

    ARRAYS_ND_BASE& operator-=(const ARRAYS_ND_BASE& v)
    {assert(Equal_Dimensions(*this,v));array-=v.array;return *this;}

    ARRAYS_ND_BASE& operator-=(const T& a)
    {array-=a;return *this;}

    template<class T2>
    ARRAYS_ND_BASE& operator*=(const ARRAYS_ND_BASE<T2,VECTOR<int,d> >& v)
    {assert(Equal_Dimensions(*this,v));array*=v.array;return *this;}

    ARRAYS_ND_BASE& operator*=(const T& a)
    {array*=a;return *this;}

    template<class T2> typename enable_if<is_convertible<T2,T>::value,ARRAYS_ND_BASE&>::type
    operator*=(const T2 a)
    {array*=(T)a;return *this;}

    template<class T2> typename enable_if<is_scalar<T2>::value && !is_convertible<T2,T>::value,ARRAYS_ND_BASE&>::type
    operator*=(const T2 a)
    {array*=a;return *this;}

    template<class T2>
    ARRAYS_ND_BASE& operator/=(const T2 a)
    {return *this*=Inverse(a);}

    T* Get_Array_Pointer()
    {return array.Get_Array_Pointer();}

    const T* Get_Array_Pointer() const
    {return array.Get_Array_Pointer();}

    int Number_True() const
    {return array.Number_True();}

    int Pack_Size() const
    {return sizeof(T);}

    template<class T_ARRAY0>
    void Pack(T_ARRAY0& buffer,typename T_ARRAY0::INDEX& position,const TV_INT& p)
    {*(T*)(&buffer(position+1))=(*this)(p);position+=sizeof(T);}
    
    template<class T_ARRAY0>
    void Unpack(const T_ARRAY0& buffer,typename T_ARRAY0::INDEX& position,const TV_INT& p)
    {(*this)(p)=*(const T*)(&buffer(position+1));position+=sizeof(T);}

    void Fill(const T& constant)
    {array.Fill(constant);}

    void Copy(const ARRAYS_ND_BASE& old_copy)
    {assert(old_copy.Domain_Indices()==Domain_Indices());array.Copy(old_copy.array);}

    template<class T2>
    void Copy(const T2 constant,const ARRAYS_ND_BASE& old_copy)
    {assert(old_copy.Domain_Indices()==Domain_Indices());array=constant*old_copy.array;}

    template<class T2>
    void Copy(const T2 c1,const ARRAYS_ND_BASE& v1,const ARRAYS_ND_BASE& v2)
    {assert(Equal_Dimensions(v1,v2)&&Equal_Dimensions(v2,*this));array=c1*v1.array+v2.array;}

    template<class T2>
    void Copy(const T2 c1,const ARRAYS_ND_BASE& v1,const T2 c2,const ARRAYS_ND_BASE& v2)
    {assert(Equal_Dimensions(v1,v2)&&Equal_Dimensions(v2,*this));array=c1*v1.array+c2*v2.array;}

    template<class T2>
    void Copy(const T2 c1,const ARRAYS_ND_BASE& v1,const T2 c2,const ARRAYS_ND_BASE& v2,const T2 c3,const ARRAYS_ND_BASE& v3)
    {assert(Equal_Dimensions(v1,v2)&&Equal_Dimensions(v2,v3)&&Equal_Dimensions(v3,*this));array=c1*v1.array+c2*v2.array+c3*v3.array;}

    void Clamp_Below(const T& value)
    {array.Clamp_Below(value);}

    T Average() const
    {return array.Average();}

    T Max() const
    {return array.Max();}

    T Max_Abs() const
    {return array.Max_Abs();}

    T Maxmag() const
    {return array.Maxmag();}

    T Min() const
    {return array.Min();}

    T Minmag() const
    {return array.Minmag();}

    T Sum() const
    {return array.Sum();}

    T Sumabs() const
    {return array.Sumabs();}

    T Componentwise_Max_Abs() const
    {return array.Componentwise_Max_Abs();}

    static T Dot_Product(const ARRAYS_ND_BASE& a1,const ARRAYS_ND_BASE& a2)
    {assert(Equal_Dimensions(a1,a2));
    return ARRAY_VIEW<T>::Dot_Product(a1.array,a2.array);}

    template<class TV2>
    static typename SCALAR_POLICY<T>::TYPE Maximum_Magnitude(const ARRAYS_ND_BASE<TV2,VECTOR<int,d> >& a)
    {STATIC_ASSERT(is_same<T,TV2>::value);
    return ARRAY_VIEW<T>::Maximum_Magnitude(a.array);}

    T& operator()(const int i,const int j,const int k)
    {STATIC_ASSERT(d==3);return (*this)(VECTOR<int,3>(i,j,k));}

    const T& operator()(const int i,const int j,const int k) const
    {STATIC_ASSERT(d==3);return (*this)(VECTOR<int,3>(i,j,k));}

    T& operator()(const int i,const int j)
    {STATIC_ASSERT(d==2);return (*this)(VECTOR<int,2>(i,j));}

    const T& operator()(const int i,const int j) const
    {STATIC_ASSERT(d==2);return (*this)(VECTOR<int,2>(i,j));}

    T& operator()(const int i)
    {STATIC_ASSERT(d==1);return (*this)(VECTOR<int,1>(i));}

    const T& operator()(const int i) const
    {STATIC_ASSERT(d==1);return (*this)(VECTOR<int,1>(i));}

    T& operator()(const TV_INT& index)
    {assert(domain.Lazy_Inside_Half_Open(index));return array(index.Dot(stride)+offset);}

    const T& operator()(const TV_INT& index) const
    {assert(domain.Lazy_Inside_Half_Open(index));return array(index.Dot(stride)+offset);}

    bool Valid_Index(const TV_INT& index) const
    {return domain.Lazy_Inside_Half_Open(index);}

    int Standard_Index(const TV_INT& index) const
    {assert(Valid_Index(index));return index.Dot(stride)+offset;}

    void Exchange(ARRAYS_ND_BASE& a)
    {exchange(domain,a.domain);exchange(stride,a.stride);exchange(offset,a.offset);array.Exchange(a.array);}

    static void Exchange(ARRAYS_ND_BASE& a,ARRAYS_ND_BASE& b)
    {a.Exchange(b);}

    TV_INT Clamp(const TV_INT& i) const
    {return domain.Clamp(i);}

    TV_INT Clamp_End_Minus_One(const TV_INT& i) const
    {return clamp(i,domain.min_corner,domain.max_corner-1);}

    TV_INT Clamp_End_Minus_Two(const TV_INT& i) const
    {return clamp(i,domain.min_corner,domain.max_corner-2);}

    TV_INT Clamp_End_Minus_Three(const TV_INT& i) const
    {return clamp(i,domain.min_corner,domain.max_corner-3);}

    TV_INT Clamp_Interior(const TV_INT& i) const
    {return clamp(i,domain.min_corner+1,domain.max_corner-1);}

    TV_INT Clamp_Interior_End_Minus_One(const TV_INT& i) const
    {return clamp(i,domain.min_corner+1,domain.max_corner-2);}

    template<class T2>
    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const ARRAYS_ND_BASE<T2,VECTOR<int,d> >& b)
    {return a.domain==b.domain;}

    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const int m_start,const int m_end,const int n_start,const int n_end,const int mn_start,const int mn_end)
    {STATIC_ASSERT(d==3);return a.domain==RANGE<TV_INT>(m_start,m_end,n_start,n_end,mn_start,mn_end);}

    // note that these functions move the *contents* of the grid, not the grid itself, being careful about the order of grid traversal.
    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const int m_start,const int m_end,const int n_start,const int n_end)
    {STATIC_ASSERT(d==2);return a.domain==RANGE<TV_INT>(m_start,m_end,n_start,n_end);}

    // note that these functions move the *contents* of the grid, not the grid itself, being careful about the order of grid traversal.
    static bool Equal_Dimensions(const ARRAYS_ND_BASE& a,const int m_start,const int m_end)
    {STATIC_ASSERT(d==1);return a.domain==RANGE<TV_INT>(m_start,m_end);}

    static void Extract_Dimension(const ARRAYS_ND_BASE& old_array,ARRAYS_ND_BASE<typename ELEMENT_OF_VECTOR<T>::TYPE,VECTOR<int,d> >& extracted_array,int dim)
    {assert(Equal_Dimensions(old_array,extracted_array));//extracted_array.Resize(old_array.domain,false,false);
    for(int i=0;i<old_array.array.m;i++) extracted_array.array(i)=old_array.array(i)(dim);}

    static void Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy)
    {if(&old_copy!=&new_copy) Put(old_copy,new_copy);}

    static void Shifted_Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy,const TV_INT& shift)
    {Shifted_Put(old_copy,new_copy,shift);}

    static void Limited_Shifted_Get(ARRAYS_ND_BASE& new_copy,const ARRAYS_ND_BASE& old_copy,const TV_INT& shift)
    {RANGE<TV_INT> new_domain=new_copy.Domain_Indices(),old_domain=old_copy.Domain_Indices();
    RANGE<TV_INT> box(TV_INT::Componentwise_Max(new_domain.min_corner,old_domain.min_corner-shift),TV_INT::Componentwise_Min(new_domain.max_corner,old_domain.max_corner-shift));
    for(RANGE_ITERATOR<TV_INT::m> it(box);it.Valid();it.Next())
        new_copy(it.index)=old_copy(it.index+shift);}

    static void Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {if(&old_copy!=&new_copy) Put(old_copy,new_copy,RANGE<TV_INT>::Intersect(old_copy.Domain_Indices(),new_copy.Domain_Indices()));}

    void Put_With_Range(RANGE<TV_INT>& range,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {if(&old_copy!=&new_copy) Put(old_copy,new_copy,range);}

    static void Shifted_Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const TV_INT& shift)
    {if(shift==TV_INT()) Put(old_copy,new_copy);
    else for(RANGE_ITERATOR<TV_INT::m> it(new_copy.Domain_Indices());it.Valid();it.Next()) new_copy(it.index)=old_copy(it.index+shift);}

    template<class T2>
    static void Put(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy)
    {Put(constant,old_copy,new_copy,old_copy.Domain_Indices());}

protected:
    template<class T2>
    static void Put(const T2 constant,const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const RANGE<TV_INT>& box)
    {assert(old_copy.Domain_Indices().Contains(box));assert(new_copy.Domain_Indices().Contains(box));
    for(RANGE_ITERATOR<TV_INT::m> it(box);it.Valid();it.Next()) new_copy(it.index)=constant*old_copy(it.index);}

    static void Put(const ARRAYS_ND_BASE& old_copy,ARRAYS_ND_BASE& new_copy,const RANGE<TV_INT>& box)
    {assert(old_copy.Domain_Indices().Contains(box));assert(new_copy.Domain_Indices().Contains(box));
    for(RANGE_ITERATOR<TV_INT::m> it(box);it.Valid();it.Next()) new_copy(it.index)=old_copy(it.index);}
public:
    // note that these functions move the *contents* of the grid, not the grid itself, being careful about the order of grid traversal.
    void Move_Contents_By_Offset(const VECTOR<int,3>& offset)
    {STATIC_ASSERT(d==3);TV_INT i,s(offset.Componentwise_Greater_Equal(TV_INT())),c(s*2-1),e(s*domain.Edge_Lengths()),a(domain.max_corner-e),b(domain.min_corner+e);
    for(i.x=a.x;i.x<=b.x;i.x+=c.x) for(i.y=a.y;i.y<=b.y;i.y+=c.y) for(i.z=a.z;i.z<=b.z;i.z+=c.z) (*this)(i)=(*this)(i+offset);}

    void Move_Contents_By_Offset(const VECTOR<int,2>& offset)
    {STATIC_ASSERT(d==2);TV_INT i,s(offset.Componentwise_Greater_Equal(TV_INT())),c(s*2-1),e(s*domain.Edge_Lengths()),a(domain.max_corner-e),b(domain.min_corner+e);
    for(i.x=a.x;i.x<=b.x;i.x+=c.x) for(i.y=a.y;i.y<=b.y;i.y+=c.y) (*this)(i)=(*this)(i+offset);}

    void Move_Contents_By_Offset(const VECTOR<int,1>& offset)
    {STATIC_ASSERT(d==1);TV_INT i,s(offset.Componentwise_Greater_Equal(TV_INT())),c(s*2-1),e(s*domain.Edge_Lengths()),a(domain.max_corner-e),b(domain.min_corner+e);
    for(i.x=a.x;i.x<=b.x;i.x+=c.x) (*this)(i)=(*this)(i+offset);}

//#####################################################################
};
}
#endif

