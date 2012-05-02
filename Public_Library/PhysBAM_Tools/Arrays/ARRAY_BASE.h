//#####################################################################
// Copyright 2004-2009, Kevin Der, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_BASE
//#####################################################################
#ifndef __ARRAY_BASE__
#define __ARRAY_BASE__

#include <PhysBAM_Tools/Arrays/ARRAY_DIFFERENCE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_LEFT_MULTIPLE.h>
#include <PhysBAM_Tools/Arrays/ARRAY_NEGATION.h>
#include <PhysBAM_Tools/Arrays/ARRAY_PLUS_SCALAR.h>
#include <PhysBAM_Tools/Arrays/ARRAY_PRODUCT.h>
#include <PhysBAM_Tools/Arrays/ARRAY_SUM.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/maxmag.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/minmag.h>
#include <PhysBAM_Tools/Utilities/STATIC_ASSERT.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
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

private:
    template<class S> struct ELEMENT_OF{typedef typename S::ELEMENT TYPE;};
    typedef typename IF<IS_VECTOR<T>::value,ELEMENT_OF<T>,FIRST<UNUSABLE> >::TYPE::TYPE ELEMENT_OF_T;
public:

    T& operator()(const ID i)
    {return Derived()(i);}

    const T& operator()(const ID i) const
    {return Derived()(i);}

    bool Valid_Index(const ID i) const
    {return (unsigned)Value(i)<Value(Size());}

    ARRAY_VIEW<ELEMENT_OF_T> Flattened() // valid only for contiguous arrays of VECTOR<T,d>
    {T_ARRAY& self=Derived();return ARRAY_VIEW<typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    ARRAY_VIEW<const ELEMENT_OF_T> Flattened() const // valid only for contiguous arrays of VECTOR<T,d>
    {const T_ARRAY& self=Derived();return ARRAY_VIEW<const typename T::ELEMENT>(T::m*self.Size(),self.Get_Array_Pointer()->begin());}

    template<class ID2>
    ARRAY_VIEW<const T,ID2> Array_View(const ID first,const ID2 length) const
    {const T_ARRAY& self=Derived();assert((unsigned)(Value(first)+Value(length))<=(unsigned)Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<const T,ID2>(length,self.Get_Array_Pointer()+Value(first));}

    template<class ID2>
    ARRAY_VIEW<T,ID2> Array_View(const ID first,const ID2 length)
    {T_ARRAY& self=Derived();assert((unsigned)(Value(first)+Value(length))<=(unsigned)Value(self.Size()) && Value(length)>=0);return ARRAY_VIEW<T,ID2>(length,(T*)self.Get_Array_Pointer()+Value(first));}

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
    {return *this*=Inverse(a);}

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
    {assert(Size()==a2.Size());return (m*(*this*a2)).Sum().Sum();}

    template<class T2,class T_ARRAY2> typename DISABLE_IF<IS_SCALAR<T2>::value,double>::TYPE
    Inner_Product_Double_Precision(const ARRAY_BASE<T2,T_ARRAY2,ID>& m,const ARRAY_BASE<T,T_ARRAY,ID>& a2) const
    {assert(Size()==a2.Size());double result(0);ID size=Size();for(ID i(0);i<size;i++) result+=m(i).Inner_Product((*this)(i),a2(i));return result;}

    T Max() const
    {const T_ARRAY& self=Derived();T result=self(ID(0));ID m=self.Size();for(ID i(1);i<m;i++) result=PhysBAM::max(result,self(i));return result;}

    T Maxabs() const
    {const T_ARRAY& self=Derived();T result=T();;ID m=self.Size();for(ID i(0);i<m;i++) result=PhysBAM::max(result,abs(self(i)));return result;}

    T Maxmag() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result=PhysBAM::maxmag(result,self(i));return result;}

    ID Argmax() const
    {const T_ARRAY& self=Derived();ID result(0),m=self.Size();for(ID i(1);i<m;i++) if(self(i)>self(result)) result=i;return result;}

    T Min() const
    {const T_ARRAY& self=Derived();T result=self(ID(0));ID m=self.Size();for(ID i(1);i<m;i++) result=PhysBAM::min(result,self(i));return result;}

    T Minmag() const
    {const T_ARRAY& self=Derived();T result=self(ID(0));ID m=self.Size();for(ID i(1);i<m;i++) result=PhysBAM::minmag(result,self(i));return result;}

    ID Argmin() const
    {const T_ARRAY& self=Derived();ID result(0),m=self.Size();for(ID i(1);i<m;i++) if(self(i)<self(result)) result=i;return result;}

    T Componentwise_Maxabs() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result=T::Componentwise_Max(result,abs(self(i)));return result;}

    T Sum() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result+=self(i);return result;}

    double Sum_Double_Precision() const
    {const T_ARRAY& self=Derived();double result=0;ID m=self.Size();for(ID i(0);i<m;i++) result+=self(i);return result;}

    T Sumabs() const
    {const T_ARRAY& self=Derived();T result=T();ID m=self.Size();for(ID i(0);i<m;i++) result+=abs(self(i));return result;}
    
    T Average() const
    {const T_ARRAY& self=Derived();return self.Size()?Sum()/typename ARRAY_BASE<T,T_ARRAY,ID>::SCALAR(self.Size()):T();}
    
    template<class T_ARRAY1>
    ELEMENT Weighted_Sum(const T_ARRAY1& weights) const
    {STATIC_ASSERT_SAME(typename T_ARRAY1::ELEMENT,SCALAR);assert(weights.Size()==Size());
    ELEMENT result((ELEMENT()));INDEX m=Size();for(INDEX i(0);i<m;i++) result+=weights(i)*(*this)(i);return result;}

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

    void Sort()
    {::PhysBAM::Sort(*this);}

    void Fill(T value)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)=value;}

    template<class T_ARRAY1,class T_ARRAY2>
    static void Copy(const T_ARRAY1& old_copy,T_ARRAY2& new_copy)
    {new_copy=old_copy;}

    template<class T2,class T_ARRAY1,class T_ARRAY2>
    static void Copy(const T2 constant,const T_ARRAY1& array,T_ARRAY2& result)
    {result=constant*array;}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T_ARRAY2& v2,T_ARRAY3& result)
    {result=c1*v1+v2;}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T2 c2,const T_ARRAY2& v2,T_ARRAY3& result)
    {result=c1*v1+c2*v2;}

    template<class T2,class T_ARRAY1,class T_ARRAY2,class T_ARRAY3,class T_ARRAY4>
    static void Copy(const T2 c1,const T_ARRAY1& v1,const T2 c2,const T_ARRAY2& v2,const T2 c3,const T_ARRAY3& v3,T_ARRAY4& result)
    {result=c1*v1+c2*v2+c3*v3;}

    static void Get(T_ARRAY& new_copy,const T_ARRAY& old_copy)
    {if(&old_copy!=&new_copy) new_copy=old_copy.Prefix(new_copy.Size());}

    static void Put(const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {if(&old_copy!=&new_copy) new_copy.Prefix(old_copy.Size())=old_copy;}

    template<class T2>
    static void Put(const T2 constant,const T_ARRAY& old_copy,T_ARRAY& new_copy)
    {new_copy.Prefix(old_copy.Size())=constant*old_copy;}

    void Clamp_Below(const T& value)
    {T_ARRAY& self=Derived();ID m=self.Size();for(ID i(0);i<m;i++) self(i)=clamp_min(self(i),value);}

    static void Find_Common_Elements(const T_ARRAY& a,const T_ARRAY& b,T_ARRAY& result)
    {assert(&a!=&result);assert(&b!=&result);result.Remove_All();
    ID m=a.Size();for(ID i(0);i<m;i++) if(b.Contains(a(i))) result.Append(a(i));}

    template<class T_ARRAY1,class T_ARRAY2>
    static bool Equal_Dimensions(const T_ARRAY1& a,const T_ARRAY2& b)
    {return a.Size()==b.Size();}

    template<class T_ARRAY1,class T_ARRAY_INT>
    static void Permute(const T_ARRAY1& source,T_ARRAY1& destination,const T_ARRAY_INT& permutation)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(ID,typename T_ARRAY_INT::ELEMENT);
    ID m=permutation.Size();for(ID i(0);i<m;i++) destination(i)=source(permutation(i));}

    template<class T_ARRAY1,class T_ARRAY_INT>
    static void Unpermute(const T_ARRAY1& source,T_ARRAY1& destination,const T_ARRAY_INT& permutation)
    {STATIC_ASSERT_SAME(T,typename T_ARRAY1::ELEMENT);
    STATIC_ASSERT_SAME(ID,typename T_ARRAY_INT::ELEMENT);
    ID m=permutation.Size();for(ID i(0);i<m;i++) destination(permutation(i))=source(i);}

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

    T* begin() // for stl
    {return Derived().Get_Array_Pointer();}

    const T* begin() const // for stl
    {return Derived().Get_Array_Pointer();}

    T* end() // for stl
    {return Derived().Get_Array_Pointer()+Derived().Size();}

    const T* end() const // for stl
    {return Derived().Get_Array_Pointer()+Derived().Size();}

    ID Binary_Search(const T& value) const// lower_bound binary search
    {return ID(std::lower_bound(begin(),end(),value)-begin());}

//#####################################################################
};
template<class T,class T_ARRAY,class ID>
inline std::ostream& operator<<(std::ostream& output,const ARRAY_BASE<T,T_ARRAY,ID>& a)
{output<<"(";
const T_ARRAY& a_=a.Derived();
ID m=a_.Size();
for(ID i(0);i<m;i++){
    output<<a_(i);
    if(i<m-1) output<<" ";}
output<<")";
return output;}
//#####################################################################
template<class T_ARRAY1,class T_ARRAY2> struct CAN_ASSIGN<T_ARRAY1,T_ARRAY2,typename ENABLE_IF<IS_ARRAY<T_ARRAY1>::value && IS_ARRAY<T_ARRAY2>::value && IS_SAME<typename T_ARRAY1::ELEMENT,typename T_ARRAY2::ELEMENT>::value && !IS_SAME<T_ARRAY1,T_ARRAY2>::value>::TYPE>
{static const bool value=true;};
}
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#endif
