//#####################################################################
// Copyright 2008-2009, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STENCIL
//#####################################################################
#ifndef __STENCIL__
#define __STENCIL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Math_Tools/Hash.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <limits.h>
#define SUPPORT_FORMATTED_STENCIL_OUTPUT

#ifdef SUPPORT_FORMATTED_STENCIL_OUTPUT
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <iomanip>
#endif

namespace PhysBAM{

template<class T,int d> class STENCIL;
template<class T,int d> class STENCIL_ITERATOR;
template<class T,int d> inline ARRAY<T,VECTOR<int,d> >& operator+=(ARRAY<T,VECTOR<int,d> >& array,const STENCIL<T,d>& stencil);
template<class T,int d> inline ARRAY<T,VECTOR<int,d> >&
Accumulation(ARRAY<T,VECTOR<int,d> >& array,const STENCIL<T,d>& stencil,const VECTOR<int,d>& shift,const T scale);
template<class T,int d> inline ARRAY<T,VECTOR<int,d> >&
Accumulation_Cropped(ARRAY<T,VECTOR<int,d> >& array,const STENCIL<T,d>& stencil,const VECTOR<int,d>& shift,const T scale);
template<class T,int d> inline std::ostream& operator<<(std::ostream& output,const STENCIL<T,d>& stencil);

template<class T,int d>
class STENCIL
{
    typedef T ELEMENT;
    typedef ARRAY<T,VECTOR<int,d> > T_ARRAYS;

public:
    ARRAY<PAIR<VECTOR<int,d>,T> > entries;
private:

    friend class STENCIL_ITERATOR<T,d>;
    friend class STENCIL_ITERATOR<const T,d>;
    friend T_ARRAYS& operator+=<>(T_ARRAYS&,const STENCIL&);
    friend T_ARRAYS& Accumulation<>(T_ARRAYS&,const STENCIL&,const VECTOR<int,d>&,const T);
    friend T_ARRAYS& Accumulation_Cropped<>(T_ARRAYS&,const STENCIL&,const VECTOR<int,d>&,const T);
    friend std::ostream& operator<< <>(std::ostream&,const STENCIL<T,d>&);
public:

    int Size() const
    {return entries.m;}

    void Remove_All()
    {entries.Remove_All();}

    void Remove_Null_Entries()
    {for(int i=entries.m;i>=1;i--) if(entries(i).y==T()) entries.Remove_Index_Lazy(i);return;}

    void Delete(const VECTOR<int,d>& key)
    {for(int i=0;i<entries.m;i++) if(entries(i).x==key){entries.Remove_Index_Lazy(i);return;}
    PHYSBAM_FATAL_ERROR("Deleting non-existent element");}

    T& Get_Or_Insert(const VECTOR<int,d>& key)
    {for(int i=0;i<entries.m;i++) if(entries(i).x==key) return entries(i).y;
    entries.Append(PAIR<VECTOR<int,d>,T>(key,T()));return entries(entries.m).y;}

    void Insert(const VECTOR<int,d>& key,const T& data)
    {if(Contains(key)) PHYSBAM_FATAL_ERROR("Attempting to insert existing element");
    entries.Append(PAIR<VECTOR<int,d>,T>(key,data));}

    bool Contains(const VECTOR<int,d>& key) const
    {for(int i=0;i<entries.m;i++) if(entries(i).x==key) return true;return false;}

    T& Get(const VECTOR<int,d>& key)
    {for(int i=0;i<entries.m;i++) if(entries(i).x==key) return entries(i).y;
    PHYSBAM_FATAL_ERROR("Getting non-existent element");}

    const T& Get(const VECTOR<int,d>& key) const
    {for(int i=0;i<entries.m;i++) if(entries(i).x==key) return entries(i).y;
    PHYSBAM_FATAL_ERROR("Getting non-existent element");}

    const T operator()(const VECTOR<int,d>& key) const
    {for(int i=0;i<entries.m;i++) if(entries(i).x==key) return entries(i).y;
    return (T)0;}

    const T operator()(const int key) const
    {STATIC_ASSERT(d==1);
    for(int i=0;i<entries.m;i++) if(entries(i).x(1)==key) return entries(i).y;
    return (T)0;}

    T operator*(const T_ARRAYS& array) const
    {T result=0;for(int i=0;i<entries.m;i++) result+=array(entries(i).x)*entries(i).y;return result;}

    T operator*(const ARRAY<T>& array) const
    {STATIC_ASSERT(d==1);
    T result=0;for(int i=0;i<entries.m;i++) result+=array(entries(i).x(1))*entries(i).y;return result;}

    T Contraction(const T_ARRAYS& array,const VECTOR<int,d>& shift) const
    {T result=0;for(int i=0;i<entries.m;i++) result+=array(entries(i).x+shift)*entries(i).y;return result;}

    T Contraction_Cropped(const T_ARRAYS& array,const VECTOR<int,d>& shift) const
    {T result=0;for(int i=0;i<entries.m;i++) if(array.Domain_Indices().Lazy_Inside(entries(i).x+shift)) result+=array(entries(i).x+shift)*entries(i).y;return result;}

    T Contraction_Scaled(const T_ARRAYS& array,const VECTOR<int,d>& shift,const T_ARRAYS& scale_array) const
    {T result=0;for(int i=0;i<entries.m;i++) result+=array(entries(i).x+shift)*entries(i).y*scale_array(entries(i).x+shift);return result;}

    T operator*(const STENCIL& stencil) const
    {T result=0;for(int i=0;i<entries.m;i++) result+=stencil(entries(i).x)*entries(i).y;return result;}

    STENCIL operator+(const STENCIL& stencil) const
    {STENCIL result(*this);for(int i=0;i<stencil.entries.m;i++) result.Get_Or_Insert(stencil.entries(i).x)+=stencil.entries(i).y;return result;}

    STENCIL operator*(const T a) const
    {STENCIL result(*this);for(int i=0;i<entries.m;i++) result.entries(i).y*=a;return result;}

    STENCIL Convolve(const STENCIL& stencil) const
    {STENCIL result;
    for(STENCIL_ITERATOR<const T,d> iterator1(*this);iterator1.Valid();iterator1.Next())
        for(STENCIL_ITERATOR<const T,d> iterator2(stencil);iterator2.Valid();iterator2.Next())
            result.Get_Or_Insert(iterator1.Key()+iterator2.Key())+=iterator1.Data()*iterator2.Data();
    return result;}

    STENCIL Elementwise_Multiply(const T_ARRAYS& array) const
    {STENCIL result(*this);for(int i=0;i<entries.m;i++) result.entries(i).y*=array(entries(i).x);return result;}

    STENCIL& operator+=(const STENCIL& stencil)
    {for(int i=0;i<stencil.entries.m;i++) Get_Or_Insert(stencil.entries(i).x)+=stencil.entries(i).y;return *this;}

    STENCIL& operator=(const STENCIL& stencil)
    {entries.Remove_All();for(int i=0;i<stencil.entries.m;i++) Get_Or_Insert(stencil.entries(i).x)=stencil.entries(i).y;return *this;}
    
    bool operator==(const STENCIL& stencil) const // order and null entries important
    {return entries==stencil.entries;}
    
    bool operator!=(const STENCIL& stencil) const // order and null entries important
    {return entries!=stencil.entries;}
    
    STENCIL Shift(const VECTOR<int,d>& shift) const
    {STENCIL result(*this);for(int i=0;i<entries.m;i++) result.entries(i).x+=shift;return result;}

    template<class ExistsFunctor>
    STENCIL Crop(ExistsFunctor exists)
    {STENCIL result(*this);for(int i=entries.m;i>=1;i--) if(!exists(entries(i).x)) result.entries.Remove_Index(i);return result;}

    T Magnitude_Squared() const
    {T result=0;for(int i=0;i<entries.m;i++) result+=sqr(entries(i).y);return result;}
    
    void Normalize()
    {T result=0;for(int i=0;i<entries.m;i++) result+=entries(i).y;
    if (result==0) return;
    for(int i=0;i<entries.m;i++) entries(i).y/=result;}
    
    void Print_Stencil(std::ostream& output) const
    {for(int i=0;i<entries.m;i++) output<<entries(i).x<<"->"<<entries(i).y<<";";output<<std::endl;}
    
    STENCIL Clamped_Shift(const VECTOR<int,d>& shift,const int n) const
    {STENCIL result;
    for(int i=0;i<entries.m;i++){
        PAIR<VECTOR<int,d>,T> new_entry(entries(i).x+shift,entries(i).y);
        VECTOR<int,2>& index=new_entry.x;
        index.x=(index.x+n-1)%n+1;
        index.y=(index.y+n-1)%n+1;
        result.entries.Append(new_entry);}
    return result;}

    void Clamp_Shift(const VECTOR<int,d>& shift,const int n)
    {for(int i=0;i<entries.m;i++){
        VECTOR<int,2>& index=entries(i).x;
        index.x=(index.x+shift.x+n-1)%n+1;
        index.y=(index.y+shift.y+n-1)%n+1;}}


    friend HASH Hash_Reduce(const STENCIL& key){return HASH(key.entries);}

//#####################################################################
};
template<class T,int d,int d2> inline T 
Contraction(const VECTOR<STENCIL<T,d>,d2>& v1,const VECTOR<ARRAY<T,VECTOR<int,d> >,d2>& v2,const VECTOR<int,d>& shift)
{T result=0;for(int i=0;i<d2;i++) result+=v1(i).Contraction(v2(i),shift);return result;}

template<class T,int d,int d2> inline T 
Contraction_Cropped(const VECTOR<STENCIL<T,d>,d2>& v1,const VECTOR<ARRAY<T,VECTOR<int,d> >,d2>& v2,const VECTOR<int,d>& shift)
{T result=0;for(int i=0;i<d2;i++) result+=v1(i).Contraction_Cropped(v2(i),shift);return result;}

template<class T,int d,int d2> inline T 
Dot_Product(const VECTOR<STENCIL<T,d>,d2>& v1,const VECTOR<ARRAY<T,VECTOR<int,d> >,d2>& v2)
{T result=0;for(int i=0;i<d2;i++) result+=v1(i)*v2(i);return result;}

template<class T,int d,int d2> inline T 
Dot_Product(const VECTOR<STENCIL<T,d>,d2>& v1,const VECTOR<STENCIL<T,d>,d2>& v2)
{T result=0;for(int i=0;i<d2;i++) result+=v1(i)*v2(i);return result;}

template<class T,int d,int d2> inline VECTOR<STENCIL<T,d>,d2>
operator*(const VECTOR<STENCIL<T,d>,d2>& v,const T a)
{VECTOR<STENCIL<T,d>,d2> result;for(int i=0;i<d2;i++) result(i)=v(i)*a;return result;}

template<class T,int d> inline ARRAY<T,VECTOR<int,d> >&
operator+=(ARRAY<T,VECTOR<int,d> >& array,const STENCIL<T,d>& stencil)
{for(int i=0;i<stencil.entries.m;i++) array(stencil.entries(i).x)+=stencil.entries(i).y;return array;}

template<class T,int d,int d2> inline VECTOR<ARRAY<T,VECTOR<int,d> >,d2>&
operator+=(VECTOR<ARRAY<T,VECTOR<int,d> >,d2>& v1,const VECTOR<STENCIL<T,d>,d2>& v2)
{for(int i=0;i<d2;i++) v1(i)+=v2(i);return v1;}

template<class T,int d> inline ARRAY<T,VECTOR<int,d> >&
Accumulation_Cropped(ARRAY<T,VECTOR<int,d> >& array,const STENCIL<T,d>& stencil,const VECTOR<int,d>& shift,const T scale)
{for(int i=0;i<stencil.entries.m;i++) if(array.Domain_Indices().Lazy_Inside(stencil.entries(i).x+shift)) array(stencil.entries(i).x+shift)+=stencil.entries(i).y*scale;return array;}

template<class T,int d,class ExistsFunctor> inline ARRAY<T,VECTOR<int,d> >&
Accumulation_Cropped(ARRAY<T,VECTOR<int,d> >& array,const STENCIL<T,d>& stencil,const VECTOR<int,d>& shift,const T scale,ExistsFunctor exists)
{for(int i=0;i<stencil.entries.m;i++) if(exists(stencil.entries(i).x+shift)) array(stencil.entries(i).x+shift)+=stencil.entries(i).y*scale;return array;}

template<class T,int d,int d2> inline VECTOR<ARRAY<T,VECTOR<int,d> >,d2>&
Accumulation_Cropped(VECTOR<ARRAY<T,VECTOR<int,d> >,d2>& v1,const VECTOR<STENCIL<T,d>,d2>& v2,const VECTOR<int,d>& shift,const T scale)
{for(int i=0;i<d2;i++) Accumulation_Cropped(v1(i),v2(i),shift,scale);return v1;}

template<class T,int d> inline ARRAY<T,VECTOR<int,d> >&
Accumulation(ARRAY<T,VECTOR<int,d> >& array,const STENCIL<T,d>& stencil,const VECTOR<int,d>& shift,const T scale)
{for(int i=0;i<stencil.entries.m;i++) array(stencil.entries(i).x+shift)+=stencil.entries(i).y*scale;return array;}

template<class T,int d,int d2> inline VECTOR<ARRAY<T,VECTOR<int,d> >,d2>&
Accumulation(VECTOR<ARRAY<T,VECTOR<int,d> >,d2>& v1,const VECTOR<STENCIL<T,d>,d2>& v2,const VECTOR<int,d>& shift,const T scale)
{for(int i=0;i<d2;i++) Accumulation(v1(i),v2(i),shift,scale);return v1;}

template<class T,int d,int d2> inline T
Magnitude_Squared(const VECTOR<STENCIL<T,d>,d2>& v)
{T result=0;for(int i=0;i<d2;i++) result+=v(i).Magnitude_Squared();return result;}

template<class T,int d> inline std::ostream&
operator<<(std::ostream& output,const STENCIL<T,d>& stencil)
{for(int i=0;i<stencil.entries.m;i++) output<<stencil.entries(i).x<<"->"<<stencil.entries(i).y<<";";return output;}

#ifdef SUPPORT_FORMATTED_STENCIL_OUTPUT

template<class T> void 
Print_Stencil(std::ostream& output,const STENCIL<T,2>& stencil,std::string prefix="DEBUG:",int imin=INT_MAX,int imax=INT_MIN,int jmin=INT_MAX,int jmax=INT_MIN)
{const int index_width=3,data_width=8,max_width=max(index_width,data_width);
if(!stencil.Size()){output<<prefix<<"<empty stencil>"<<std::endl;return;}
for(STENCIL_ITERATOR<const T,2> iterator(stencil);iterator.Valid();iterator.Next()){
    int i,j;iterator.Key().Get(i,j);imin=min(i,imin);imax=max(i,imax);jmin=min(j,jmin);jmax=max(j,jmax);}
output<<prefix;for(int i=0;i<index_width;i++) output<<" ";
for(int i=imin;i<=imax;i++) output<<" | "<<std::setw(max_width)<<i;
output<<std::endl;
output<<prefix;for(int i=1;i<=(imax-imin+1)*(max_width+3)+index_width;i++) output<<"-";
output<<std::endl;
for(int j=jmax;j>=jmin;j--){
    output<<prefix<<std::setw(index_width)<<j;
    for(int i=imin;i<=imax;i++){
        output<<" | ";
        if(stencil.Contains(VECTOR<int,2>(i,j))) output<<std::setw(max_width)<<stencil.Get(VECTOR<int,2>(i,j));
        else for(int k=0;k<max_width;k++) output<<" ";}
    output<<std::endl;}}

template<class T> void 
Print_Stencil(std::ostream& output,const STENCIL<T,3>& stencil)
{if(!stencil.Size()){output<<"<empty stencil>"<<std::endl;return;}
HASHTABLE<int,STENCIL<T,2> > slices;
RANGE<VECTOR<int,3> > box(INT_MAX,INT_MIN,INT_MAX,INT_MIN,INT_MAX,INT_MIN);
for(STENCIL_ITERATOR<const T,3> iterator(stencil);iterator.Valid();iterator.Next()){
    const VECTOR<int,3>& index=iterator.Key();
    box.Enlarge_To_Include_Point(index);
    slices.Get_Or_Insert(index.z).Insert(index.Remove_Index(3),iterator.Data());}
for(int z=box.max_corner.z;z>=box.min_corner.z;z--){
    output<<"z="<<z<<std::endl;
    Print_Stencil(output,slices.Get_Default(z),"  ",
        box.min_corner.x,box.max_corner.x,
        box.min_corner.y,box.max_corner.y
);}}
#endif

}
#endif
