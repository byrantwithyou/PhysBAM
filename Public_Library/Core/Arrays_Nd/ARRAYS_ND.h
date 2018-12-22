//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Duc Nguyen, Craig Schroeder, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY
//#####################################################################
#ifndef __ARRAYS_ND__
#define __ARRAYS_ND__

#include <Core/Arrays_Nd/ARRAYS_ND_BASE.h>
#include <Core/Math_Tools/FIXED_NUMBER.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Utilities/EXCEPTIONS.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int d>
class ARRAY<T,VECTOR<int,d> >:public ARRAYS_ND_BASE<T,VECTOR<int,d> >
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef VECTOR<int,d> TV_INT;
    enum WORKAROUND1 {dimension=d};
    template<class T2> struct REBIND{typedef ARRAY<T2,TV_INT> TYPE;};
    typedef T ELEMENT;typedef TV_INT INDEX;typedef T& RESULT_TYPE;typedef const T& CONST_RESULT_TYPE;
    typedef ARRAYS_ND_BASE<T,VECTOR<int,d> > BASE;

    using BASE::array; // one-dimensional data storage
    using BASE::domain;using BASE::Exchange;using BASE::stride;using BASE::offset;
    using BASE::Calculate_Acceleration_Constants;

    ARRAY() = default;

    ARRAY(const RANGE<TV_INT>& domain_input)
    {Initialize(domain_input,true,T());}

    ARRAY(const RANGE<TV_INT>& domain_input,NO_INIT)
    {Initialize(domain_input,false,T());}

    ARRAY(const RANGE<TV_INT>& domain_input,USE_INIT,const T& initialization_value=T())
    {Initialize(domain_input,true,initialization_value);}

    ARRAY(const TV_INT& size_input)
        :ARRAY(RANGE<TV_INT>(TV_INT(),size_input))
    {}

    ARRAY(const TV_INT& size_input,NO_INIT)
        :ARRAY(RANGE<TV_INT>(TV_INT(),size_input),no_init)
    {}

    ARRAY(const TV_INT& size_input,USE_INIT,const T& initialization_value=T())
        :ARRAY(RANGE<TV_INT>(TV_INT(),size_input),use_init,initialization_value)
    {}

    ARRAY(const ARRAY& old_array)
    {Initialize(old_array.domain,false,T());array=old_array.array;}

    ARRAY(const ARRAY& old_array,NO_INIT)
    {Initialize(old_array.domain,false,T());}

    ARRAY(ARRAY&&) = default;

    ~ARRAY()
    {delete[] array.Get_Array_Pointer();}

protected:
    void Initialize(const RANGE<TV_INT>& box,const bool initialize_using_initialization_value,const T& initialization_value)
    {TV_INT counts_new=box.Edge_Lengths();
    assert(counts_new.Min()>=0);
    int size_new=counts_new.Product();
    Calculate_Acceleration_Constants(box);
    array.Set(size_new,new T[size_new]);
    if(initialize_using_initialization_value) array.Fill(initialization_value);}
public:

    void Clean_Memory()
    {Resize(RANGE<TV_INT>::Empty_Box(),no_init);}

    void Delete_Pointers_And_Clean_Memory() // only valid if T is a pointer type
    {for(int i=0;i<array.Size();i++) delete array(i);Clean_Memory();}

    ARRAY& operator=(const ARRAY& source)
    {if(this==&source) return *this;
    Resize_In_Place(source.Domain_Indices());
    array=source.array;return *this;}

    ARRAY& operator=(ARRAY&& a)
    {domain=a.domain;stride=a.stride;offset=a.offset;array=std::move(a.array);return *this;}

    template<class T_ARRAY1>
    ARRAY& operator=(const T_ARRAY1& source)
    {STATIC_ASSERT(is_same<ELEMENT,typename T_ARRAY1::ELEMENT>::value);
    Resize_In_Place(source.Domain_Indices());
    ARRAY_BASE<T,BASE,TV_INT>::operator=(source);return *this;}

    void Resize(const RANGE<TV_INT>& box)
    {if(box==domain) return;
    ARRAY new_array(box);
    this->Put(*this,new_array,domain.Intersect(box));
    Exchange(new_array);}

    void Resize(const RANGE<TV_INT>& box,NO_INIT)
    {if(box==domain) return;
    ARRAY new_array(box,no_init);
    Exchange(new_array);}

    void Resize(const RANGE<TV_INT>& box,USE_INIT,const T& initialization_value)
    {if(box==domain) return;
    ARRAY new_array(box,use_init,initialization_value);
    this->Put(*this,new_array,domain.Intersect(box));
    Exchange(new_array);}

    void Resize(const RANGE<TV_INT>& box,INIT_ALL,const T& initialization_value)
    {Resize(box,no_init);array.Fill(initialization_value);}

    void Resize(const TV_INT& corner)
    {Resize(RANGE<TV_INT>(TV_INT(),corner));}

    void Resize(const TV_INT& corner,NO_INIT)
    {Resize(RANGE<TV_INT>(TV_INT(),corner),no_init);}

    void Resize(const TV_INT& corner,USE_INIT,const T& initialization_value)
    {Resize(RANGE<TV_INT>(TV_INT(),corner),use_init,initialization_value);}

    void Resize(const TV_INT& corner,INIT_ALL,const T& initialization_value)
    {Resize(RANGE<TV_INT>(TV_INT(),corner),init_all,initialization_value);}

    void Reallocate_In_Place(const RANGE<TV_INT>& box)
    {TV_INT counts_new(box.Edge_Lengths());int size_new=counts_new.Product();Calculate_Acceleration_Constants(box);delete [] array.Get_Array_Pointer();array.Set(size_new,new T[size_new]);}

    void Resize_In_Place(const RANGE<TV_INT>& box)
    {TV_INT counts_new(box.Edge_Lengths());
    if(array.Size()>=counts_new.Product()) Calculate_Acceleration_Constants(box);
    else Reallocate_In_Place(box);}

    template<class RW>
    void Read(std::istream& input)
    {Read_With_Length<RW>(input,1);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_With_Length<RW>(output,1);}

protected:
    template<class RW>
    void Read_With_Length(std::istream& input,const int length2)
    {int read_length;RANGE<TV_INT> domain_temp;Read_Binary<RW>(input,read_length,domain_temp);
    Resize(domain_temp,no_init);
    Read_Binary_Array<RW>(input,array.Get_Array_Pointer(),array.Size());}

    template<class RW>
    void Write_With_Length(std::ostream& output,const int length2) const
    {Write_Binary<RW>(output,length2,domain);Write_Binary_Array<RW>(output,array.Get_Array_Pointer(),array.Size());}
//#####################################################################
};

template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const ARRAYS_ND_BASE<T,VECTOR<int,1> >& a)
{
    for(int i=a.domain.min_corner.x;i<a.domain.max_corner.x;i++) output_stream<<a(i)<<" ";
    output_stream<<std::endl;
    return output_stream;
}

template<class T> inline std::ostream& operator<<(std::ostream& output,const ARRAYS_ND_BASE<T,VECTOR<int,2> >& a)
{
    for(int i=a.domain.min_corner.x;i<a.domain.max_corner.x;i++){
        for(int j=a.domain.min_corner.y;j<a.domain.max_corner.y;j++)
            output<<a(i,j)<<" ";
        output<<std::endl;}
    return output;
}

template<class T> inline std::ostream& operator<<(std::ostream& output,const ARRAYS_ND_BASE<T,VECTOR<int,3> >& a)
{
    for(int i=a.domain.min_corner.x;i<a.domain.max_corner.x;i++){
        for(int j=a.domain.min_corner.y;j<a.domain.max_corner.y;j++){
            for(int k=a.domain.min_corner.z;k<a.domain.max_corner.z;k++)
                output<<a(i,j,k)<<" ";
            output<<std::endl;}
    output<<"------------------------------------------"<<std::endl;}
    return output;
}
}
#endif
