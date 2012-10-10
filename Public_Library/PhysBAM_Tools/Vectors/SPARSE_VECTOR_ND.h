//#####################################################################
// Copyright 2003-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_VECTOR_ND
//#####################################################################
#ifndef __SPARSE_VECTOR_ND__
#define __SPARSE_VECTOR_ND__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
namespace PhysBAM{

template<class T>
class SPARSE_VECTOR_ND
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    int n;
    ARRAY<int> indices;
    ARRAY<T> x;

     SPARSE_VECTOR_ND(const int n_input)
        :n(n_input)
    {
    }

    SPARSE_VECTOR_ND(const SPARSE_VECTOR_ND& vector)
        :n(vector.n),indices(vector.indices),x(vector.x)
    {
    }

    ~SPARSE_VECTOR_ND() {}

private:
    void operator=(const SPARSE_VECTOR_ND&);
public:

    bool Element_Present_And_Location(const int i,int& location)
    {assert((unsigned)i<(unsigned)n);
    location=indices.Binary_Search(i);
    return location<indices.m && indices(location)==i;}

    const T operator()(const int i) const
    {int k;if(Element_Present_And_Location(i,k)) return x(k);return T();}

    int Get_Or_Add_Index(const int i)
    {int k;if(!Element_Present_And_Location(i,k)){indices.Insert(i,k);x.Insert(T(),k);}return k;}

    void Set_Element(const int i,const T& element)
    {int k=Get_Or_Add_Index(i);x(k)=element;}

    void Add_Element(const int i,const T& element)
    {int k=Get_Or_Add_Index(i);x(k)+=element;}

    bool Element_Present(const int i)
    {int k;return Element_Present_And_Location(i,k);}

    void Clear()
    {indices.Clean_Memory();x.Clean_Memory();}

    T Dot_Product(const ARRAY<T>& vector)
    {T sum=T();for(int i=0;i<indices.m;i++) sum+=x(i)*vector(indices(i));return sum;}

    void Negate()
    {x=-x;}

    SPARSE_VECTOR_ND<T>& operator*=(const T a)
    {x*=a;}

    void Write_Internal_Arrays(std::ostream& output_stream)
    {for(int i=0;i<indices.m;i++) output_stream<<indices(i)<<", "<<x(i)<<std::endl;}

    // template<class RW> void Read(std::istream& input)
    // {Read_Binary<RW>(input,n,number_of_active_indices);
    // delete[] indices;delete[] x;indices=new int[number_of_active_indices+1];x=new T[number_of_active_indices+1];
    // for(int i=0;i<number_of_active_indices;i++){Read_Binary<RW>(input,indices[i]);Read_Binary<RW>(input,x[i]);}}

    // template<class RW> void Write(std::ostream& output) const
    // {Write_Binary<RW>(output,n,number_of_active_indices);
    // for(int i=0;i<number_of_active_indices;i++){Write_Binary<RW>(output,indices[i]);Write_Binary<RW>(output,x[i]);}}
//#####################################################################
};
template<class T> inline std::ostream& operator<<(std::ostream& output_stream,const SPARSE_VECTOR_ND<T>& v)
{for(int i=0;i<v.n;i++)output_stream<<v(i)<<" ";output_stream<<std::endl;return output_stream;}
}
#endif
