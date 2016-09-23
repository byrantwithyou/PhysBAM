//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STATIC_TENSOR
//#####################################################################
#ifndef __STATIC_TENSOR__
#define __STATIC_TENSOR__

#include <Core/Utilities/STATIC_ASSERT.h>
#include <Core/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int rank,int dim>
class STATIC_TENSOR
{
    struct UNUSABLE{};
public:
    typedef T SCALAR;
    typedef T ELEMENT;
    typedef STATIC_TENSOR<T,rank-1,dim> LOWER;
    typedef VECTOR<int,rank> TV_INT;
    typedef TV_INT INDEX;

    VECTOR<LOWER,dim> x;

    T& operator()(const int* index)
    {assert((unsigned)*index<(unsigned)dim);return x[*index](index+1);}

    const T& operator()(const int* index) const
    {assert((unsigned)*index<(unsigned)dim);return x[*index](index+1);}

    T& operator()(const VECTOR<int,rank>& index)
    {return (*this)(&index(0));}

    const T& operator()(const VECTOR<int,rank>& index) const
    {return (*this)(&index(0));}

    template<int dim2>
    void Set_Subset(const STATIC_TENSOR<T,rank,dim2>& st)
    {for(int i=0;i<dim2;i++) x[i].Set_Subset(st.x[i]);}

    template<int dim2>
    void Add_Subset(const STATIC_TENSOR<T,rank,dim2>& st)
    {for(int i=0;i<dim2;i++) x[i].Add_Subset(st.x[i]);}

    template<int dim2>
    void Subtract_Subset(const STATIC_TENSOR<T,rank,dim2>& st)
    {for(int i=0;i<dim2;i++) x[i].Subtract_Subset(st.x[i]);}
};

template<class T,int dim>
class STATIC_TENSOR<T,0,dim>
{
    struct UNUSABLE{};
public:
    typedef T SCALAR;
    typedef T ELEMENT;
    enum WORKAROUND {rank=0};
    typedef VECTOR<int,rank> TV_INT;
    typedef TV_INT INDEX;

    T x;

    STATIC_TENSOR()
        :x(T())
    {}

    T& operator()(const int* index)
    {return x;}

    const T& operator()(const int* index) const
    {return x;}

    T& operator()(const VECTOR<int,rank>& index)
    {return x;}

    const T& operator()(const VECTOR<int,rank>& index) const
    {return x;}

    template<int dim2>
    void Set_Subset(const STATIC_TENSOR<T,0,dim2>& st)
    {x=st.x;}

    template<int dim2>
    void Add_Subset(const STATIC_TENSOR<T,0,dim2>& st)
    {x+=st.x;}

    template<int dim2>
    void Subtract_Subset(const STATIC_TENSOR<T,0,dim2>& st)
    {x-=st.x;}
};
}
#endif
