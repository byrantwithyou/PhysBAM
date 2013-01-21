//#####################################################################
// Copyright 2007, Jon Gretarsson, Avi Robinson-Mosher, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_FLAT_MXN
//#####################################################################
#ifndef __SPARSE_MATRIX_FLAT_MXN__
#define __SPARSE_MATRIX_FLAT_MXN__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
namespace PhysBAM{

template<class T>
class SPARSE_MATRIX_FLAT_MXN
{
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    typedef T SCALAR;
    int m,n; // size of the m by n matrix
    ARRAY<int> offsets;
    ARRAY<SPARSE_MATRIX_ENTRY<T> > A;
    SPARSE_MATRIX_FLAT_MXN<T>* Q;
    SPARSE_MATRIX_FLAT_MXN<T>* L;

    SPARSE_MATRIX_FLAT_MXN();

    SPARSE_MATRIX_FLAT_MXN(const SPARSE_MATRIX_FLAT_MXN<T>& matrix);

    ~SPARSE_MATRIX_FLAT_MXN();

    SPARSE_MATRIX_FLAT_MXN& operator=(const SPARSE_MATRIX_FLAT_MXN& matrix)
    {m=matrix.m;n=matrix.n;offsets=matrix.offsets;A=matrix.A;
    return *this;}

    SPARSE_MATRIX_FLAT_MXN& operator=(const SPARSE_MATRIX_FLAT_NXN<T>& matrix)
    {m=matrix.n;n=matrix.n;offsets=matrix.offsets;A=matrix.A;
    return *this;}

    const T operator()(const int i,const int j) const
    {int index=Find_Index(i,j);assert(A(index).j==j);return A(index).a;}

    T& operator()(const int i,const int j);

    void Set_Element(const int i,const int j,const T element)
    {(*this)(i,j)=element;}

    void Add_Element(const int i,const int j,const T element)
    {(*this)(i,j)+=element;}

    template<class RW> void Read(std::istream& input)
    {Read_Binary<RW>(input,m,n,offsets,A);}

    template<class RW> void Write(std::ostream& output)
    {Write_Binary<RW>(output,m,n,offsets,A);}

//#####################################################################
    SPARSE_MATRIX_FLAT_NXN<T>* Create_Submatrix(const INTERVAL<int>& rows);
    void Set_Row_Lengths(ARRAY_VIEW<int> lengths);
    int Find_Index(const int i,const int j) const;
    int Find_Index_Exists(const int i,const int j) const;
    bool Element_Present(const int i,const int j) const;
    void Times(const ARRAY<T>& x,ARRAY<T>& result) const;
    void Transpose_Times(const ARRAY<T>& x,ARRAY<T>& result) const;
    void Times_Add(const ARRAY<T>& x,ARRAY<T>& result) const;
    void Times_Add_Row(const ARRAY<T>& x,ARRAY<T>& result,const int row) const;
    void Times_Subtract(const ARRAY<T>& x,ARRAY<T>& result) const;
    void Transpose_Times_Add(const ARRAY<T>& x,ARRAY<T>& result) const;
    void Transpose_Times_Subtract(const ARRAY<T>& x,ARRAY<T>& result) const;
    void Negate();
    SPARSE_MATRIX_FLAT_MXN<T>& operator*=(const T a);
    SPARSE_MATRIX_FLAT_MXN<T>& operator/=(const T a);
    SPARSE_MATRIX_FLAT_MXN<T>& operator+=(const T a);
    SPARSE_MATRIX_FLAT_MXN<T>& operator-=(const T a);
    void Compress(SPARSE_MATRIX_FLAT_MXN<T>& compressed);
    void Transpose(SPARSE_MATRIX_FLAT_MXN<T>& A_transpose) const;
    SPARSE_MATRIX_FLAT_MXN<T> Times_Transpose(const SPARSE_MATRIX_FLAT_MXN<T>& rhs);
    SPARSE_MATRIX_FLAT_MXN<T> Times_Diagonal_Times(const ARRAY<T> diagonal,const SPARSE_MATRIX_FLAT_MXN<T>& rhs); // (*this) * diagonal * (rhs)
    SPARSE_MATRIX_FLAT_MXN<T> Scale_Rows(const ARRAY<T>& d) const;
    SPARSE_MATRIX_FLAT_NXN<T> operator+(const SPARSE_MATRIX_FLAT_NXN<T>& A_rhs) const;
    SPARSE_MATRIX_FLAT_MXN<T> operator+(const SPARSE_MATRIX_FLAT_MXN<T>& A_rhs) const;
    SPARSE_MATRIX_FLAT_MXN<T> operator-(const SPARSE_MATRIX_FLAT_MXN<T>& A_rhs) const;
    SPARSE_MATRIX_FLAT_MXN<T> operator*(const SPARSE_MATRIX_FLAT_MXN<T>& rhs) const;
    SPARSE_MATRIX_FLAT_NXN<T> Create_NXN_Matrix();
    void Set_Times_Diagonal(const ARRAY<T>& D);
    void Set_Diagonal_Times(const ARRAY<T>& D);
    void Write_Row_Lengths();
    void Print_Row(const int row);
    void Reset(const int c);
    void Append_Entry_To_Current_Row(const int c,const T a);
    void Finish_Row();
    void Sort_Entries();
    void Get_Row(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q_i,int row);
    void Construct_Incomplete_LQ_Factorization(const int p_l=10,const int p_q=10,const T zero_tolerance=1e-8,const T zero_replacement=1e-8);
    void Fast_Sparse_Multiply(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q,ARRAY<SPARSE_MATRIX_ENTRY<T> >& l);
//#####################################################################
};
template<class T> std::ostream& operator<<(std::ostream& output_stream,const SPARSE_MATRIX_FLAT_MXN<T>& A);
}
#endif
