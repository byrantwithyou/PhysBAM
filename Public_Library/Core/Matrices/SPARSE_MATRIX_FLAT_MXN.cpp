//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Frederic Gibou, Jon Gretarsson, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPARSE_MATRIX_FLAT_MXN
//#####################################################################
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/Inverse.h>
#include <Core/Math_Tools/sqr.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>::
SPARSE_MATRIX_FLAT_MXN()
    :m(0),n(0),Q(0),L(0),C(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>::
SPARSE_MATRIX_FLAT_MXN(const SPARSE_MATRIX_FLAT_MXN<T>& matrix)
    :m(matrix.m),n(matrix.n),offsets(matrix.offsets),A(matrix.A),Q(0),L(0),C(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>::
~SPARSE_MATRIX_FLAT_MXN()
{
    delete Q;
    delete L;
    delete C;
}
//#####################################################################
// Function Create_Submatrix
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>* SPARSE_MATRIX_FLAT_MXN<T>::
Create_Submatrix(const INTERVAL<int>& rows)
{
    assert(rows.Size()==m);
    int entries=0;for(int index=offsets(0);index<offsets(m);index++)if(rows.Lazy_Inside_Half_Open(A(index).j)) entries++;
    SPARSE_MATRIX_FLAT_MXN<T>* submatrix=new SPARSE_MATRIX_FLAT_MXN<T>();
    submatrix->n=rows.Size();
    submatrix->offsets.Resize(submatrix->n+1);
    submatrix->A.Resize(entries);
    int next_index=0;
    for(int i=0;i<submatrix->n;i++){
        submatrix->offsets(i)=next_index;
        for(int old_index=offsets(i);old_index<offsets(i+1);old_index++)if(rows.Lazy_Inside_Half_Open(A(old_index).j)){
            submatrix->A(next_index).j=A(old_index).j-rows.min_corner;submatrix->A(next_index).a=A(old_index).a;next_index++;}}
    submatrix->offsets(submatrix->n)=next_index;
    return submatrix;
}
//#####################################################################
// Function Set_Row_Lengths
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Set_Row_Lengths(ARRAY_VIEW<int> lengths)
{
    Reset(lengths.m);
    m=lengths.m;offsets.Resize(m+1,false,false);offsets(0)=0;
    for(int i=0;i<m;i++){offsets(i+1)=offsets(i)+lengths(i);}
    A.Resize(offsets(m));
}
//#####################################################################
// Function Find_Index
//#####################################################################
template<class T> int SPARSE_MATRIX_FLAT_MXN<T>::
Find_Index(const int i,const int j) const
{
    assert(A.m);assert((unsigned)i<(unsigned)m);assert((unsigned)j<(unsigned)n);
    int index=offsets(i);while(A(index).j>=0 && A(index).j<j)index++;
    assert(index<offsets(i+1));return index;
}
//#####################################################################
// Function Find_Index_Return
//#####################################################################
template<class T> int SPARSE_MATRIX_FLAT_MXN<T>::
Find_Index_Exists(const int i,const int j) const
{
    assert(A.m);assert((unsigned)i<(unsigned)m);assert((unsigned)j<(unsigned)n);
    int index=offsets(i);while(A(index).j>=0 && A(index).j<j)index++;
    if(index<offsets(i+1) && A(index).j==j)
        return index;
    else
        return -1;
}
//#####################################################################
// Function operator()
//#####################################################################
template<class T> T& SPARSE_MATRIX_FLAT_MXN<T>::
operator()(const int i,const int j)
{
    int index=Find_Index(i,j);
    if(A(index).j!=j){ // need to add entry
        if(A(index).j>=0){ // shift entries over to make room
            assert(A(offsets(i+1)-1).j<0);
            for(int jj=offsets(i+1)-1;jj>index;jj--) A(jj)=A(jj-1);}
        A(index).j=j;A(index).a=0;}
    return A(index).a;
}
//#####################################################################
// Function Element_Present
//#####################################################################
template<class T> bool SPARSE_MATRIX_FLAT_MXN<T>::
Element_Present(const int i,const int j) const
{
    assert((unsigned)i<(unsigned)m);assert((unsigned)j<(unsigned)n);
    for(int index=offsets(i);index<offsets(i+1);index++)if(A(index).j==j) return true;
    return false;
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Times_Add(ARRAY_VIEW<const T> x,ARRAY_VIEW<T> result) const
{
    int index=offsets(0);
    for(int i=0;i<m;i++){
        int end=offsets(i+1);T sum=(T)0;
        for(;index<end;index++) sum+=A(index).a*x(A(index).j);
        result(i)+=sum;}
}
//#####################################################################
// Function Times_Add
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Times_Add_Row(ARRAY_VIEW<const T> x,ARRAY_VIEW<T> result,const int row) const
{
    int index=offsets(row);
    int end=offsets(row+1);
    T sum=(T)0;
    for(;index<end;index++) sum+=A(index).a*x(A(index).j);
    result(row)+=sum;
}
//#####################################################################
// Function Times_Subtract
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Times_Subtract(ARRAY_VIEW<const T> x,ARRAY_VIEW<T> result) const
{
    int index=offsets(0);
    for(int i=0;i<m;i++){
        int end=offsets(i+1);T sum=(T)0;
        for(;index<end;index++) sum+=A(index).a*x(A(index).j);
        result(i)-=sum;}
}
//#####################################################################
// Function Times
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Times(ARRAY_VIEW<const T> x,ARRAY_VIEW<T> result) const
{
    result.Fill(0);
    Times_Add(x,result);
}
//#####################################################################
// Function Transpose_Times_Add
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose_Times_Add(ARRAY_VIEW<const T> x,ARRAY_VIEW<T> result) const
{
    int index=offsets(0);
    for(int i=0;i<m;i++){
        int end=offsets(i+1);T y=x(i);
        for(;index<end;index++) result(A(index).j)+=A(index).a*y;}
}
//#####################################################################
// Function Transpose_Times_Subtract
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose_Times_Subtract(ARRAY_VIEW<const T> x,ARRAY_VIEW<T> result) const
{
    int index=offsets(0);
    for(int i=0;i<m;i++){
        int end=offsets(i+1);T y=x(i);
        for(;index<end;index++) result(A(index).j)-=A(index).a*y;}
}
//#####################################################################
// Function Times
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose_Times(ARRAY_VIEW<const T> x,ARRAY_VIEW<T> result) const
{
    result.Fill(0);
    Transpose_Times_Add(x,result);
}
//#####################################################################
// Function Negate
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Negate()
{
    for(int index=0;index<A.m;index++) A(index).a=-A(index).a;
}
//#####################################################################
// Function operator*=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator*=(const T a)
{
    for(int index=0;index<A.m;index++) A(index).a*=a;
    return *this;
}
//#####################################################################
// Function operator/=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator/=(const T a)
{
    return *this*=1/a;
}
//#####################################################################
// Function operator+=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator+=(const T a)
{
    int index=offsets(0);
    for(int i=0;i<m;i++){
        int end=offsets(i+1),found=0;
        for(;index<end;index++){if(A(index).j==i){found=1;A(index).a+=a;}}
        PHYSBAM_ASSERT(found);}
    return *this;
}
//#####################################################################
// Function operator-=
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T>& SPARSE_MATRIX_FLAT_MXN<T>::
operator-=(const T a)
{
    return *this+=-a;
}
//#####################################################################
// Function Compress
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Compress(SPARSE_MATRIX_FLAT_MXN<T>& compressed)
{
    ARRAY<int> row_lengths(m);
    int index=offsets(0);
    for(int i=0;i<m;i++){
        int end=offsets(i+1);
        while(index<end){if(A(index).j<0) break;
            index++;row_lengths(i)++;}
        index=end;}
    compressed.Set_Row_Lengths(row_lengths);
    index=offsets(0);
    int compressed_index=0;
    for(int i=0;i<m;i++){
        int end=offsets(i+1);
        while(index<end){if(A(index).j<0) break;
            compressed.A(compressed_index++)=A(index);index++;}
        index=end;}
}
//#####################################################################
// Function Transpose
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Transpose(SPARSE_MATRIX_FLAT_MXN<T>& A_transpose) const
{
    ARRAY<int> row_lengths(n);for(int index=0;index<A.m;index++) row_lengths(A(index).j)++;
    A_transpose.Set_Row_Lengths(row_lengths);A_transpose.n=m;
    for(int i=0;i<m;i++) for(int index=offsets(i);index<offsets(i+1);index++) A_transpose(A(index).j,i)=A(index).a;
}
//#####################################################################
// Function Times_Transpose
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
Times_Transpose(const SPARSE_MATRIX_FLAT_MXN<T>& rhs)
{
    assert(rhs.n==n);
    SPARSE_MATRIX_FLAT_MXN<T> result;
    const int columns=rhs.m;const int rows=m;
    ARRAY<int> row_lengths(rows);
    ARRAY<int> right_columns_present,left_rows_present;
    for(int i=0;i<m;i++) if(offsets(i+1)-offsets(i)) left_rows_present.Append(i);
    for(int i=0;i<rhs.m;i++) if(rhs.offsets(i+1)-rhs.offsets(i)) right_columns_present.Append(i);
    
    for(int i=0;i<left_rows_present.m;i++){
        int row=left_rows_present(i);
        int start=offsets(row),end=offsets(row+1);
        for(int j=0;j<right_columns_present.m;j++){
            int col=right_columns_present(j);
            int right_index=rhs.offsets(col);int right_end=rhs.offsets(col+1);int index=start;
            while(index<end && right_index<right_end){
                if(A(index).j==rhs.A(right_index).j){row_lengths(row)++;break;}
                else if(A(index).j>rhs.A(right_index).j) right_index++;
                else index++;}}}
    result.Set_Row_Lengths(row_lengths);result.n=columns;
    for(int i=0;i<left_rows_present.m;i++){
        int row=left_rows_present(i);
        int start=offsets(row),end=offsets(row+1);
        for(int j=0;j<right_columns_present.m;j++){
            int col=right_columns_present(j);
            int right_index=rhs.offsets(col);int right_end=rhs.offsets(col+1);int index=start;
            while(index<end && right_index<right_end){
                if(A(index).j==rhs.A(right_index).j){result.Add_Element(row,col,A(index).a*rhs.A(right_index).a);index++;right_index++;}
                else if(A(index).j>rhs.A(right_index).j) right_index++;
                else index++;}}}
    return result;
}
//#####################################################################
// Function Times_Diagonal_Times
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
Times_Diagonal_Times(ARRAY_VIEW<const T> diagonal,const SPARSE_MATRIX_FLAT_MXN<T>& rhs)
{
    assert(rhs.m==n && diagonal.m==n);
    SPARSE_MATRIX_FLAT_MXN result;
    const int columns=rhs.n,rows=m;
    ARRAY<int> row_lengths(rows);
    ARRAY<int> column_indices;
    for(int row=0;row<rows;row++){
        int start=offsets(row),end=offsets(row+1);
        column_indices.Remove_All();
        for(int i=start;i<end;i++){ // insertion sort to get actual row counts, worst-case cost C^2
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                column_indices.Append_Unique(rhs.A(j).j);
        }
        row_lengths(row)=column_indices.m;
    }
    
    result.Set_Row_Lengths(row_lengths);
    result.n=columns;

    for(int row=0;row<rows;row++){
        int start=offsets(row),end=offsets(row+1);
        for(int i=start;i<end;i++){
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                result.Add_Element(row,rhs.A(j).j,A(i).a*rhs.A(j).a*diagonal(col));
        }
    }
    return result;
}
//#####################################################################
// Function Scale_Rows
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
Scale_Rows(ARRAY_VIEW<const T> d) const
{
    SPARSE_MATRIX_FLAT_MXN<T> result=*this;
    for(int i=0;i<result.m;i++)
        for(int j=result.offsets(i);j<result.offsets(i+1);j++)
            result.A(j).a*=d(i);
    return result;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
operator+(const SPARSE_MATRIX_FLAT_MXN& A_rhs) const
{
    assert(m==A_rhs.m && n==A_rhs.n);
    ARRAY<int> row_lengths(m);
    int left_index=offsets(0),right_index=A_rhs.offsets(0);
    for(int i=0;i<m;i++){
        int left_end=offsets(i+1),right_end=A_rhs.offsets(i+1);
        while(left_index<left_end && right_index<right_end){
            if(A(left_index).j==A_rhs.A(right_index).j){row_lengths(i)++;left_index++;right_index++;}
            else if(A(left_index).j>A_rhs.A(right_index).j){right_index++;row_lengths(i)++;}
            else{left_index++;row_lengths(i)++;}}
        row_lengths(i)+=left_end-left_index+right_end-right_index;
        left_index=left_end;right_index=right_end;}
    SPARSE_MATRIX_FLAT_MXN<T> result;
    result.Set_Row_Lengths(row_lengths);
    result.n=n;

    int index=offsets(0);
    for(int i=0;i<m;i++){int end=offsets(i+1);
        for(;index<end;index++) result(i,A(index).j)=A(index).a;}

    index=A_rhs.offsets(0);
    for(int i=0;i<m;i++){int end=A_rhs.offsets(i+1);
        for(;index<end;index++) result(i,A_rhs.A(index).j)+=A_rhs.A(index).a;}
    return result;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
operator-(const SPARSE_MATRIX_FLAT_MXN& A_rhs) const
{
    assert(m==A_rhs.m && n==A_rhs.n);
    ARRAY<int> row_lengths(m);
    int left_index=offsets(0),right_index=A_rhs.offsets(0);
    for(int i=0;i<m;i++){
        int left_end=offsets(i+1),right_end=A_rhs.offsets(i+1);
        while(left_index<left_end && right_index<right_end){
            if(A(left_index).j==A_rhs.A(right_index).j){row_lengths(i)++;left_index++;right_index++;}
            else if(A(left_index).j>A_rhs.A(right_index).j){right_index++;row_lengths(i)++;}
            else{left_index++;row_lengths(i)++;}}
        row_lengths(i)+=left_end-left_index+right_end-right_index;
        left_index=left_end;right_index=right_end;}
    SPARSE_MATRIX_FLAT_MXN<T> result;
    result.Set_Row_Lengths(row_lengths);result.n=n;

    int index=offsets(0);
    for(int i=0;i<m;i++){int end=offsets(i+1);
        for(;index<end;index++) result(i,A(index).j)=A(index).a;}

    index=A_rhs.offsets(0);
    for(int i=0;i<m;i++){int end=A_rhs.offsets(i+1);
        for(;index<end;index++) result(i,A_rhs.A(index).j)-=A_rhs.A(index).a;}
    return result;
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> SPARSE_MATRIX_FLAT_MXN<T> SPARSE_MATRIX_FLAT_MXN<T>::
operator*(const SPARSE_MATRIX_FLAT_MXN& rhs) const
{
    assert(rhs.m==n);
    SPARSE_MATRIX_FLAT_MXN result;
    const int columns=rhs.n,rows=m;
    ARRAY<int> row_lengths(rows);
    ARRAY<int> column_indices;
    for(int row=0;row<rows;row++){
        int start=offsets(row),end=offsets(row+1);
        column_indices.Remove_All();
        for(int i=start;i<end;i++){ // insertion sort to get actual row counts, worst-case cost C^2
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                column_indices.Append_Unique(rhs.A(j).j);
        }
        row_lengths(row)=column_indices.m;
    }

    result.Set_Row_Lengths(row_lengths);
    result.n=columns;

    for(int row=0;row<rows;row++){
        int start=offsets(row),end=offsets(row+1);
        for(int i=start;i<end;i++){
            int col=A(i).j;
            for(int j=rhs.offsets(col);j<rhs.offsets(col+1);j++)
                result.Add_Element(row,rhs.A(j).j,A(i).a*rhs.A(j).a);
        }
    }
    return result;
}
//#####################################################################
// Function Set_Times_Diagonal
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Set_Times_Diagonal(ARRAY_VIEW<const T> D)
{
    assert(n==D.m);
    for(int i=0;i<A.m;i++)
        A(i).a*=D(A(i).j);
}
//#####################################################################
// Function Set_Diagonal_Times
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Set_Diagonal_Times(ARRAY_VIEW<const T> D)
{
    assert(m==D.m);
    for(int row=0;row<m;row++){
        int start=offsets(row),end=offsets(row+1);
        for(int i=start;i<end;i++)
            A(i).a*=D(row);}
}
//#####################################################################
// Function Write_Row_Lengths
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Write_Row_Lengths()
{
    for(int i=0;i<m;i++) LOG::cout<<offsets(i+1)-offsets(i)<<" ";
    LOG::cout<<std::endl;
}
//#####################################################################
// Function Print_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Print_Row(const int row)
{
    for(int i=0;i<n;i++) if(Element_Present(row,i)) LOG::cout<<"Col: "<<i<<" Val: "<<(*this)(row,i)<<",  ";
    LOG::cout<<std::endl;
}
//#####################################################################
// Function operator<<
//#####################################################################
template<class T> std::ostream&
operator<<(std::ostream& output_stream,const SPARSE_MATRIX_FLAT_MXN<T>& A)
{
    for(int i=0;i<A.m;i++){
        for(int j=0;j<A.n;j++)output_stream<<(A.Element_Present(i,j)?A(i,j):0)<<" ";
        output_stream<<std::endl;}
    return output_stream;
}
//#####################################################################
// Function Reset
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Reset(const int c)
{
    n=c;
    m=0;
    offsets.Remove_All();
    A.Remove_All();
    offsets.Append(0);
    delete C;
    C=0;
    delete L;
    L=0;
    delete Q;
    Q=0;
    diagonal_index.Remove_All();
}
//#####################################################################
// Function Append_Entry_To_Current_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Append_Entry_To_Current_Row(const int c,const T a)
{
    A.Append(SPARSE_MATRIX_ENTRY<T>(c,a));
}
//#####################################################################
// Function Finish_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Finish_Row()
{
    offsets.Append(A.m);
    m++;
}
//#####################################################################
// Function Sort_Entries
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Sort_Entries()
{
    for(int i=0;i<m;i++)
        A.Array_View(offsets(i),offsets(i+1)-offsets(i)).Sort();
}
//#####################################################################
// Function Get_Row
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Get_Row(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q_i,int row)
{    
    q_i.Resize(offsets(row+1)-offsets(row));
    q_i=A.Array_View(offsets(row),offsets(row+1)-offsets(row));
}
//#####################################################################
// Function Construct_Incomplete_QR_Factorization
//#####################################################################
template<class T> void
Check_Sorted(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q)
{
    for(int i=0;i<q.m-1;i++)
        assert(q(i).j<q(i+1).j);
}
//#####################################################################
// Function Keep_Largest_N
//#####################################################################
template<class T> void
Keep_Largest_N(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q,int n)
{
    static ARRAY<SPARSE_MATRIX_ENTRY<T> > elements;
    elements.Resize(0);
    int j;
    for(int i=0;i<q.m;i++){
        T element_magnitude=abs(q(i).a);
        for(j=0;j<elements.m;j++)
            if(element_magnitude>abs(elements(j).a)){
                elements.Insert(q(i),j);
                break;}
        if(j>elements.m && elements.m<n) elements.Append(q(i));
        if(elements.m>n) elements.Resize(n);}
    q=elements;
    q.Sort();
    Check_Sorted(q);
}
//#####################################################################
// Function Subtract_C_Times
//#####################################################################
template<class T> void
Subtract_C_Times(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q,T c,const ARRAY_VIEW<SPARSE_MATRIX_ENTRY<T> >& view)
{
    int q_index=0;
    int current_q_column=0;
    for(int view_index=0;view_index<view.m;view_index++){
        while(q_index<q.m && current_q_column<view(view_index).j){q_index++;current_q_column=q(q_index).j;}
        if(current_q_column==view(view_index).j){
            q(q_index).a-=c*view(view_index).a;
            view_index++;
        }
        else if(current_q_column>view(view_index).j)
            q.Insert(SPARSE_MATRIX_ENTRY<T>(view(view_index).j,-c*view(view_index).a),q_index++);
        else
            q.Insert(SPARSE_MATRIX_ENTRY<T>(view(view_index).j,-c*view(view_index).a),q_index+1);
    }
    Check_Sorted(q);
}
//#####################################################################
// Function Two_Norm
//#####################################################################
template<class T> T
Two_Norm(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q)
{
    double total=0;
    for(int i=0;i<q.m;i++) total+=sqr(q(i).a);
    return (T)sqrt(total);
}
//#####################################################################
// Function Construct_Incomplete_LQ_Factorization
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Construct_Incomplete_LQ_Factorization(const int p_l,const int p_q,const T zero_tolerance,const T zero_replacement)
{
    ARRAY<SPARSE_MATRIX_ENTRY<T> > q_i,l_i;
    delete Q;delete L;Q=new SPARSE_MATRIX_FLAT_MXN<T>();L=new SPARSE_MATRIX_FLAT_MXN<T>();
    Q->Reset(n);L->Reset(m);
    for(int i=0;i<m;i++){
        Get_Row(q_i,i);
        Q->Fast_Sparse_Multiply(q_i,l_i);
        Keep_Largest_N(l_i,p_l);
        for(int j=0;j<l_i.m;j++)
            Subtract_C_Times(q_i,l_i(j).a,Q->A.Array_View(Q->offsets(l_i(j).j),Q->offsets(l_i(j).j+1)-Q->offsets(l_i(j).j)));
        Keep_Largest_N(q_i,p_q);
        T magnitude=Two_Norm(q_i);
        if(magnitude<zero_tolerance){PHYSBAM_FATAL_ERROR("Zero pivot in LQ factorization");magnitude=zero_replacement;}
        T one_over_magnitude=Inverse(magnitude);
        for(int j=0;j<q_i.m;j++)
            Q->Append_Entry_To_Current_Row(q_i(j).j,q_i(j).a*one_over_magnitude);
        Q->Finish_Row();
        for(int j=0;j<l_i.m;j++)
            L->Append_Entry_To_Current_Row(l_i(j).j,l_i(j).a);
        L->Append_Entry_To_Current_Row(i,magnitude);
        L->Finish_Row();
    }
}
//#####################################################################
// Function Fast_Sparse_Multiply
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Fast_Sparse_Multiply(ARRAY<SPARSE_MATRIX_ENTRY<T> >& q,ARRAY<SPARSE_MATRIX_ENTRY<T> >& l)
{
    SPARSE_MATRIX_FLAT_MXN<T> transpose;
    Transpose(transpose);
    l.Resize(0);
    for(int i=0;i<q.m;i++)
        Subtract_C_Times(l,-q(i).a,transpose.A.Array_View(transpose.offsets(q(i).j),transpose.offsets(q(i).j+1)-transpose.offsets(q(i).j)));
}
//#####################################################################
// Function Row_Subset
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Row_Subset(const ARRAY<int>& rows)
{
    delete Q;
    delete L;
    delete C;
    Q=0;
    L=0;
    C=0;
    ARRAY<int> new_offsets;
    ARRAY<SPARSE_MATRIX_ENTRY<T> > new_A;
    new_offsets.Append(0);
    for(int i=0;i<rows.m;i++){
        INTERVAL<int> I(offsets(rows(i)),offsets(rows(i)+1));
        new_offsets.Append(new_offsets.Last()+I.Size());
        new_A.Append_Elements(A.Array_View(I));}
    A.Exchange(new_A);
    offsets.Exchange(new_offsets);
    m=rows.m;
}
//#####################################################################
// Function Column_Subset
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Column_Subset(const ARRAY<int>& cols)
{
    SPARSE_MATRIX_FLAT_MXN<T> tmp;
    Transpose(tmp);
    tmp.Row_Subset(cols);
    tmp.Transpose(*this);
}
//#####################################################################
// Function Initialize_Diagonal_Index
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Initialize_Diagonal_Index()
{
    diagonal_index.Resize(n,false,false);
    for(int i=0;i<n;i++){diagonal_index(i)=Find_Index(i,i);assert(A(diagonal_index(i)).j==i);}
}
//#####################################################################
// Function Solve_Forward_Substitution
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Solve_Forward_Substitution(ARRAY_VIEW<const T> b,ARRAY_VIEW<T> x,const bool diagonal_is_identity,const bool diagonal_is_inverted) const
{
    if(diagonal_is_identity) for(int i=0;i<n;i++){
        T sum=0;for(int index=offsets(i);index<diagonal_index(i);index++)sum+=A(index).a*x(A(index).j);
        x(i)=b(i)-sum;}
    else if(!diagonal_is_inverted) for(int i=0;i<n;i++){
        T sum=0;for(int index=offsets(i);index<diagonal_index(i);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)/A(diagonal_index(i)).a;}
    else for(int i=0;i<n;i++){
        T sum=0;for(int index=offsets(i);index<diagonal_index(i);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)*A(diagonal_index(i)).a;}
}
//#####################################################################
// Function Solve_Backward_Substitution
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Solve_Backward_Substitution(ARRAY_VIEW<const T> b,ARRAY_VIEW<T> x,const bool diagonal_is_identity,const bool diagonal_is_inverted) const
{
    if(diagonal_is_identity) for(int i=n-1;i>=0;i--){
        T sum=0;for(int index=diagonal_index(i)+1;index<offsets(i+1);index++)sum+=A(index).a*x(A(index).j);
        x(i)=b(i)-sum;}
    else if(!diagonal_is_inverted) for(int i=n-1;i>=0;i--){
        T sum=0;for(int index=diagonal_index(i)+1;index<offsets(i+1);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)/A(diagonal_index(i)).a;}
    else for(int i=n-1;i>=0;i--){
        T sum=0;for(int index=diagonal_index(i)+1;index<offsets(i+1);index++)sum+=A(index).a*x(A(index).j);
        x(i)=(b(i)-sum)*A(diagonal_index(i)).a;}
}
//#####################################################################
// Function Construct_Incomplete_Cholesky_Factorization
//#####################################################################
// actually an LU saving square roots, with an inverted diagonal saving divides
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Construct_Incomplete_Cholesky_Factorization(const bool modified_version,const T modified_coefficient,const T zero_tolerance,const T zero_replacement)
{
    delete C;C=new SPARSE_MATRIX_FLAT_MXN<T>(*this);
    C->In_Place_Incomplete_Cholesky_Factorization(modified_version,modified_coefficient,zero_tolerance,zero_replacement);
}
//#####################################################################
// Function In_Place_Incomplete_Cholesky_Factorization
//#####################################################################
// actually an LU saving square roots, with an inverted diagonal saving divides
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
In_Place_Incomplete_Cholesky_Factorization(const bool modified_version,const T modified_coefficient,const T zero_tolerance,const T zero_replacement)
{
    Initialize_Diagonal_Index();
    for(int i=0;i<n;i++){ // for each row
        int row_diagonal_index=diagonal_index(i),row_end=offsets(i+1)-1;T sum=0;
        for(int k_bar=offsets(i);k_bar<row_diagonal_index;k_bar++){ // for all the entries before the diagonal element
            int k=A(k_bar).j;int row2_diagonal_index=diagonal_index(k),row2_end=offsets(k+1)-1;
            A(k_bar).a*=A(row2_diagonal_index).a; // divide by the diagonal element (which has already been inverted)
            int j_bar=k_bar+1; // start with the next element in the row, when subtracting the dot product
            for(int i_bar=row2_diagonal_index+1;i_bar<=row2_end;i_bar++){ // run through the rest of the elements in the row2
                int i=A(i_bar).j;T dot_product_term=A(k_bar).a*A(i_bar).a;
                while(j_bar<row_end && A(j_bar).j<i) j_bar++; // gets j_bar such that j_bar>=i
                if(A(j_bar).j==i) A(j_bar).a-=dot_product_term;else if(modified_version) sum+=dot_product_term;}}
        T denominator=A(row_diagonal_index).a-modified_coefficient*sum;
        if(i==n && denominator<=zero_tolerance) A(row_diagonal_index).a=1/zero_replacement; // ensure last diagonal element is not zero
        else A(row_diagonal_index).a=1/denominator;} // finally, store the diagonal element in inverted form
}
//#####################################################################
// Function Gauss_Seidel_Single_Iteration
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Gauss_Seidel_Single_Iteration(ARRAY_VIEW<T> x,ARRAY_VIEW<const T> b)
{
    assert(x.m==b.m && x.m==n);
    for(int i=0;i<n;i++){
        T rho=0;T diagonal_entry=0;
        for(int index=offsets(i);index<offsets(i+1);index++){
            if(A(index).j==i) diagonal_entry=A(index).a;
            else rho+=A(index).a*x(A(index).j);}
        x(i)=(b(i)-rho)/diagonal_entry;}
}
//#####################################################################
// Function Gauss_Seidel_Solve
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Gauss_Seidel_Solve(ARRAY_VIEW<T> x,ARRAY_VIEW<const T> b,const T tolerance,const int max_iterations)
{
    assert(x.m==b.m && x.m==n);
    ARRAY<T> last_x(x);
    for(int k=0;k<max_iterations;k++){
        Gauss_Seidel_Single_Iteration(x,b);
        T residual=0;for(int j=0;j<n;j++){residual+=sqr(last_x(j)-x(j));last_x(j)=x(j);}if(residual < tolerance) return;}
}
//#####################################################################
// Function Positive_Diagonal_And_Nonnegative_Row_Sum
//#####################################################################
template<class T> bool SPARSE_MATRIX_FLAT_MXN<T>::
Positive_Diagonal_And_Nonnegative_Row_Sum(const T tolerance) const
{
    bool return_value=true;
    for(int i=0;i<n;i++){
        if((*this)(i,i)<=0){
            LOG::cout<<"diagonal entry "<<i<<" contains nonpositive element: "<<(*this)(i,i)<<std::endl;
            return false;}
        T sum=0;for(int index=offsets(i);index<offsets(i+1);index++)sum+=A(index).a;
        if(sum<-tolerance){
            LOG::cout<<"sum of row "<<i<<" is negative: "<<sum<<std::endl;
            return_value=false;}}
    return return_value;
}
//#####################################################################
// Function Conjugate_With_Diagonal_Matrix
//#####################################################################
template<class T> void SPARSE_MATRIX_FLAT_MXN<T>::
Conjugate_With_Diagonal_Matrix(ARRAY_VIEW<T> x)
{
    int index=offsets(0);
    for(int i=0;i<n;i++){
        int end=offsets(i+1);
        for(;index<end;index++) A(index).a*=x(i)*x(A(index).j);}
}
//#####################################################################
template class SPARSE_MATRIX_FLAT_MXN<float>;
template std::ostream& operator<<(std::ostream&,const SPARSE_MATRIX_FLAT_MXN<float>&);
template class SPARSE_MATRIX_FLAT_MXN<double>;
template std::ostream& operator<<(std::ostream&,const SPARSE_MATRIX_FLAT_MXN<double>&);
}
