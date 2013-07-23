//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_BASE
//#####################################################################
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_BASE.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Vectors/VECTOR.h>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Function In_Place_Cholesky_Inverse
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Cholesky_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse)
{
    assert(Rows()==Columns());
    inverse.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    LEFT_VECTOR b((INITIAL_SIZE)Rows()); // holds piece of the identity matrix
    In_Place_Cholesky_Factorization();
    for(int j=0;j<Columns();j++){
        b(j)=1;
        inverse.Set_Column(j,Transpose_Upper_Triangular_Solve(Lower_Triangular_Solve(b)));
        b(j)=0;}
}
//#####################################################################
// Function In_Place_PLU_Inverse
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_PLU_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse)
{
    assert(Rows()==Columns());
    inverse.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    COLUMN_PERMUTATION p;
    T_MATRIX L;
    In_Place_PLU_Factorization(L,p);
    RIGHT_VECTOR b((INITIAL_SIZE)Columns()); // used for piece of the identity matrix
    for(int j=0;j<Columns();j++){
        b(j)=0;
        inverse.Set_Column(j,Upper_Triangular_Solve(L.Lower_Triangular_Solve(b.Permute(p))));
        b(j)=0;}
}
//#####################################################################
// Function In_Place_LU_Inverse
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_LU_Inverse(MATRIX_BASE<T,T_MATRIX2>& inverse) // don't assume ARRAY
{
    assert(Rows()==Columns());
    inverse.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    T_MATRIX L;
    RIGHT_VECTOR b(Columns()); // used forpiece of the identity matrix
    In_Place_LU_Factorization(L);
    for(int j=0;j<Columns();j++){
        b(j)=1;
        inverse.Set_Column(j,Upper_Triangular_Solve(L.Lower_Triangular_Solve(b)));
        b(j)=0;}
}
//#####################################################################
// Function Gram_Schmidt_QR_Factorization
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Gram_Schmidt_QR_Factorization(MATRIX_BASE<T,T_MATRIX2>& R)
{
    R.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=0;j<Columns();j++){ // for each column
        for(int i=0;i<Rows();i++) R(j,j)+=sqr((*this)(i,j));R(j,j)=sqrt(R(j,j)); // compute the L2 norm
        T one_over_Rjj=1/R(j,j);
        for(int i=0;i<Rows();i++) (*this)(i,j)*=one_over_Rjj; // orthogonalize the column
        for(int k=j+1;k<Columns();k++){ // subtract this columns contributution from the rest of the columns
            for(int i=0;i<Rows();i++) R(j,k)+=(*this)(i,j)*(*this)(i,k);
            for(int i=0;i<Rows();i++) (*this)(i,k)-=R(j,k)*(*this)(i,j);}}
}
//#####################################################################
// Function Householder_QR_Factorization
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2,class T_MATRIX3> void MATRIX_BASE<T,T_MATRIX>::
Householder_QR_Factorization(MATRIX_BASE<T,T_MATRIX2>& V,MATRIX_BASE<T,T_MATRIX3>& R)
{
    V.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());R.Derived()=T_MATRIX3(Columns(),Columns());T_MATRIX temp(*this);LEFT_VECTOR a(Rows()),v,new_a;
    for(int j=0;j<Columns();j++){ // for each column
        for(int i=0;i<Rows();i++) a(i)=temp(i,j);
        v=a.Householder_Vector(j);for(int i=0;i<Rows();i++) V(i,j)=v(i); // store the v's in V
        for(int k=j;k<Columns();k++){ // Householder transform each column
            for(int i=0;i<Rows();i++) a(i)=temp(i,k);
            new_a=a.Householder_Transform(v);for(int i=0;i<Rows();i++) temp(i,k)=new_a(i);}}
    for(int i=0;i<Columns();i++) for(int j=0;j<Columns();j++) R(i,j)=temp(i,j); // store R
}
//#####################################################################
// Function Robust_Householder_QR_Solve
//#####################################################################
template<class T,class T_MATRIX> template<class T_VECTOR1,class T_VECTOR2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Robust_Householder_QR_Solve(ARRAY_BASE<T,T_VECTOR1>& b,ARRAY_BASE<int,T_VECTOR2>& p)
{
    assert(Rows()==b.Size() && Columns()==p.Size());
    ARRAY<T> a((INITIAL_SIZE)Rows());for(int i=0;i<Columns();i++) p(i)=i; // TODO: This should not assume ARRAY.
    ARRAY<T> column_norm(Columns());for(int j=0;j<Columns();j++) for(int i=0;i<Rows();i++) column_norm(j)+=sqr((*this)(i,j));
    for(int j=0;j<Columns();j++){
        int max_column=0;T max_column_norm=0;for(int k=j;k<Columns();k++) if(column_norm(k)>max_column_norm){max_column_norm=column_norm(k);max_column=k;}
        if(max_column_norm<FLT_MIN) return;
        if(max_column!=j){exchange(column_norm(j),column_norm(max_column));exchange(p(j),p(max_column));for(int i=0;i<Rows();i++) exchange((*this)(i,j),(*this)(i,max_column));}
        if(j==Rows()) return;
        Get_Column(j,a);ARRAY<T> v=a.Householder_Vector(j);T two_over_v_dot_v=(T)2/v.Magnitude_Squared();
        if((*this)(j,j)>=0)(*this)(j,j)=-sqrt(max_column_norm);else (*this)(j,j)=sqrt(max_column_norm);for(int i=j+1;i<Rows();i++)(*this)(i,j)=(T)0;
        for(int k=j+1;k<Columns();k++){
            T v_dot_a=0;for(int i=j;i<Rows();i++) v_dot_a+=v(i)*(*this)(i,k);T coefficient=v_dot_a*two_over_v_dot_v;for(int i=j;i<Rows();i++) (*this)(i,k)-=coefficient*v(i);}
        T v_dot_b=0;for(int i=j;i<Rows();i++) v_dot_b+=v(i)*b(i);T coefficient=v_dot_b*two_over_v_dot_v;for(int i=j;i<Rows();i++) b(i)-=coefficient*v(i);
        for(int k=j+1;k<Columns();k++) column_norm(k)-=sqr((*this)(j,k));}
}
//#####################################################################
// Function Robust_Householder_QR_Solve
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_PLU_Factorization(MATRIX_BASE<T,T_MATRIX2>& L,COLUMN_PERMUTATION& p)
{
    assert((INITIAL_SIZE)Rows()==(INITIAL_SIZE)Columns());
    L.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    p=COLUMN_PERMUTATION(INITIAL_SIZE(Columns()));for(int i=0;i<Columns();i++) p(i)=i; // initialize p
    for(int j=0;j<Columns();j++){ // for each column
        // find the largest element and switch rows
        int row=j;T value=abs((*this)(j,j));
        for(int i=j+1;i<Columns();i++) if(abs((*this)(i,j))>value){row=i;value=abs((*this)(i,j));}
        if(row!=j){ // need to switch rows
            exchange(p(j),p(row)); // update permutation matrix
            for(int k=0;k<j;k++) exchange(L(j,k),L(row,k)); // update L
            for(int k=j;k<Columns();k++) exchange((*this)(j,k),(*this)(row,k));} // update U
        // standard LU factorization steps
        T diagonal_inverse=1/(*this)(j,j);for(int i=j;i<Columns();i++) L(i,j)=(*this)(i,j)*diagonal_inverse; // fill in the column for L
        for(int i=j+1;i<Columns();i++) for(int k=j;k<Columns();k++) (*this)(i,k)-=L(i,j)*(*this)(j,k);} // sweep across each row below row j  TODO: can order be changed?
}
//#####################################################################
// Function In_Place_Cholesky_Factorization
//#####################################################################
template<class T,class T_MATRIX> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Cholesky_Factorization()
{
    assert(Rows()==Columns());
    for(int j=0;j<Columns();j++){ // for each column
        for(int k=0;k<j;k++) for(int i=j;i<Rows();i++) (*this)(i,j)-=(*this)(j,k)*(*this)(i,k); // subtract off the known stuff in previous columns
        (*this)(j,j)=sqrt((*this)(j,j));T diagonal_inverse=1/(*this)(j,j);for(int i=j+1;i<Columns();i++) (*this)(i,j)*=diagonal_inverse;} // update L
    for(int i=0;i<Rows();i++) for(int j=i+1;j<Columns();j++) (*this)(i,j)=0; // zero out upper triangular part  TODO: Loop the other way around
}
//#####################################################################
// Function In_Place_Cholesky_Factorization
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_LU_Factorization(MATRIX_BASE<T,T_MATRIX2>& L)
{
    assert(Rows()==Columns());
    L.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=0;j<Columns();j++){ // for each column
        T diagonal_inverse=1/(*this)(j,j);for(int i=j;i<Columns();i++) L(i,j)=(*this)(i,j)*diagonal_inverse; // fill in the column for L
        for(int i=j+1;i<Columns();i++) for(int k=j;k<Columns();k++) (*this)(i,k)-=L(i,j)*(*this)(j,k);} // sweep across each row below row j  TODO: can the order of these loops be swapped?
}
//####################################################################################
// Function Number_Of_Nonzero_Rows
//####################################################################################
template<class T,class T_MATRIX> int MATRIX_BASE<T,T_MATRIX>::
Number_Of_Nonzero_Rows(const T threshold) const
{
    T threshold_squared=sqr(threshold);int nonzero_rows=0;
    for(int i=0;i<Rows();i++){
        T row_norm_squared=0;for(int j=0;j<Columns();j++) row_norm_squared+=sqr((*this)(i,j));
        if(row_norm_squared>threshold_squared) nonzero_rows++;}
    return nonzero_rows;
}
namespace{
template<class T> static void 
Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(const MATRIX_MXN<T>& A,ARRAY<PAIR<int,T> >& max_off_diagonal_element_of_row,const int i)
{max_off_diagonal_element_of_row(i)=PAIR<int,T>();
for(int j=0;j<A.n;j++) if(i!=j) 
    if(abs(A(i,j))>max_off_diagonal_element_of_row(i).y){max_off_diagonal_element_of_row(i).y=abs(A(i,j));max_off_diagonal_element_of_row(i).x=j;}}
template<class T> static void 
Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(const MATRIX_MXN<T>& A,ARRAY<PAIR<int,T> >& max_off_diagonal_element_of_row,const int j)
{for(int i=0;i<A.m;i++) if(i!=j){
    if(abs(A(i,j))>max_off_diagonal_element_of_row(i).y){max_off_diagonal_element_of_row(i).y=abs(A(i,j));max_off_diagonal_element_of_row(i).x=j;}
    else if(max_off_diagonal_element_of_row(i).x==j) Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(A,max_off_diagonal_element_of_row,i);}}
}
//#####################################################################
// Function Jacobi_Singular_Value_Decomposition
//#####################################################################
template<class T,class T_MATRIX> void MATRIX_BASE<T,T_MATRIX>::
Jacobi_Singular_Value_Decomposition(ARRAY<VECTOR<int,2> >& left_givens_pairs,ARRAY<VECTOR<T,2> >& left_givens_coefficients,
    ARRAY<VECTOR<int,2> >& right_givens_pairs,ARRAY<VECTOR<T,2> >& right_givens_coefficients,const T tolerance,const int max_iterations)
{
    assert(Rows()>=2 && Columns()>=2);
    ARRAY<PAIR<int,T> > max_off_diagonal_element_of_row(Rows());
    for(int i=0;i<Rows();i++) Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,i);
    left_givens_pairs.Remove_All();left_givens_coefficients.Remove_All();
    right_givens_pairs.Remove_All();right_givens_coefficients.Remove_All();
    for(int iteration=0;iteration<max_iterations;iteration++){
        T max_off_diagonal_element=0;int i_max=0,j_max=0;
        for(int i=0;i<Rows();i++)
            if(max_off_diagonal_element_of_row(i).y>max_off_diagonal_element){max_off_diagonal_element=max_off_diagonal_element_of_row(i).y;i_max=i;j_max=max_off_diagonal_element_of_row(i).x;}
        if(max_off_diagonal_element<tolerance) return;
        if(i_max>Columns()){
            int i=j_max,j=i_max;T c,s;VECTOR<T,2>((*this)(i,i),(*this)(j,i)).Normalized().Get(c,s);
            Left_Givens_Rotation(i,j,c,-s);left_givens_pairs.Append(VECTOR<int,2>(i,j));left_givens_coefficients.Append(VECTOR<T,2>(c,s));
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,j);}
        else if(j_max>Rows()){
            int i=i_max,j=j_max;T c,s;VECTOR<T,2>((*this)(i,i),(*this)(i,j)).Normalized().Get(c,s);
            Right_Givens_Rotation(i,j,c,-s);right_givens_pairs.Append(VECTOR<int,2>(i,j));right_givens_coefficients.Append(VECTOR<T,2>(c,s));
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,j);}
        else{
            int i=(i_max<j_max)?i_max:j_max,j=(i_max<j_max)?j_max:i_max;
            MATRIX<T,2> B((*this)(i,i),(*this)(j,i),(*this)(i,j),(*this)(j,j)),U,V;DIAGONAL_MATRIX<T,2> sigma;B.Fast_Singular_Value_Decomposition(U,sigma,V);
            T c_left,s_left,c_right,s_right;U.Column(1).Get(c_left,s_left);V.Column(1).Get(c_right,s_right);
            Left_Givens_Rotation(i,j,c_left,-s_left);left_givens_pairs.Append(VECTOR<int,2>(i,j));left_givens_coefficients.Append(VECTOR<T,2>(c_left,s_left));
            Right_Givens_Rotation(i,j,c_right,-s_right);right_givens_pairs.Append(VECTOR<int,2>(i,j));right_givens_coefficients.Append(VECTOR<T,2>(c_right,s_right));
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Row_Change(*this,max_off_diagonal_element_of_row,j);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,i);
            Update_Max_Off_Diagonal_Element_Of_Row_After_Column_Change(*this,max_off_diagonal_element_of_row,j);}}
}
//#####################################################################
namespace PhysBAM{
template void MATRIX_BASE<float,MATRIX_MXN<float> >::In_Place_Robust_Householder_QR_Solve<ARRAY<float>,ARRAY<int> >(ARRAY_BASE<float,ARRAY<float> >&,ARRAY_BASE<int,ARRAY<int> >&);
template void MATRIX_BASE<float,MATRIX<float,6,6> >::In_Place_Gram_Schmidt_QR_Factorization<MATRIX<float,6,6> >(MATRIX_BASE<float,MATRIX<float,6,6> >&);
template void MATRIX_BASE<float,MATRIX_MXN<float> >::In_Place_LU_Factorization<MATRIX_MXN<float> >(MATRIX_BASE<float,MATRIX_MXN<float> >&);
template int MATRIX_BASE<float,MATRIX_MXN<float> >::Number_Of_Nonzero_Rows(const float threshold) const;
template void MATRIX_BASE<float,MATRIX_MXN<float> >::Householder_QR_Factorization<MATRIX_MXN<float>,MATRIX_MXN<float> >(MATRIX_BASE<float,MATRIX_MXN<float> >&,MATRIX_BASE<float,MATRIX_MXN<float> >&);
template void MATRIX_BASE<double,MATRIX_MXN<double> >::In_Place_Robust_Householder_QR_Solve<ARRAY<double>,ARRAY<int> >(ARRAY_BASE<double,ARRAY<double> >&,ARRAY_BASE<int,ARRAY<int> >&);
template void MATRIX_BASE<double,MATRIX<double,6,6> >::In_Place_Gram_Schmidt_QR_Factorization<MATRIX<double,6,6> >(MATRIX_BASE<double,MATRIX<double,6,6> >&);
template void MATRIX_BASE<double,MATRIX_MXN<double> >::In_Place_LU_Factorization<MATRIX_MXN<double> >(MATRIX_BASE<double,MATRIX_MXN<double> >&);
template int MATRIX_BASE<double,MATRIX_MXN<double> >::Number_Of_Nonzero_Rows(const double threshold) const;
template void MATRIX_BASE<double,MATRIX_MXN<double> >::Householder_QR_Factorization<MATRIX_MXN<double>,MATRIX_MXN<double> >(MATRIX_BASE<double,MATRIX_MXN<double> >&,MATRIX_BASE<double,MATRIX_MXN<double> >&);

template void MATRIX_BASE<float,MATRIX<float,0,0> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<float,MATRIX<float,1,1> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<float,MATRIX<float,2,2> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<float,MATRIX<float,3,3> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<float,MATRIX<float,4,4> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<float,MATRIX<float,5,5> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<float,MATRIX<float,6,6> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<float,MATRIX_MXN<float> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX<double,0,0> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX<double,1,1> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX<double,2,2> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX<double,3,3> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX<double,4,4> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX<double,5,5> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX<double,6,6> >::In_Place_Cholesky_Factorization();
template void MATRIX_BASE<double,MATRIX_MXN<double> >::In_Place_Cholesky_Factorization();

template void MATRIX_BASE<float,MATRIX<float,0,0> >::In_Place_PLU_Factorization<MATRIX<float,0,0> >(MATRIX_BASE<float,MATRIX<float,0,0> >&,VECTOR<int,0>&);
template void MATRIX_BASE<float,MATRIX<float,1,1> >::In_Place_PLU_Factorization<MATRIX<float,1,1> >(MATRIX_BASE<float,MATRIX<float,1,1> >&,VECTOR<int,1>&);
template void MATRIX_BASE<float,MATRIX<float,2,2> >::In_Place_PLU_Factorization<MATRIX<float,2,2> >(MATRIX_BASE<float,MATRIX<float,2,2> >&,VECTOR<int,2>&);
template void MATRIX_BASE<float,MATRIX<float,3,3> >::In_Place_PLU_Factorization<MATRIX<float,3,3> >(MATRIX_BASE<float,MATRIX<float,3,3> >&,VECTOR<int,3>&);
template void MATRIX_BASE<float,MATRIX<float,4,4> >::In_Place_PLU_Factorization<MATRIX<float,4,4> >(MATRIX_BASE<float,MATRIX<float,4,4> >&,VECTOR<int,4>&);
template void MATRIX_BASE<float,MATRIX<float,5,5> >::In_Place_PLU_Factorization<MATRIX<float,5,5> >(MATRIX_BASE<float,MATRIX<float,5,5> >&,VECTOR<int,5>&);
template void MATRIX_BASE<float,MATRIX<float,6,6> >::In_Place_PLU_Factorization<MATRIX<float,6,6> >(MATRIX_BASE<float,MATRIX<float,6,6> >&,VECTOR<int,6>&);
template void MATRIX_BASE<float,MATRIX_MXN<float> >::In_Place_PLU_Factorization<MATRIX_MXN<float> >(MATRIX_BASE<float,MATRIX_MXN<float> >&,ARRAY<int>&);
template void MATRIX_BASE<double,MATRIX<double,0,0> >::In_Place_PLU_Factorization<MATRIX<double,0,0> >(MATRIX_BASE<double,MATRIX<double,0,0> >&,VECTOR<int,0>&);
template void MATRIX_BASE<double,MATRIX<double,1,1> >::In_Place_PLU_Factorization<MATRIX<double,1,1> >(MATRIX_BASE<double,MATRIX<double,1,1> >&,VECTOR<int,1>&);
template void MATRIX_BASE<double,MATRIX<double,2,2> >::In_Place_PLU_Factorization<MATRIX<double,2,2> >(MATRIX_BASE<double,MATRIX<double,2,2> >&,VECTOR<int,2>&);
template void MATRIX_BASE<double,MATRIX<double,3,3> >::In_Place_PLU_Factorization<MATRIX<double,3,3> >(MATRIX_BASE<double,MATRIX<double,3,3> >&,VECTOR<int,3>&);
template void MATRIX_BASE<double,MATRIX<double,4,4> >::In_Place_PLU_Factorization<MATRIX<double,4,4> >(MATRIX_BASE<double,MATRIX<double,4,4> >&,VECTOR<int,4>&);
template void MATRIX_BASE<double,MATRIX<double,5,5> >::In_Place_PLU_Factorization<MATRIX<double,5,5> >(MATRIX_BASE<double,MATRIX<double,5,5> >&,VECTOR<int,5>&);
template void MATRIX_BASE<double,MATRIX<double,6,6> >::In_Place_PLU_Factorization<MATRIX<double,6,6> >(MATRIX_BASE<double,MATRIX<double,6,6> >&,VECTOR<int,6>&);
template void MATRIX_BASE<double,MATRIX_MXN<double> >::In_Place_PLU_Factorization<MATRIX_MXN<double> >(MATRIX_BASE<double,MATRIX_MXN<double> >&,ARRAY<int>&);
template void MATRIX_BASE<double,MATRIX<double,1,1> >::In_Place_Cholesky_Inverse<MATRIX<double,1,1> >(MATRIX_BASE<double,MATRIX<double,1,1> >&);
template void MATRIX_BASE<double,MATRIX<double,3,3> >::In_Place_Cholesky_Inverse<MATRIX<double,3,3> >(MATRIX_BASE<double,MATRIX<double,3,3> >&);
template void MATRIX_BASE<double,MATRIX<double,6,6> >::In_Place_Cholesky_Inverse<MATRIX<double,6,6> >(MATRIX_BASE<double,MATRIX<double,6,6> >&);
template void MATRIX_BASE<double,MATRIX_MXN<double> >::In_Place_Cholesky_Inverse<MATRIX_MXN<double> >(MATRIX_BASE<double,MATRIX_MXN<double> >&);
template void MATRIX_BASE<double,MATRIX_MXN<double> >::In_Place_LU_Inverse<MATRIX_MXN<double> >(MATRIX_BASE<double,MATRIX_MXN<double> >&);
template void MATRIX_BASE<float,MATRIX<float,1,1> >::In_Place_Cholesky_Inverse<MATRIX<float,1,1> >(MATRIX_BASE<float,MATRIX<float,1,1> >&);
template void MATRIX_BASE<float,MATRIX<float,3,3> >::In_Place_Cholesky_Inverse<MATRIX<float,3,3> >(MATRIX_BASE<float,MATRIX<float,3,3> >&);
template void MATRIX_BASE<float,MATRIX<float,6,6> >::In_Place_Cholesky_Inverse<MATRIX<float,6,6> >(MATRIX_BASE<float,MATRIX<float,6,6> >&);
template void MATRIX_BASE<float,MATRIX_MXN<float> >::In_Place_Cholesky_Inverse<MATRIX_MXN<float> >(MATRIX_BASE<float,MATRIX_MXN<float> >&);
template void MATRIX_BASE<float,MATRIX_MXN<float> >::In_Place_LU_Inverse<MATRIX_MXN<float> >(MATRIX_BASE<float,MATRIX_MXN<float> >&);
}
