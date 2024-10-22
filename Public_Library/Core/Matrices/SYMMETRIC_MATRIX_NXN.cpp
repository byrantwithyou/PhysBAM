//#####################################################################
// Copyright 2005-2007, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_NXN
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_NXN.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T>::
SYMMETRIC_MATRIX_NXN(const int n_input)
    :n(n_input),x((n*n+n)/2)
{
}
//#####################################################################
// Function Outer_Product
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T> SYMMETRIC_MATRIX_NXN<T>::
Outer_Product(const ARRAY<T>& u)
{
    SYMMETRIC_MATRIX_NXN<T> result(u.m);
    for(int i=0;i<u.m;i++) for(int j=0;j<i;j++) result(i,j)=u(i)*u(j);
    return result;
}
//#####################################################################
// Function Sqr
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T> SYMMETRIC_MATRIX_NXN<T>::
Sqr() const
{
    SYMMETRIC_MATRIX_NXN<T> result(n);
    for(int j=0;j<n;j++) for(int i=j;i<n;i++){
        for(int k=0;k<j;k++) result(i,j)+=x(((2*n-k-1)*k>>1)+i)*x(((2*n-k-1)*k>>1)+j);
        for(int k=j;k<=i;k++) result(i,j)+=x(((2*n-k-1)*k>>1)+i)*x(((2*n-j-1)*j>>1)+k);
        for(int k=i+1;k<n;k++) result(i,j)+=x(((2*n-i-1)*i>>1)+k)*x(((2*n-j-1)*j>>1)+k);}
    return result;
}
//#####################################################################
// Function Givens_Conjugate
//#####################################################################
template<class T> void SYMMETRIC_MATRIX_NXN<T>::
Givens_Conjugate(const int i,const int j,const T c,const T s)
{
    if(i>j){Givens_Conjugate(j,i,c,-s);return;}
    assert(0<=i && i<j && j<n);
    for(int k=0;k<i;k++){T u=(*this)(i,k),v=(*this)(j,k);(*this)(i,k)=c*u-s*v;(*this)(j,k)=s*u+c*v;}
    for(int k=i+1;k<j;k++){T u=(*this)(k,i),v=(*this)(j,k);(*this)(k,i)=c*u-s*v;(*this)(j,k)=s*u+c*v;}
    for(int k=j+1;k<n;k++){T u=(*this)(k,i),v=(*this)(k,j);(*this)(k,i)=c*u-s*v;(*this)(k,j)=s*u+c*v;}
    T u=(*this)(i,i),v=(*this)(j,i),w=(*this)(j,j);
    (*this)(i,i)=c*c*u-(T)2*c*s*v+s*s*w;
    (*this)(j,i)=(c*c-s*s)*v+c*s*(u-w);
    (*this)(j,j)=c*c*w+(T)2*c*s*v+s*s*u;
}
//#####################################################################
// Function Jacobi_Solve_Eigenproblem
//#####################################################################
template<class T> void SYMMETRIC_MATRIX_NXN<T>::
Jacobi_Solve_Eigenproblem(ARRAY<VECTOR<int,2> >& givens_pairs,ARRAY<VECTOR<T,2> >& givens_coefficients,const T tolerance,const int max_iterations)
{
    assert(n>=2);
    givens_pairs.Resize(0);
    givens_coefficients.Resize(0);
    for(int iteration=0;iteration<max_iterations;iteration++){
        T max_off_diagonal_element=0;
        int i_max=0,j_max=0;
        for(int j=0;j<n;j++) for(int i=j+1;i<n;i++) if(abs((*this)(i,j))>max_off_diagonal_element){max_off_diagonal_element=abs((*this)(i,j));i_max=i;j_max=j;}
        if(max_off_diagonal_element<tolerance) return;
        T q=(T).5*((*this)(i_max,i_max)-(*this)(j_max,j_max))/(*this)(i_max,j_max),t,c,s;
        t=(q>0)?(T)1/(q+sqrt(q*q+(T)1)):(T)1/(q-sqrt(q*q+(T)1));
        c=(T)1/sqrt(t*t+(T)1);
        s=c*t;
        Givens_Conjugate(j_max,i_max,c,s);
        (*this)(j_max,i_max)=0;
        givens_pairs.Append(VECTOR<int,2>(j_max,i_max));
        givens_coefficients.Append(VECTOR<T,2>(c,s));}
}
//#####################################################################
// Function Maximum_Eigenvalue_Eigenvector_Pair
//#####################################################################
template<class T>
void SYMMETRIC_MATRIX_NXN<T>::
Maximum_Eigenvalue_Eigenvector_Pair(T& max_eigenvalue,ARRAY<T>& max_eigenvector,RANDOM_NUMBERS<T>* random_numbers,const T tolerance,
    const T randomization_decay_factor,const int max_iterations)
{
    ARRAY<T> last_eigenvector(n);
    T randomization_factor=(T)1,tolerance_squared=sqr(tolerance),tolerance_scaled=tolerance/sqrt((T)n);
    max_eigenvector.Resize(n); // Must provide initial guess if no randomization is used
    for(int iteration=0;iteration<max_iterations;iteration++){
        last_eigenvector=max_eigenvector;
        if(random_numbers) for(int i=0;i<n;i++) last_eigenvector(i)+=random_numbers->Get_Uniform_Number(-randomization_factor,randomization_factor);
        last_eigenvector/=ARRAY<T>::Dot_Product_Double_Precision(last_eigenvector,last_eigenvector);
        max_eigenvector=(*this)*last_eigenvector;
        max_eigenvalue=ARRAY<T>::Dot_Product_Double_Precision(max_eigenvector,max_eigenvector);
        max_eigenvalue=sqrt(max_eigenvalue);
        if(ARRAY<T>::Dot_Product(last_eigenvector,max_eigenvector)<0) max_eigenvalue=-max_eigenvalue;
        max_eigenvector/=max_eigenvalue;
        last_eigenvector-=max_eigenvector;
        if((!random_numbers || randomization_factor<tolerance_scaled) && last_eigenvector.Magnitude_Squared()<tolerance_squared) return;
        if(random_numbers) randomization_factor*=randomization_decay_factor;}
}
//#####################################################################
// Function In_Place_Cholesky_Factorization
//#####################################################################
template<class T> void SYMMETRIC_MATRIX_NXN<T>::
In_Place_Cholesky_Factorization(MATRIX_MXN<T>& L)
{
    L=MATRIX_MXN<T>(n);
    for(int j=0;j<n;j++){ // for each column
        for(int k=0;k<j-1;k++) for(int i=j;i<n;i++) Element_Lower(i,j)-=L(j,k)*L(i,k); // subtract off the known stuff in previous columns
        L(j,j)=sqrt(Element_Lower(j,j));
        T diagonal_inverse=1/L(j,j);
        for(int i=j+1;i<n;i++) L(i,j)=Element_Lower(i,j)*diagonal_inverse;} // update L
}
//#####################################################################
// Function operator+=
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T>& SYMMETRIC_MATRIX_NXN<T>::
operator+=(const SYMMETRIC_MATRIX_NXN<T>& A)
{
    x+=A.x;
    return *this;
}
//#####################################################################
// Function operator-=
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T>& SYMMETRIC_MATRIX_NXN<T>::
operator-=(const SYMMETRIC_MATRIX_NXN<T>& A)
{
    x-=A.x;
    return *this;
}
//#####################################################################
// Function operator*=
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T>& SYMMETRIC_MATRIX_NXN<T>::
operator*=(const T a)
{
    x*=a;
    return *this;
}
//#####################################################################
// Function operator+
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T> SYMMETRIC_MATRIX_NXN<T>::
operator+(const SYMMETRIC_MATRIX_NXN<T>& A) const
{
    assert(A.n==n);
    SYMMETRIC_MATRIX_NXN<T> result(n);
    result.x=x+A.x;
    return result;
}
//#####################################################################
// Function operator-
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T> SYMMETRIC_MATRIX_NXN<T>::
operator-(const SYMMETRIC_MATRIX_NXN<T>& A) const
{
    assert(A.n==n);
    SYMMETRIC_MATRIX_NXN<T> result(n);
    result.x=x-A.x;
    return result;
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T> SYMMETRIC_MATRIX_NXN<T>::
operator*(const T a) const
{
    SYMMETRIC_MATRIX_NXN<T> result(n);
    result.x=x*a;
    return result;
}
//#####################################################################
// Function operator*
//#####################################################################
template<class T> ARRAY<T> SYMMETRIC_MATRIX_NXN<T>::
operator*(const ARRAY<T>& y) const
{
    assert(y.m==n);
    ARRAY<T> result(n);
    for(int i=0;i<n;i++){
        for(int j=0;j<=i;j++)
            result(i)+=x(((2*n-j-1)*j>>1)+i)*y(j);
        for(int j=i+1;j<n;j++)
            result(i)+=x(((2*n-i-1)*i>>1)+j)*y(j);}
    return result;
}
//#####################################################################
// Function Set_Identity_Matrix
//#####################################################################
template<class T> void SYMMETRIC_MATRIX_NXN<T>::
Set_Identity_Matrix()
{
    Set_Zero_Matrix();
    for(int i=0,j=n;j>0;i+=j--) x(i)=1;
}
//#####################################################################
// Function Set_Zero_Matrix
//#####################################################################
template<class T> void SYMMETRIC_MATRIX_NXN<T>::
Set_Zero_Matrix()
{
    x.Fill(0);
}
//#####################################################################
// Function Identity_Matrix
//#####################################################################
template<class T> SYMMETRIC_MATRIX_NXN<T> SYMMETRIC_MATRIX_NXN<T>::
Identity_Matrix(const int n)
{
    SYMMETRIC_MATRIX_NXN<T> A(n);
    for(int i=0,j=n;j>0;i+=j--) A.x(i)=1;
    return A;
}
//#####################################################################
// Function Trace
//#####################################################################
template<class T> T SYMMETRIC_MATRIX_NXN<T>::
Trace() const
{
    T trace=0;
    for(int i=0;i<n;i++) trace+=(*this)(i,i);
    return trace;
}
//#####################################################################
// Function operator<<
//#####################################################################
template<class T> std::ostream&
operator<<(std::ostream& output_stream,const SYMMETRIC_MATRIX_NXN<T>& A)
{
    for(int i=0;i<A.n;i++){for(int j=0;j<A.n;j++) output_stream<<A(i,j)<<" ";output_stream<<std::endl;}
    return output_stream;
}
//#####################################################################
namespace PhysBAM{
template class SYMMETRIC_MATRIX_NXN<float>;
template class SYMMETRIC_MATRIX_NXN<double>;
}
