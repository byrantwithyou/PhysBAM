#include "POISSON_PROJECTION_SYSTEM.h"

template<class TV> POISSON_PROJECTION_SYSTEM<TV>::
POISSON_PROJECTION_SYSTEM()
    :KRYLOV_SYSTEM_BASE<T>(true,true)
{}
template<class TV> POISSON_PROJECTION_SYSTEM<TV>::
~POISSON_PROJECTION_SYSTEM()
{}
template<class TV> void POISSON_PROJECTION_SYSTEM<TV>::
Initialize()
{
    if(!poisson.A.m){
        if(gradient.A.m && !neg_divergence.A.m) gradient.Transpose(neg_divergence);
        if(!gradient.A.m && neg_divergence.A.m) neg_divergence.Transpose(gradient);
        SPARSE_MATRIX_FLAT_MXN<T> tmp(neg_divergence);
        poisson=tmp.Times_Diagonal_Times(beta_inverse,gradient).Create_NXN_Matrix();}
    poisson.Construct_Incomplete_Cholesky_Factorization();
    temp.Resize(poisson.n);
}
template<class TV> void POISSON_PROJECTION_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
    VECTOR_T& vresult=dynamic_cast<VECTOR_T&>(result);
    poisson.Times(vx.v,vresult.v);
}
template<class TV> double POISSON_PROJECTION_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x),vy=dynamic_cast<const VECTOR_T&>(y);
    return vx.v.Dot_Product_Double_Precision(vx.v,vy.v);
}
template<class TV> typename TV::SCALAR POISSON_PROJECTION_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
    return vx.v.Maximum_Magnitude();
}
template<class TV> void POISSON_PROJECTION_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& vx=dynamic_cast<VECTOR_T&>(x);
    for(int i=0;i<projections.m;i++) vx.v-=projections(i)*vx.v.Dot_Product(vx.v,projections(i));
}
template<class TV> void POISSON_PROJECTION_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
template<class TV> void POISSON_PROJECTION_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project(x);
}
template<class TV> void POISSON_PROJECTION_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    const VECTOR_T& vr=dynamic_cast<const VECTOR_T&>(r);
    VECTOR_T& vz=dynamic_cast<VECTOR_T&>(z);
    poisson.C->Solve_Forward_Substitution(vr.v,temp,true);
    poisson.C->Solve_Backward_Substitution(temp,vz.v,false,true);
}
template struct POISSON_PROJECTION_SYSTEM<VECTOR<float,1> >;
template struct POISSON_PROJECTION_SYSTEM<VECTOR<float,2> >;
template struct POISSON_PROJECTION_SYSTEM<VECTOR<float,3> >;
template struct POISSON_PROJECTION_SYSTEM<VECTOR<double,1> >;
template struct POISSON_PROJECTION_SYSTEM<VECTOR<double,2> >;
template struct POISSON_PROJECTION_SYSTEM<VECTOR<double,3> >;
