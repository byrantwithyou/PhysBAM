//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Math_Tools/RANGE_ITERATOR.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS_SPLINE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int degree> PARTICLE_GRID_WEIGHTS_SPLINE<TV,degree>::
PARTICLE_GRID_WEIGHTS_SPLINE(const GRID<TV>& grid,int threads)
    :BASE(threads),grid(grid)
{
    this->use_gradient_transfer=(degree==1);
    this->constant_scalar_inertia_tensor=(degree!=1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int degree> PARTICLE_GRID_WEIGHTS_SPLINE<TV,degree>::
~PARTICLE_GRID_WEIGHTS_SPLINE()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV,int degree> void PARTICLE_GRID_WEIGHTS_SPLINE<TV,degree>::
Compute(int p,typename PARTICLE_GRID_ITERATOR<TV>::SCRATCH& scratch,bool want_gradient) const
{
    const PRECOMPUTE_DATA& pd=precompute_data(p);
    int pow_n=pow<TV::m,int>(n);
    scratch.index.Resize(pow_n);
    scratch.weight.Resize(pow_n);
    if(want_gradient) scratch.gradient.Resize(pow_n);

    int k=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+n));it.Valid();it.Next()){
        scratch.index(k)=pd.base+it.index;
        TV all_w;
        for(int i=0;i<TV::m;i++)
            all_w(i)=pd.w(it.index(i))(i);
        scratch.weight(k)=all_w.Product();
        if(want_gradient){
            TV dw(TV()+1);
            for(int j=0;j<TV::m;j++){
                T dj=pd.dw(it.index(j))(j);
                for(int i=0;i<TV::m;i++)
                    dw(i)*=(i==j)?dj:all_w(j);}
            scratch.gradient(k)=dw;}
        k++;}
}
//#####################################################################
// Function Compute_Weights
//#####################################################################
template<class TV> static void
Compute_Weights(VECTOR<TV,2>& w,VECTOR<TV,2>& dw,TV x,TV one_over_dX)
{
    typedef typename TV::SCALAR T;
    for(int k=0;k<TV::m;++k)
        if(x(k)<1e-24) x(k)=1e-10;
    w(0)=(T)1-x;
    dw(0)=-one_over_dX;
    x-=1;
    w(1)=(T)1+x;
    dw(1)=one_over_dX;
}
//#####################################################################
// Function Compute_Weights
//#####################################################################
template<class TV> static void
Compute_Weights(VECTOR<TV,3>& w,VECTOR<TV,3>& dw,TV x,TV one_over_dX)
{
    typedef typename TV::SCALAR T;
    w(0)=(T).5*sqr(x-(T)1.5);
    dw(0)=(x-(T)1.5)*one_over_dX;
    x-=(T)1;
    w(1)=-x*x+(T).75;
    dw(1)=(T)(-2)*x*one_over_dX;
    x-=(T)1;
    w(2)=(T).5*sqr(x+(T)1.5);
    dw(2)=(x+(T)1.5)*one_over_dX;
}
//#####################################################################
// Function Compute_Weights
//#####################################################################
template<class TV> static void
Compute_Weights(VECTOR<TV,4>& w,VECTOR<TV,4>& dw,TV x,TV one_over_dX)
{
    typedef typename TV::SCALAR T;
    x=clamp(x,TV()+1,TV()+2);
    TV x2,x3,z=(T)2-x,z2=z*z;
    w(0)=((T)1/6)*z2*z;
    dw(0)=-(T).5*one_over_dX*z2;
    x-=(T)1;
    x2=x*x;x3=x2*x;
    w(1)=(T).5*x3-x2+(T)2/3;
    dw(1)=((T)1.5*x2-(T)2*x)*one_over_dX;
    x-=(T)1;
    x2=x*x;x3=x2*x;
    w(2)=-(T).5*x3-x2+(T)2/3;
    dw(2)=(-(T)1.5*x2-(T)2*x)*one_over_dX;
    x-=(T)1;
    z=(T)2+x;z2=z*z;
    w(3)=((T)1/6)*z2*z;
    dw(3)=(T).5*one_over_dX*z2;
}
//#####################################################################
// Function Update
//#####################################################################
template<class TV,int degree> void PARTICLE_GRID_WEIGHTS_SPLINE<TV,degree>::
Update(const ARRAY_VIEW<TV>& X)
{
    precompute_data.Resize(X.m);
    for(int p=0;p<X.m;p++){
        PRECOMPUTE_DATA& pd=precompute_data(p);
        pd.base=grid.Cell(X(p)-(T).5*degree*grid.DX(),degree+1);
        TV X_eval=X(p)-grid.Center(pd.base);
        Compute_Weights(pd.w,pd.dw,X_eval*grid.one_over_dX,grid.one_over_dX);}
}
//#####################################################################
// Function Constant_Scalar_Dp
//#####################################################################
template<class TV,int degree> typename TV::SCALAR PARTICLE_GRID_WEIGHTS_SPLINE<TV,degree>::
Constant_Scalar_Inverse_Dp() const
{
    PHYSBAM_ASSERT(degree>1);
    return (6-degree)*sqr(grid.one_over_dX(0));
}
//#####################################################################
// Function Dp
//#####################################################################
template<class TV,int degree> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> PARTICLE_GRID_WEIGHTS_SPLINE<TV,degree>::
Dp(const TV& X) const
{
    if(degree>1) return SYMMETRIC_MATRIX<T,TV::m>()+sqr(grid.dX(0))/(6-degree);
    TV Z=X-grid.Center(grid.Index(X));
    return DIAGONAL_MATRIX<T,TV::m>(Z*(grid.dX-Z));
}
//#####################################################################
// Function Order
//#####################################################################
template<class TV,int degree> int PARTICLE_GRID_WEIGHTS_SPLINE<TV,degree>::
Order() const
{
    return degree;
}
//#####################################################################
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,1>,1>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,2>,1>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,3>,1>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,1>,1>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,2>,1>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,3>,1>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,1>,2>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,2>,2>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,3>,2>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,1>,2>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,2>,2>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,3>,2>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,1>,3>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,2>,3>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<float,3>,3>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,1>,3>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,2>,3>;
template class PARTICLE_GRID_WEIGHTS_SPLINE<VECTOR<double,3>,3>;
}
