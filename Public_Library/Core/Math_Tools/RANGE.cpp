//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Duc Nguyen, Avi Robinson-Mosher, Craig Schroeder,
//      Andrew Selle, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Vectors/VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV RANGE<TV>::
Normal(const int aggregate) const
{
    assert(aggregate>=0 && aggregate<d*2);
    TV direction=TV::Axis_Vector(aggregate>>1);
    return (aggregate&1)?direction:-direction;
}
//#####################################################################
// Function Surface
//#####################################################################
template<class TV> TV RANGE<TV>::
Surface(const TV& X) const
{
    if(!Lazy_Inside(X)) return clamp(X,min_corner,max_corner);
    TV min_side=X-min_corner,max_side=max_corner-X,result=X;
    int min_argmin=min_side.Arg_Min(),max_argmin=max_side.Arg_Min();
    if(min_side(min_argmin)<=max_side(max_argmin)) result(min_argmin)=min_corner(min_argmin);
    else result(max_argmin)=max_corner(max_argmin);
    return result;
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class TV> typename TV::SCALAR RANGE<TV>::
Signed_Distance(const TV& X) const
{
    TV lengths=Edge_Lengths();
    if(lengths.Contains(0)) PHYSBAM_NOT_IMPLEMENTED();
    TV phi=abs(X-Center())-(T).5*lengths;
    if(!phi.All_Less_Equal(TV())) return TV::Componentwise_Max(phi,TV()).Magnitude();
    return phi.Max();
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV RANGE<TV>::
Normal(const TV& X) const
{
    if(Lazy_Inside(X)){
        TV phis_max=X-max_corner,phis_min=min_corner-X;
        int axis=TV::Componentwise_Max(phis_min,phis_max).Arg_Max();
        return T(phis_max[axis]>phis_min[axis]?1:-1)*TV::Axis_Vector(axis);}
    else{
        TV phis_max=X-max_corner,phis_min=min_corner-X;
        TV normal;
        for(int i=0;i<d;i++){
            T phi=max(phis_min(i),phis_max(i));
            normal(i)=phi>0?(phis_max(i)>phis_min(i)?phi:-phi):0;}
        return normal.Normalized();}
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> RANGE<TV>::
Hessian(const TV& X) const
{
    if(Lazy_Inside(X)) return SYMMETRIC_MATRIX<T,TV::m>();
    TV phis_max=X-max_corner,phis_min=min_corner-X;
    TV normal,D_normal_diag;

    for(int i=0;i<d;i++){
        T phi=max(phis_min(i),phis_max(i));
        if(phi>0){
            normal(i)=phis_max(i)>phis_min(i)?phi:-phi;
            D_normal_diag(i)=1;}}
    T norm=normal.Normalize();
    return (DIAGONAL_MATRIX<T,TV::m>(D_normal_diag)-SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(normal))/norm;
}
//#####################################################################
// Function Name
//#####################################################################
template<class TV> std::string RANGE<TV>::
Name()
{
    return LOG::sprintf("RANGE<VECTOR<T,%d>",d);
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template std::string RANGE<VECTOR<T,d> >::Name(); \
    template VECTOR<T,d> RANGE<VECTOR<T,d> >::Normal(const int) const; \
    template VECTOR<T,d> RANGE<VECTOR<T,d> >::Normal(const VECTOR<T,d>&) const; \
    template VECTOR<T,d> RANGE<VECTOR<T,d> >::Surface(const VECTOR<T,d>&) const; \
    template SYMMETRIC_MATRIX<T,d> RANGE<VECTOR<T,d> >::Hessian(const VECTOR<T,d>&) const; \
    template VECTOR<T,d>::SCALAR RANGE<VECTOR<T,d> >::Signed_Distance(const VECTOR<T,d>&) const;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
