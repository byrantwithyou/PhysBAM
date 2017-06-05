//#####################################################################
// Copyright 2002-2006, Ron Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/constants.h>
#include <Grid_Tools/Fourier_Transforms/FFT_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Function Symmetry_Helper
//#####################################################################
template<class T> static void
Symmetry_Helper(std::complex<T>& a,std::complex<T>& b)
{
    std::complex<T> c=(T).5*(a+conj(b));
    a=c;
    b=conj(c);
}
//#####################################################################
// Function Enforce_Real_Valued_Symmetry
//#####################################################################
// enforce symmetry so that the inverse transform is a real valued function
template<class T,int d> void
Enforce_Real_Valued_Symmetry(ARRAY<std::complex<T>,VECTOR<int,d> >& u_hat)
{
    typedef VECTOR<int,d> TV_INT;
    TV_INT counts=u_hat.domain.Edge_Lengths();
    for(RANGE_ITERATOR<d> it(counts);it.Valid();it.Next()){
        TV_INT b=-it.index;
        for(int i=0;i<d;i++) if(b(i)<0) b(i)+=counts(i);
        Symmetry_Helper(u_hat(it.index),u_hat(b));}
}
//#####################################################################
// Function First_Derivatives
//#####################################################################
template<class T,int d> void
First_Derivatives(const ARRAY<std::complex<T>,VECTOR<int,d> >& u_hat,
    VECTOR<ARRAY<std::complex<T>,VECTOR<int,d> >,d>& du_hat,
    const VECTOR<T,d>& domain_size)
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    TV_INT counts=u_hat.domain.Edge_Lengths(),hi=counts/2,lo=hi-counts;
    TV coefficients=(T)(2*pi)/domain_size;
    for(RANGE_ITERATOR<d> it(counts);it.Valid();it.Next()){
        TV k=coefficients*TV(wrap(it.index,lo,hi));
        std::complex<T> w=u_hat(it.index)*std::complex<T>(0,1);
        for(int i=0;i<d;i++) du_hat(i)(it.index)=w*k(i);}
    for(int i=0;i<d;i++) Enforce_Real_Valued_Symmetry(du_hat(i));
}
//#####################################################################
// Make_Divergence_Free
//#####################################################################
// the constant fields are unchanged
template<class T,int d> void
Make_Divergence_Free(VECTOR<ARRAY<std::complex<T>,VECTOR<int,d> >,d>& u_hat,
    const VECTOR<T,d>& domain_size)
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    TV_INT counts=u_hat(0).domain.Edge_Lengths(),hi=counts/2,lo=hi-counts;
    TV coefficients=(T)(2*pi)/domain_size;
    for(RANGE_ITERATOR<d> it(counts);it.Valid();it.Next()){
        if(it.index==TV_INT()) continue;
        TV k=coefficients*TV(wrap(it.index,lo,hi));
        std::complex<T> w=0;
        for(int i=0;i<d;i++) w+=k(i)*u_hat(i)(it.index);
        w/=k.Magnitude_Squared();
        for(int i=0;i<d;i++) u_hat(i)(it.index)-=w*k(i);}
    for(int i=0;i<d;i++) Enforce_Real_Valued_Symmetry(u_hat(i));
}
//#####################################################################
// Function Filter_High_Frequencies
//#####################################################################
// Lanczos filter - doesn't change (0,0) frequency
template<class T,int d> void
Filter_High_Frequencies(ARRAY<std::complex<T>,VECTOR<int,d> >& u_hat,T scale)
{
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    TV_INT counts=u_hat.domain.Edge_Lengths(),hi=counts/2,lo=hi-counts;
    T coefficient=2*(T)pi/sqrt((T)counts.Magnitude_Squared())*scale;
    for(RANGE_ITERATOR<d> it(counts);it.Valid();it.Next()){
        TV freq(wrap(it.index,lo,hi));
        u_hat(it.index)*=sinc(coefficient*freq.Magnitude());}
    Enforce_Real_Valued_Symmetry(u_hat);
}
//#####################################################################
template void Enforce_Real_Valued_Symmetry<double,1>(ARRAY<std::complex<double>,VECTOR<int,1> >&);
template void Enforce_Real_Valued_Symmetry<double,2>(ARRAY<std::complex<double>,VECTOR<int,2> >&);
template void Enforce_Real_Valued_Symmetry<double,3>(ARRAY<std::complex<double>,VECTOR<int,3> >&);
template void Enforce_Real_Valued_Symmetry<float,1>(ARRAY<std::complex<float>,VECTOR<int,1> >&);
template void Enforce_Real_Valued_Symmetry<float,2>(ARRAY<std::complex<float>,VECTOR<int,2> >&);
template void Enforce_Real_Valued_Symmetry<float,3>(ARRAY<std::complex<float>,VECTOR<int,3> >&);
template void Make_Divergence_Free<double,1>(VECTOR<ARRAY<std::complex<double>,VECTOR<int,1> >,1>&,VECTOR<double,1> const&);
template void Make_Divergence_Free<double,2>(VECTOR<ARRAY<std::complex<double>,VECTOR<int,2> >,2>&,VECTOR<double,2> const&);
template void Make_Divergence_Free<double,3>(VECTOR<ARRAY<std::complex<double>,VECTOR<int,3> >,3>&,VECTOR<double,3> const&);
template void Make_Divergence_Free<float,1>(VECTOR<ARRAY<std::complex<float>,VECTOR<int,1> >,1>&,VECTOR<float,1> const&);
template void Make_Divergence_Free<float,2>(VECTOR<ARRAY<std::complex<float>,VECTOR<int,2> >,2>&,VECTOR<float,2> const&);
template void Make_Divergence_Free<float,3>(VECTOR<ARRAY<std::complex<float>,VECTOR<int,3> >,3>&,VECTOR<float,3> const&);
template void First_Derivatives<double,1>(ARRAY<std::complex<double>,VECTOR<int,1> > const&,VECTOR<ARRAY<std::complex<double>,VECTOR<int,1> >,1>&,VECTOR<double,1> const&);
template void First_Derivatives<double,2>(ARRAY<std::complex<double>,VECTOR<int,2> > const&,VECTOR<ARRAY<std::complex<double>,VECTOR<int,2> >,2>&,VECTOR<double,2> const&);
template void First_Derivatives<double,3>(ARRAY<std::complex<double>,VECTOR<int,3> > const&,VECTOR<ARRAY<std::complex<double>,VECTOR<int,3> >,3>&,VECTOR<double,3> const&);
template void First_Derivatives<float,1>(ARRAY<std::complex<float>,VECTOR<int,1> > const&,VECTOR<ARRAY<std::complex<float>,VECTOR<int,1> >,1>&,VECTOR<float,1> const&);
template void First_Derivatives<float,2>(ARRAY<std::complex<float>,VECTOR<int,2> > const&,VECTOR<ARRAY<std::complex<float>,VECTOR<int,2> >,2>&,VECTOR<float,2> const&);
template void First_Derivatives<float,3>(ARRAY<std::complex<float>,VECTOR<int,3> > const&,VECTOR<ARRAY<std::complex<float>,VECTOR<int,3> >,3>&,VECTOR<float,3> const&);
}
