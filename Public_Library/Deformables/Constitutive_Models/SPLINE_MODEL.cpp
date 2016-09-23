//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Deformables/Constitutive_Models/SPLINE_MODEL.h>
namespace PhysBAM{
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> SPLINE_MODEL<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T id_mu=Mu(id),id_lambda=Lambda(id);
    DIAGONAL_MATRIX<T,d> strain=F-1,strain_abs=strain.Abs(),strain_sign=strain.Sign();
    DIAGONAL_MATRIX<T,d> D;
    int index=hardening_deformation.m>1?id:1;
    T hardening_deformation_=hardening_deformation(index),hardening_strength_=hardening_strength(index),coefficient_=coefficient(index),base_=base(index);
    for(int i=0;i<d;i++)
        D.x(i)=strain_abs.x(i)<hardening_deformation_?strain_abs.x(i)*(1+coefficient_*sqr(strain_abs.x(i))):hardening_strength_*strain_abs.x(i)-base_;
    return 2*id_mu*strain_sign*D+id_lambda*strain.Trace();
}
template class SPLINE_MODEL<float,2>;
template class SPLINE_MODEL<float,3>;
template class SPLINE_MODEL<double,2>;
template class SPLINE_MODEL<double,3>;
}
