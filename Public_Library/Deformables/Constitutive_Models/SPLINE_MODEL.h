//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPLINE_MODEL
//#####################################################################
#ifndef __SPLINE_MODEL__
#define __SPLINE_MODEL__

#include <Tools/Math_Tools/constants.h>
#include <Deformables/Constitutive_Models/COROTATED.h>
namespace PhysBAM{

template<class T,int d>
class SPLINE_MODEL:public COROTATED<T,d>
{
    typedef VECTOR<T,d> TV;
public:
    typedef COROTATED<T,d> BASE;
    using BASE::constant_lambda;using BASE::constant_mu;using BASE::lambda;using BASE::mu;

private:
    ARRAY<T> hardening_deformation,hardening_strength;
    ARRAY<T> coefficient,base;
public:

    template<class T_FIELD0,class T_FIELD1>
    SPLINE_MODEL(const T_FIELD0& youngs_modulus_input=3e6,const T_FIELD0& poissons_ratio_input=.475,const T_FIELD1& hardening_deformation_input=.5,
        const T_FIELD1& hardening_strength_input=7,const T_FIELD0& Rayleigh_coefficient=.05)
    {
        Set_Parameters(youngs_modulus_input,poissons_ratio_input,hardening_deformation_input,hardening_strength_input,Rayleigh_coefficient);
    }

    template<class T_FIELD0,class T_FIELD1>
    void Set_Parameters(const T_FIELD0& youngs_modulus_input,const T_FIELD0& poissons_ratio_input,const T_FIELD1& hardening_deformation_input,const T_FIELD1& hardening_strength_input,
        const T_FIELD0& Rayleigh_coefficient)
    {COROTATED<T,d>::Set_Parameters(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
    Set(hardening_deformation,hardening_deformation_input);
    Set(hardening_strength,hardening_strength_input);
    coefficient=(hardening_strength-1)/(3*(hardening_deformation*hardening_deformation));
    base=hardening_deformation*(hardening_strength-1)*((T)2/3);}

private:
    static void Set(ARRAY<T>& parameter,const T value)
    {parameter.Resize(1);parameter(0)=value;}

    static void Set(ARRAY<T>& parameter,ARRAY_VIEW<const T> value)
    {parameter=value;}
public:

    DIAGONAL_MATRIX<T,2> P_From_Strain(const DIAGONAL_MATRIX<T,2>& F,const T scale,const int simplex) const
    {DIAGONAL_MATRIX<T,2> strain=F-1,strain_abs=strain.Abs(),strain_sign=strain.Sign();
    DIAGONAL_MATRIX<T,2> D;
    int index=hardening_deformation.m>1?simplex:1;
    T hardening_deformation_=hardening_deformation(index),hardening_strength_=hardening_strength(index),coefficient_=coefficient(index),base_=base(index);
    D.x.x=strain_abs.x.x<hardening_deformation_?strain_abs.x.x*(1+coefficient_*sqr(strain_abs.x.x)):hardening_strength_*strain_abs.x.x-base_;
    D.x.y=strain_abs.x.y<hardening_deformation_?strain_abs.x.y*(1+coefficient_*sqr(strain_abs.x.y)):hardening_strength_*strain_abs.x.y-base_;
    if(!mu.m) return 2*scale*constant_mu*strain_sign*D+scale*constant_lambda*strain.Trace();
    else return 2*scale*mu(simplex)*strain_sign*D+scale*lambda(simplex)*strain.Trace();}

    DIAGONAL_MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const T scale,const int simplex) const
    {DIAGONAL_MATRIX<T,3> strain=F-1,strain_abs=strain.Abs(),strain_sign=strain.Sign();
    DIAGONAL_MATRIX<T,3> D;
    int index=hardening_deformation.m>1?simplex:1;
    T hardening_deformation_=hardening_deformation(index),hardening_strength_=hardening_strength(index),coefficient_=coefficient(index),base_=base(index);
    D.x.x=strain_abs.x.x<hardening_deformation_?strain_abs.x.x*(1+coefficient_*sqr(strain_abs.x.x)):hardening_strength_*strain_abs.x.x-base_;
    D.x.y=strain_abs.x.y<hardening_deformation_?strain_abs.x.y*(1+coefficient_*sqr(strain_abs.x.y)):hardening_strength_*strain_abs.x.y-base_;
    D.x.z=strain_abs.x.z<hardening_deformation_?strain_abs.x.z*(1+coefficient_*sqr(strain_abs.x.z)):hardening_strength_*strain_abs.x.z-base_;
    if(!mu.m) return 2*scale*constant_mu*strain_sign*D+scale*constant_lambda*strain.Trace();
    else return 2*scale*mu(simplex)*strain_sign*D+scale*lambda(simplex)*strain.Trace();}

//#####################################################################
};
}
#endif
