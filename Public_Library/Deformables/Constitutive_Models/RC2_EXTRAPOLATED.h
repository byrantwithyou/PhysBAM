//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RC2_EXTRAPOLATED
//#####################################################################
#ifndef __RC2_EXTRAPOLATED__
#define __RC2_EXTRAPOLATED__

#include <Deformables/Constitutive_Models/GENERAL_ENERGY.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

using ::std::log;

template<class T,int d> class DIAGONAL_MATRIX;
template<class T,int d>
class RC2_EXTRAPOLATED:public ISOTROPIC_CONSTITUTIVE_MODEL<T,d>
{
    typedef VECTOR<T,d> TV;

public:

    typedef ISOTROPIC_CONSTITUTIVE_MODEL<T,d> BASE;
    
    using BASE::enforce_definiteness;
    using BASE::constant_lambda;using BASE::constant_mu;using BASE::constant_alpha;using BASE::constant_beta;
    using BASE::Alpha;using BASE::Beta;using BASE::Lambda;using BASE::Mu;

    struct HELPER
    {
        TV q,u,g,Hu,uTu;
        MATRIX<T,d> dq,du,dg;
        VECTOR<SYMMETRIC_MATRIX<T,d>,d> ddq,ddu,ddg,TT;
        VECTOR<VECTOR<SYMMETRIC_MATRIX<T,d>,d>,d> A;
        T m,h,phi,E,z,xi,s,b,c;
        TV dm,dh,dphi,dE,dz,dxi,ds,db,dc;
        SYMMETRIC_MATRIX<T,d> ddm,ddh,ddphi,ddE,ddz,ddxi,dds,H,Tu,ddb,ddc;

        bool Compute_E(const GENERAL_ENERGY<T>& base,T extrapolation_cutoff,const TV& f,const int id);
        void Compute_dE(const GENERAL_ENERGY<T>& base,const TV& f,const int id);
        void Compute_ddE(const GENERAL_ENERGY<T>& base,const TV& f,const int id);
        VECTOR<T,TV::SPIN::m> g_it,dE_it,u_it,H_iitt,uTu_it,T_iiittt,T_itit;
        VECTOR<T,(TV::m==3)*3> T_xiitt,T_xxit,H_xit;
    };

    GENERAL_ENERGY<T>& base;

    T youngs_modulus,poissons_ratio;
    T extrapolation_cutoff;
    T panic_threshold;

    RC2_EXTRAPOLATED(GENERAL_ENERGY<T>& ge_input,const T youngs_modulus_input=3e6,const T poissons_ratio_input=.475,const T Rayleigh_coefficient=.05,
        const T extrapolation_cutoff=.3,const T extra_force_coefficient=20);
    virtual ~RC2_EXTRAPOLATED();

    void Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input) override;
    T Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const override;
    DIAGONAL_MATRIX<T,d> P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const override;
    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dP_dF,const int id) const override;
    
    MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const override;
    int P_From_Strain_Rate_Forces_Size() const override;
    void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const override;
    MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const int id) const override;
    void Test_Model() const;
};
}
#endif
