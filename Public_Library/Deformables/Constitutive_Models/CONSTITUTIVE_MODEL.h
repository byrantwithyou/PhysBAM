//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis, Russell Howes, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTITUTIVE_MODEL
//##################################################################### 
#ifndef __CONSTITUTIVE_MODEL__
#define __CONSTITUTIVE_MODEL__

#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Deformables/Constitutive_Models/CONSTITUTIVE_MODELS_FORWARD.h>
namespace PhysBAM{

template<class TV> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;
template<class T,int d>
class CONSTITUTIVE_MODEL
{
    typedef VECTOR<T,d> TV;
public:
    bool enforce_definiteness;
    T constant_lambda,constant_mu; // Lame coefficients (used by almost all derived models)
    T constant_alpha,constant_beta; // isotropic damping parameters (used by all current derived models)
    const ARRAY_VIEW<T> *lambda,*mu; // spatially varying Lame coefficients
    const ARRAY_VIEW<T> *alpha,*beta; // spatially varying damping parameters

private:
    CONSTITUTIVE_MODEL();
    CONSTITUTIVE_MODEL(const CONSTITUTIVE_MODEL&) = delete;
    void operator=(const CONSTITUTIVE_MODEL&) = delete;

    // all constitutive models should derive from one of these
    friend class ISOTROPIC_CONSTITUTIVE_MODEL<T,d>;
    friend class ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>;
public:

    virtual ~CONSTITUTIVE_MODEL();

    virtual T Maximum_Elastic_Stiffness(const int id) const; // for elastic CFL computation
    virtual T Maximum_Damping_Stiffness(const int id) const; // for damping CFL computation

    T Mu(const int id) const {return mu?(*mu)(id):constant_mu;}
    T Lambda(const int id) const {return lambda?(*lambda)(id):constant_lambda;}
    T Alpha(const int id) const {return alpha?(*alpha)(id):constant_alpha;}
    T Beta(const int id) const {return beta?(*beta)(id):constant_beta;}

//#####################################################################
    virtual MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const=0;
    virtual void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dPi_dF,const int id) const;
    virtual int P_From_Strain_Rate_Forces_Size() const;
    virtual void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const;
    virtual MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const int id) const;
    virtual void Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input);
//#####################################################################
};
}
#endif
