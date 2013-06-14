//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_PARTICLES
//#####################################################################
#ifndef __MPM_PARTICLES__
#define __MPM_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
#include "MPM_PARTICLES_FORWARD.h"

namespace PhysBAM{

template<class TV>
class MPM_PARTICLES:public CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> > // X, V, mass
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef CLONEABLE<MPM_PARTICLES<TV>,DEFORMABLE_PARTICLES<TV> > BASE;
public:
    using BASE::X;using BASE::V;using BASE::mass;using BASE::Store_Velocity;using BASE::Store_Mass;using BASE::Add_Array;using BASE::Remove_Array;using BASE::Get_Attribute_Index;using BASE::Remove_Array_Using_Index;using BASE::Resize;

    ARRAY_VIEW<TV> Xm;
    ARRAY_VIEW<T> volume;
    ARRAY_VIEW<VECTOR<TV,TV::m+1> > particle_domain;
    ARRAY_VIEW<MATRIX<T,TV::m> > Fe;
    ARRAY_VIEW<MATRIX<T,TV::m> > Fp;
    ARRAY_VIEW<T> mu,lambda,mu0,lambda0;
    ARRAY_VIEW<T> pressure;
    ARRAY_VIEW<bool> compress;
    ARRAY_VIEW<bool> use_plasticity_yield,use_plasticity_clamp;
    ARRAY_VIEW<T> yield_min,yield_max,clamp_min,clamp_max;
    ARRAY_VIEW<bool> use_visco_plasticity;
    ARRAY_VIEW<T> visco_nu,visco_tau,visco_kappa;

    RANDOM_NUMBERS<T> rand_generator;

    MPM_PARTICLES();
    virtual ~MPM_PARTICLES();

    template<class T_OBJECT> void Add_Randomly_Sampled_Object(const T_OBJECT& object,const T exclude_radius=(T)999)
    {Add_Randomly_Sampled_Implicit_Object(ANALYTIC_IMPLICIT_OBJECT<T_OBJECT>(object),exclude_radius);}
    void Add_Randomly_Sampled_Implicit_Object(const IMPLICIT_OBJECT<TV>& object,const T exclude_radius=(T)999);
    template<class T_OBJECT> void Delete_Particles_In_Object_Lazy(const T_OBJECT& object)
    {
        ARRAY<TV> new_Xm;
        for(int i=0;i<this->number;i++) if(!object.Lazy_Inside(Xm(i))) new_Xm.Append(Xm(i));
        this->Resize(new_Xm.m);
        X=Xm=new_Xm;
    }
    void Set_Material_Properties(int start_index,int count,T mass_in,T mu_in,T lambda_in,bool compress_in,bool pressure_in);
    void Set_Initial_State(int start_index,int count,MATRIX<T,TV::m> Fe_in,MATRIX<T,TV::m> Fp_in,TV V_in);
    void Set_Plasticity(int start_index,int count,bool use_plasticity_yield_in,T yield_min_in,T yield_max_in,bool use_plasticity_clamp_in,T clamp_min_in,T clamp_max_in);
    void Set_Visco_Plasticity(int start_index,int count,bool use_visco_plasticity_in,T visco_nu_in,T visco_tau_in,T visco_kappa_in);
//#####################################################################
};
}
#endif

