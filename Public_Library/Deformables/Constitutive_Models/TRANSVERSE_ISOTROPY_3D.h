//#####################################################################
// Copyright 2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRANSVERSE_ISOTROPY_3D
//##################################################################### 
#ifndef __TRANSVERSE_ISOTROPY_3D__
#define __TRANSVERSE_ISOTROPY_3D__

#include <Tools/Arrays/ARRAY.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T>
class TRANSVERSE_ISOTROPY_3D:public ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ANISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;
    using BASE::enforce_definiteness;using BASE::constant_alpha;using BASE::constant_beta;using BASE::use_isotropic_component_of_stress_derivative_only;
    
    T failure_threshold;
    ARRAY<VECTOR<T,3> > fiber_field;

    TRANSVERSE_ISOTROPY_3D(const T failure_threshold_input)
        :failure_threshold(failure_threshold_input)
    {
        use_isotropic_component_of_stress_derivative_only=false;
    }

    void Initialize_Fiber_Field_From_Current_State(const STRAIN_MEASURE<TV,3>& strain_measure,const ARRAY<VECTOR<T,3> >& fiber_field_input)
    {int n=strain_measure.tetrahedron_mesh.tetrahedrons.m;fiber_field.Resize(n);for(int t=0;t<n;t++) fiber_field(t)=strain_measure.F(t).Transpose_Times(fiber_field_input(t));}

    inline int Hessian_Index(const int m,const int n) const
    {assert(m>=n);return m+(n-1)*(10-n)/2;}

    void Invariants(ARRAY<T>& invariants,const DIAGONAL_MATRIX<T,3>& C,const VECTOR<T,3> V_fiber) const
    {invariants(0)=C.Trace();invariants(1)=(C*C).Trace();invariants(2)=C.Determinant();
    invariants(3)=VECTOR<T,3>::Dot_Product(V_fiber,C*V_fiber);invariants(4)=(C*V_fiber).Magnitude_Squared();}

    SYMMETRIC_MATRIX<T,3> S(const DIAGONAL_MATRIX<T,3>& C,const VECTOR<T,3>& V_fiber,const ARRAY<T>& invariants,const ARRAY<T> energy_gradient) const
    {SYMMETRIC_MATRIX<T,3> result;
    if(energy_gradient(0)) result+=(T)2*energy_gradient(0);
    if(energy_gradient(1)) result+=(T)4*energy_gradient(1)*C;
    if(energy_gradient(2)) result+=(T)2*energy_gradient(2)*invariants(2)*C.Inverse();
    if(energy_gradient(3)) result+=(T)2*energy_gradient(3)*SYMMETRIC_MATRIX<T,3>::Outer_Product(V_fiber);
    if(energy_gradient(4)) result+=(T)2*energy_gradient(4)*(C*SYMMETRIC_MATRIX<T,3>::Outer_Product(V_fiber)).Twice_Symmetric_Part();
    return result;}

    MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold),C=F_threshold*F_threshold;VECTOR<T,3> V_fiber=V.Transpose_Times(fiber_field(tetrahedron));
    ARRAY<T> invariants(4),energy_gradient(4);SYMMETRIC_MATRIX<T,3> result;Invariants(invariants,C,V_fiber);Energy_Gradient(energy_gradient,invariants,tetrahedron);
    return scale*F_threshold*S(C,V_fiber,invariants,energy_gradient);}

    MATRIX<T,3> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {SYMMETRIC_MATRIX<T,3> strain_rate=F_dot.Symmetric_Part(); // Use linear damping by default
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();}

    void Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,DIAGONALIZED_STRESS_DERIVATIVE<T,3>& dP_dF,const int tetrahedron=0) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold),C=F_threshold*F_threshold,C_inverse=C.Inverse();VECTOR<T,3> V_fiber=V.Transpose_Times(fiber_field(tetrahedron));
    SYMMETRIC_MATRIX<T,3> C_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR<T,3>(C_inverse.x00,C_inverse.x11,C_inverse.x22));
    SYMMETRIC_MATRIX<T,3> V_fiber_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(V_fiber);
    ARRAY<T> invariants(4),energy_gradient(4),energy_hessian(15);Invariants(invariants,C,V_fiber);
    Energy_Gradient(energy_gradient,invariants,tetrahedron);Energy_Hessian(energy_hessian,invariants,tetrahedron);
    ARRAY<VECTOR<T,3> > J_d(4),J_s(4);
    J_d(0)=VECTOR<T,3>((T)1,(T)1,(T)1);
    J_d(1)=(T)2*VECTOR<T,3>(C.x00,C.x11,C.x22);
    J_d(2)=invariants(2)*VECTOR<T,3>(C_inverse.x00,C_inverse.x11,C_inverse.x22);
    J_d(3)=VECTOR<T,3>(V_fiber_outer.x00,V_fiber_outer.x11,V_fiber_outer.x22);
    J_d(4)=(T)2*VECTOR<T,3>(C.x00*V_fiber_outer.x00,C.x11*V_fiber_outer.x11,C.x22*V_fiber_outer.x22);
    J_s(3)=(T)root_two*VECTOR<T,3>(V_fiber_outer.x10,V_fiber_outer.x20,V_fiber_outer.x21);
    J_s(4)=(T)root_two*VECTOR<T,3>((C.x00+C.x11)*V_fiber_outer.x10,(C.x00+C.x22)*V_fiber_outer.x20,(C.x11+C.x22)*V_fiber_outer.x21);
    dP_dF.dSdC_d=dP_dF.dSdC_s=SYMMETRIC_MATRIX<T,3>();dP_dF.dSdC_ds=MATRIX<T,3>();dP_dF.F=F_threshold;dP_dF.S=S(C,V_fiber,invariants,energy_gradient);
    for(int n=0;n<5;n++) for(int m=n;m<=5;m++){
        int hessian_index=Hessian_Index(m,n);if(!energy_hessian(hessian_index))continue;
        if(m==n){
            dP_dF.dSdC_d+=(T)2*energy_hessian(hessian_index)*SYMMETRIC_MATRIX<T,3>::Outer_Product(J_d(m));
            dP_dF.dSdC_s+=(T)2*energy_hessian(hessian_index)*SYMMETRIC_MATRIX<T,3>::Outer_Product(J_s(m));
            dP_dF.dSdC_ds+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_d(m),J_s(m));}
        else{
            dP_dF.dSdC_d+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_d(m),J_d(n)).Twice_Symmetric_Part();
            dP_dF.dSdC_s+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_s(m),J_s(n)).Twice_Symmetric_Part();
            dP_dF.dSdC_ds+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_d(m),J_s(n)).Twice_Symmetric_Part();}}
    if(energy_gradient(1)){
        dP_dF.dSdC_d+=(T)4*energy_gradient(1);
        dP_dF.dSdC_s+=(T)4*energy_gradient(1);}
    if(energy_gradient(2)){
        dP_dF.dSdC_d+=(T)2*energy_gradient(2)*invariants(2)*SYMMETRIC_MATRIX<T,3>(0,C_inverse_outer.x10,C_inverse_outer.x20,0,C_inverse_outer.x21,0);
        dP_dF.dSdC_s-=(T)2*energy_gradient(2)*invariants(2)*DIAGONAL_MATRIX<T,3>(C_inverse_outer.x10,C_inverse_outer.x20,C_inverse_outer.x21);}
    if(energy_gradient(4)){
        dP_dF.dSdC_d+=(T)4*energy_gradient(4)*DIAGONAL_MATRIX<T,3>(V_fiber_outer.x00,V_fiber_outer.x11,V_fiber_outer.x22);
        dP_dF.dSdC_s+=(T)2*energy_gradient(4)*SYMMETRIC_MATRIX<T,3>(V_fiber_outer.x00+V_fiber_outer.x11,V_fiber_outer.x21,V_fiber_outer.x20,
                                                                      V_fiber_outer.x00+V_fiber_outer.x22,V_fiber_outer.x10,V_fiber_outer.x11+V_fiber_outer.x22);
        dP_dF.dSdC_ds+=(T)2*(T)root_two*energy_gradient(4)*MATRIX<T,3>(V_fiber_outer.x10,V_fiber_outer.x10,0,V_fiber_outer.x20,0,V_fiber_outer.x20,0,V_fiber_outer.x21,V_fiber_outer.x21);}
    if(enforce_definiteness)dP_dF.Enforce_Definiteness();}

//#####################################################################
    virtual void Energy_Gradient(ARRAY<T>& energy_gradient,const ARRAY<T>& invariants,const int tetrahedron) const=0;
    virtual void Energy_Hessian(ARRAY<T>& energy_hessian,const ARRAY<T>& invariants,const int tetrahedron) const=0;
//#####################################################################
};
}
#endif
