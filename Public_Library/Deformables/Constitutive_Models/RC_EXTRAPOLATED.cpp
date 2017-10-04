 //#####################################################################
// Copyright 2011, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/cube.h>
#include <Core/Math_Tools/pow.h>
#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX_2X2.h>
#include <Core/Matrices/MATRIX_3X3.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Core/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Tools/Polynomials/CUBIC.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/RC_EXTRAPOLATED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> RC_EXTRAPOLATED<T,d>::
RC_EXTRAPOLATED(GENERAL_ENERGY<T>& ge_input,const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input,
    const T extra_force_coefficient_input)
    :base(ge_input),youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),
    extrapolation_cutoff(extrapolation_cutoff_input),extra_force_coefficient(extra_force_coefficient_input),
    panic_threshold((T)1e-6)
{
    Update_Lame_Constants(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
    base.Initialize(constant_mu,constant_lambda);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> RC_EXTRAPOLATED<T,d>::
~RC_EXTRAPOLATED()
{
    delete &base;
}
//#####################################################################
// Update Lame Constants
//#####################################################################
template<class T,int d> void RC_EXTRAPOLATED<T,d>::
Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input)
{
    assert(poissons_ratio>-1&&poissons_ratio<.5);
    constant_lambda=youngs_modulus_input*poissons_ratio_input/((1+poissons_ratio_input)*(1-2*poissons_ratio_input));
    constant_mu=youngs_modulus_input/(2*(1+poissons_ratio_input));
    constant_alpha=Rayleigh_coefficient_input*constant_lambda;
    constant_beta=Rayleigh_coefficient_input*constant_mu;
    youngs_modulus=youngs_modulus_input; poissons_ratio=poissons_ratio_input;
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T RC_EXTRAPOLATED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T J=F.To_Vector().Product();
    if(J<extrapolation_cutoff){
        HELPER helper;
        bool b=helper.Compute_E(base,extra_force_coefficient*youngs_modulus,extrapolation_cutoff,F.To_Vector(),id);
        if(b) return helper.E;}
    return base.E(F.To_Vector(),id);
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> RC_EXTRAPOLATED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T J=F.To_Vector().Product();
    if(J<extrapolation_cutoff){
        HELPER helper;
        bool b=helper.Compute_E(base,extra_force_coefficient*youngs_modulus,extrapolation_cutoff,F.To_Vector(),id);
        if(b){
            helper.Compute_dE(base,extra_force_coefficient*youngs_modulus,F.To_Vector(),id);
//            Add_Debug_Particle(F.To_Vector(),VECTOR<T,3>(1,0,0));
//            Debug_Particle_Set_Attribute<TV>("V",-helper.dE/youngs_modulus);
            return DIAGONAL_MATRIX<T,d>(helper.dE);}}
//    Add_Debug_Particle(F.To_Vector(),VECTOR<T,3>(0,1,0));
//    Debug_Particle_Set_Attribute<TV>("V",-base.dE(F.To_Vector(),id)/youngs_modulus);
    return DIAGONAL_MATRIX<T,d>(base.dE(F.To_Vector(),id));
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const RC_EXTRAPOLATED<T,2>& re,const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,2> >& dP_dF,const int id)
{
    T J=F.To_Vector().Product();
    if(J<re.extrapolation_cutoff){
        typename RC_EXTRAPOLATED<T,2>::HELPER helper;
        bool b=helper.Compute_E(re.base,re.extra_force_coefficient*re.youngs_modulus,re.extrapolation_cutoff,F.To_Vector(),id);
        if(b){
            helper.Compute_dE(re.base,re.extra_force_coefficient*re.youngs_modulus,F.To_Vector(),id);
            helper.Compute_ddE(re.base,re.extra_force_coefficient*re.youngs_modulus,F.To_Vector(),id);
            dP_dF.H(0,0)=helper.ddE.x00;
            dP_dF.H(1,0)=helper.ddE.x10;
            dP_dF.H(1,1)=helper.ddE.x11;
            T ss1=sqr(F.x.x),ss2=sqr(F.x.y);
            T s01=ss1-ss2;
            if(fabs(s01)<re.panic_threshold) s01=s01<0?-re.panic_threshold:re.panic_threshold;
            dP_dF.C(0)=(-helper.dE.y*F.x.x+helper.dE.x*F.x.y)/s01;
            dP_dF.B(0)=(-helper.dE.y*F.x.y+helper.dE.x*F.x.x)/s01;
            return;}}
    T x = F.x.x, y = F.x.y, xpy = x+y;
    dP_dF.H(0,0) = re.base.Exx(x,y,id);
    dP_dF.H(1,0) = re.base.Exy(x,y,id);
    dP_dF.H(1,1) = re.base.Eyy(x,y,id);
    if(fabs(xpy)<re.panic_threshold) xpy=xpy<0?-re.panic_threshold:re.panic_threshold;
    T S=re.P_From_Strain(F,id).Trace()/xpy, D=re.base.Ex_Ey_x_y(x,y,id);
    dP_dF.C(0) = (D-S)/2;
    dP_dF.B(0) = (D+S)/2;
}
///#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const RC_EXTRAPOLATED<T,3>& re,const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<VECTOR<T,3> >& dP_dF,const int id)
{
    T J=F.To_Vector().Product();
    if(J<re.extrapolation_cutoff){
        typename RC_EXTRAPOLATED<T,3>::HELPER helper;
        bool b=helper.Compute_E(re.base,re.extra_force_coefficient*re.youngs_modulus,re.extrapolation_cutoff,F.To_Vector(),id);
        if(b){
            helper.Compute_dE(re.base,re.extra_force_coefficient*re.youngs_modulus,F.To_Vector(),id);
            helper.Compute_ddE(re.base,re.extra_force_coefficient*re.youngs_modulus,F.To_Vector(),id);
            dP_dF.H(0,0)=helper.ddE.x00;
            dP_dF.H(1,1)=helper.ddE.x11;
            dP_dF.H(2,2)=helper.ddE.x22;
            dP_dF.H(1,0)=helper.ddE.x10;
            dP_dF.H(2,0)=helper.ddE.x20;
            dP_dF.H(2,1)=helper.ddE.x21;
            T ss1=sqr(F.x.x),ss2=sqr(F.x.y),ss3=sqr(F.x.z);
            T s01=ss1-ss2,s02=ss1-ss3,s12=ss2-ss3;
            if(fabs(s01)<re.panic_threshold) s01=s01<0?-re.panic_threshold:re.panic_threshold;
            if(fabs(s02)<re.panic_threshold) s02=s02<0?-re.panic_threshold:re.panic_threshold;
            if(fabs(s12)<re.panic_threshold) s12=s12<0?-re.panic_threshold:re.panic_threshold;
            dP_dF.C(2)=(-helper.dE.y*F.x.x+helper.dE.x*F.x.y)/s01;
            dP_dF.B(2)=(-helper.dE.y*F.x.y+helper.dE.x*F.x.x)/s01;
            dP_dF.C(1)=(-helper.dE.z*F.x.x+helper.dE.x*F.x.z)/s02;
            dP_dF.B(1)=(-helper.dE.z*F.x.z+helper.dE.x*F.x.x)/s02;
            dP_dF.C(0)=(-helper.dE.z*F.x.y+helper.dE.y*F.x.z)/s12;
            dP_dF.B(0)=(-helper.dE.z*F.x.z+helper.dE.y*F.x.y)/s12;
            return;}}
    T x = F.x.x, y = F.x.y, z = F.x.z, xpy = x+y, xpz = x+z, ypz = y+z;
    dP_dF.H(0,0) = re.base.Exx(x,y,z,id);
    dP_dF.H(1,1) = re.base.Eyy(x,y,z,id);
    dP_dF.H(2,2) = re.base.Ezz(x,y,z,id);
    dP_dF.H(1,0) = re.base.Exy(x,y,z,id);
    dP_dF.H(2,0) = re.base.Exz(x,y,z,id);
    dP_dF.H(2,1) = re.base.Eyz(x,y,z,id);
    if(fabs(xpy)<re.panic_threshold) xpy=xpy<0?-re.panic_threshold:re.panic_threshold;
    VECTOR<T,3> P=re.P_From_Strain(F,id).To_Vector();
    T S01=(P.x+P.y)/xpy, D01=re.base.Ex_Ey_x_y(x,y,z,id);
    dP_dF.C(2) = (D01-S01)/2;
    dP_dF.B(2) = (D01+S01)/2;
    if(fabs(xpz)<re.panic_threshold) xpz=xpz<0?-re.panic_threshold:re.panic_threshold;
    T S02=(P.x+P.z)/xpz, D02=re.base.Ex_Ez_x_z(x,y,z,id);
    dP_dF.C(1) = (D02-S02)/2;
    dP_dF.B(1) = (D02+S02)/2;
    if(fabs(ypz)<re.panic_threshold) ypz=ypz<0?-re.panic_threshold:re.panic_threshold;
    T S12=(P.y+P.z)/ypz, D12=re.base.Ey_Ez_y_z(x,y,z,id);
    dP_dF.C(0) = (D12-S12)/2;
    dP_dF.B(0) = (D12+S12)/2;
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void RC_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<TV>& dP_dF,const int id) const
{
    return Isotropic_Stress_Derivative_Helper(*this,F,dP_dF,id);
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> RC_EXTRAPOLATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int RC_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void RC_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::m;
    T sa=sqrt(id_alpha/TV::m+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> RC_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const int id) const
{
    T id_alpha=Alpha(id),id_beta=Beta(id);
    SYMMETRIC_MATRIX<T,d> strain_rate=(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::m;
    T sa=sqrt(id_alpha/TV::m+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Compute_s
//#####################################################################
template<class T> static T
Compute_s(const VECTOR<T,2>& f,const int id,T extrapolation_cutoff)
{
    VECTOR<T,2> fm1=f-1;
    QUADRATIC<T> quadratic(fm1.Product(),fm1.Sum(),1-extrapolation_cutoff);
    quadratic.Compute_Roots_In_Interval(0,1);
    if(quadratic.roots==0) return -1;
    PHYSBAM_ASSERT(quadratic.roots==1);
    return quadratic.root1;
}
//#####################################################################
// Function Compute_s
//#####################################################################
template<class T> static T
Compute_s(const VECTOR<T,3>& f,const int id,T extrapolation_cutoff)
{
    VECTOR<T,3> fm1=f-1;
    CUBIC<T> cubic(fm1.Product(),DIAGONAL_MATRIX<T,3>(fm1).Cofactor_Matrix().Trace(),fm1.Sum(),1-extrapolation_cutoff);
    cubic.Compute_Roots_In_Interval(0,1);
    if(cubic.roots==0) return -1;
    PHYSBAM_ASSERT(cubic.roots==1);
    return cubic.root1;
}
//#####################################################################
// Function Compute_E
//#####################################################################
template<class T,int d> bool RC_EXTRAPOLATED<T,d>::HELPER::
Compute_E(const GENERAL_ENERGY<T>& base,T k,T extrapolation_cutoff,const TV& f,const int id)
{
    TV fm1=f-1;
    s=Compute_s(f,id,extrapolation_cutoff);
    if(s==-1) return false;
    Q=fm1*s+1;
    z=(fm1/Q).Sum();
    xi=1/z;
    phi=base.E(Q,id);
    m=1/fm1.Magnitude();
    u=m*fm1;
    g=base.dE(Q,id);
    h=TV::Dot_Product(f-Q,u);
    E=phi+h*TV::Dot_Product(g,u)+(T).5*k*sqr(h);
    return true;
}
//#####################################################################
// Function Compute_dE
//#####################################################################
template<class T,int d> void RC_EXTRAPOLATED<T,d>::HELPER::
Compute_dE(const GENERAL_ENERGY<T>& base,T k,const TV& f,const int id)
{
    TV fm1=f-1;
    ds=-s*xi/Q;
    dQ=s+MATRIX<T,d>::Outer_Product(fm1,ds);
    dz=(T)1/Q-dQ.Transpose_Times(fm1/sqr(Q));
    dxi=-sqr(xi)*dz;
    H=base.ddE(Q,id);
    dphi=dQ.Transpose_Times(g);
    dm=-cube(m)*fm1;
    du=MATRIX<T,d>::Outer_Product(fm1,dm)+m;
    dg=H*dQ;
    dh=((T)1-dQ).Transpose_Times(u)+du.Transpose_Times(f-Q);
    dE=dphi+dh*TV::Dot_Product(g,u)+h*dg.Transpose_Times(u)+h*du.Transpose_Times(g)+k*h*dh;
}
//#####################################################################
// Function Compute_ddE
//#####################################################################
template<class T,int d> void RC_EXTRAPOLATED<T,d>::HELPER::
Compute_ddE(const GENERAL_ENERGY<T>& base,T k,const TV& f,const int id)
{
    TV fm1=f-1;
    T m2=sqr(m),m3=m*m2;
    dds.From_Matrix(-xi*MATRIX<T,d>::Outer_Product((T)1/Q,ds)-MATRIX<T,d>::Outer_Product(s/Q,dxi)+DIAGONAL_MATRIX<T,d>((T)s*xi/sqr(Q))*dQ);
    for(int i=0;i<d;i++){
        MATRIX<T,d> t;
        for(int j=0;j<d;j++){t(i,j)+=ds(j);t(j,i)+=ds(j);}
        ddQ(i)=fm1(i)*dds+t.Symmetric_Part();}
    ddz+=SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dQ,DIAGONAL_MATRIX<T,d>((T)2*fm1/cube(Q)));
    ddz-=(DIAGONAL_MATRIX<T,d>((T)2/sqr(Q))*dQ).Symmetric_Part();
    for(int i=0;i<d;i++) ddz-=fm1(i)/sqr(Q(i))*ddQ(i);
    ddxi=2*cube(xi)*SYMMETRIC_MATRIX<T,d>::Outer_Product(dz)-sqr(xi)*ddz;
    base.dddE(Q,id,&TT(1));
    ddphi=SYMMETRIC_MATRIX<T,d>::Transpose_Times_With_Symmetric_Result(dg,dQ);
    for(int i=0;i<d;i++) ddphi+=g(i)*ddQ(i);
    ddm=3*m2*m3*SYMMETRIC_MATRIX<T,d>::Outer_Product(fm1)-m3;
    for(int i=0;i<d;i++){
        MATRIX<T,d> t;
        for(int j=0;j<d;j++){t(i,j)+=dm(j);t(j,i)+=dm(j);}
        ddu(i)=fm1(i)*ddm+t.Symmetric_Part();}
    for(int i=0;i<d;i++){
        ddg(i)=SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dQ,TT(i));
        for(int j=0;j<d;j++) ddg(i)+=H(i,j)*ddQ(j);}
    ddh=((T)1-dQ).Transpose_Times(du*2).Symmetric_Part();
    for(int i=0;i<d;i++) ddh+=(f(i)-Q(i))*ddu(i)-ddQ(i)*u(i);
    ddE=ddphi+TV::Dot_Product(g,u)*ddh+MATRIX<T,d>::Outer_Product(dh*2,dg.Transpose_Times(u)).Symmetric_Part();
    ddE+=MATRIX<T,d>::Outer_Product(dh*2,du.Transpose_Times(g)).Symmetric_Part();
    ddE+=2*h*dg.Transpose_Times(du).Symmetric_Part()+k*SYMMETRIC_MATRIX<T,d>::Outer_Product(dh)+k*h*ddh;
    for(int i=0;i<d;i++) ddE+=h*g(i)*ddu(i)+h*ddg(i)*u(i);
}
namespace PhysBAM{
template class RC_EXTRAPOLATED<float,2>;
template class RC_EXTRAPOLATED<float,3>;
template class RC_EXTRAPOLATED<double,2>;
template class RC_EXTRAPOLATED<double,3>;
}
