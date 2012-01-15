//#####################################################################
// Copyright 2011, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Tools/Polynomials/CUBIC.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/RC2_EXTRAPOLATED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> RC2_EXTRAPOLATED<T,d>::
RC2_EXTRAPOLATED(GENERAL_ENERGY<T>& ge_input,const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input,
    const T extra_force_coefficient_input)
    :base(ge_input),youngs_modulus(youngs_modulus_input),poissons_ratio(poissons_ratio_input),extrapolation_cutoff(extrapolation_cutoff_input),panic_threshold((T)1e-6)
{
    Update_Lame_Constants(youngs_modulus_input,poissons_ratio_input,Rayleigh_coefficient);
    base.Initialize(constant_mu,constant_lambda);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> RC2_EXTRAPOLATED<T,d>::
~RC2_EXTRAPOLATED()
{
    delete &base;
}
//#####################################################################
// Update Lame Constants
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::
Update_Lame_Constants(const T youngs_modulus_input, const T poissons_ratio_input,const T Rayleigh_coefficient_input)
{
    assert(poissons_ratio>-1&&poissons_ratio<.5);
    constant_lambda=youngs_modulus_input*poissons_ratio_input/((1+poissons_ratio_input)*(1-2*poissons_ratio_input));
    constant_mu=youngs_modulus_input/(2*(1+poissons_ratio_input));
    constant_alpha=Rayleigh_coefficient_input*constant_lambda;
    constant_beta=Rayleigh_coefficient_input*constant_mu;
    base.Initialize(constant_mu,constant_lambda);
    youngs_modulus=youngs_modulus_input; poissons_ratio=poissons_ratio_input;
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T RC2_EXTRAPOLATED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{

    T J=F.To_Vector().Product();
    if(J<extrapolation_cutoff){
        HELPER helper;
        bool b=helper.Compute_E(base,extrapolation_cutoff,F.To_Vector(),simplex);
        if(b) return helper.E;}
    return base.E(F.To_Vector(),simplex);
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> RC2_EXTRAPOLATED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const T scale,const int simplex) const
{
    T J=F.To_Vector().Product();
    if(J<extrapolation_cutoff){
        HELPER helper;
        bool b=helper.Compute_E(base,extrapolation_cutoff,F.To_Vector(),simplex);
        if(b){
            helper.Compute_dE(base,F.To_Vector(),simplex);
//            Add_Debug_Particle(F.To_Vector(),VECTOR<T,3>(1,0,0));
//            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,-helper.dE/youngs_modulus);
            return scale*DIAGONAL_MATRIX<T,d>(helper.dE);}}
//    Add_Debug_Particle(F.To_Vector(),VECTOR<T,3>(0,1,0));
//    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,-base.dE(F.To_Vector(),simplex)/youngs_modulus);

    return scale*DIAGONAL_MATRIX<T,d>(base.dE(F.To_Vector(),simplex));
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const RC2_EXTRAPOLATED<T,2>& re,const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int simplex)
{
    T J=F.To_Vector().Product();
    if(J<re.extrapolation_cutoff){
        typename RC2_EXTRAPOLATED<T,2>::HELPER helper;
        bool b=helper.Compute_E(re.base,re.extrapolation_cutoff,F.To_Vector(),simplex);
        if(b){
            helper.Compute_dE(re.base,F.To_Vector(),simplex);
            helper.Compute_ddE(re.base,F.To_Vector(),simplex);
            dP_dF.x1111=helper.ddE.x11;
            dP_dF.x2211=helper.ddE.x21;
            dP_dF.x2222=helper.ddE.x22;
            T ss1=sqr(F.x11),ss2=sqr(F.x22);
            T s12=ss1-ss2;
            if(fabs(s12)<re.panic_threshold) s12=s12<0?-re.panic_threshold:re.panic_threshold;
            dP_dF.x2112=(-helper.dE.y*F.x11+helper.dE.x*F.x22)/s12;
            dP_dF.x2121=(-helper.dE.y*F.x22+helper.dE.x*F.x11)/s12;
            return;}}
    T x = F.x11, y = F.x22, xpy = x+y;
    dP_dF.x1111 = re.base.Exx(x,y,simplex);
    dP_dF.x2211 = re.base.Exy(x,y,simplex);
    dP_dF.x2222 = re.base.Eyy(x,y,simplex);
    if(fabs(xpy)<re.panic_threshold) xpy=xpy<0?-re.panic_threshold:re.panic_threshold;
    T S=re.P_From_Strain(F,1,simplex).Trace()/xpy, D=re.base.Ex_Ey_x_y(x,y,simplex);
    dP_dF.x2112 = (D-S)/2;
    dP_dF.x2121 = (D+S)/2;
}
///#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const RC2_EXTRAPOLATED<T,3>& re,const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int simplex)
{
    T J=F.To_Vector().Product();
    if(J<re.extrapolation_cutoff){
        typename RC2_EXTRAPOLATED<T,3>::HELPER helper;
        bool b=helper.Compute_E(re.base,re.extrapolation_cutoff,F.To_Vector(),simplex);
        if(b){
            helper.Compute_dE(re.base,F.To_Vector(),simplex);
            helper.Compute_ddE(re.base,F.To_Vector(),simplex);
            dP_dF.x1111=helper.ddE.x11;
            dP_dF.x2222=helper.ddE.x22;
            dP_dF.x3333=helper.ddE.x33;
            dP_dF.x2211=helper.ddE.x21;
            dP_dF.x3311=helper.ddE.x31;
            dP_dF.x3322=helper.ddE.x32;
            T ss1=sqr(F.x11),ss2=sqr(F.x22),ss3=sqr(F.x33);
            T s12=ss1-ss2,s13=ss1-ss3,s23=ss2-ss3;
            if(fabs(s12)<re.panic_threshold) s12=s12<0?-re.panic_threshold:re.panic_threshold;
            if(fabs(s13)<re.panic_threshold) s13=s13<0?-re.panic_threshold:re.panic_threshold;
            if(fabs(s23)<re.panic_threshold) s23=s23<0?-re.panic_threshold:re.panic_threshold;
            dP_dF.x2112=(-helper.dE.y*F.x11+helper.dE.x*F.x22)/s12;
            dP_dF.x2121=(-helper.dE.y*F.x22+helper.dE.x*F.x11)/s12;
            dP_dF.x3113=(-helper.dE.z*F.x11+helper.dE.x*F.x33)/s13;
            dP_dF.x3131=(-helper.dE.z*F.x33+helper.dE.x*F.x11)/s13;
            dP_dF.x3223=(-helper.dE.z*F.x22+helper.dE.y*F.x33)/s23;
            dP_dF.x3232=(-helper.dE.z*F.x33+helper.dE.y*F.x22)/s23;
            return;}}
    T x = F.x11, y = F.x22, z = F.x33, xpy = x+y, xpz = x+z, ypz = y+z;
    dP_dF.x1111 = re.base.Exx(x,y,z,simplex);
    dP_dF.x2222 = re.base.Eyy(x,y,z,simplex);
    dP_dF.x3333 = re.base.Ezz(x,y,z,simplex);
    dP_dF.x2211 = re.base.Exy(x,y,z,simplex);
    dP_dF.x3311 = re.base.Exz(x,y,z,simplex);
    dP_dF.x3322 = re.base.Eyz(x,y,z,simplex);
    if(fabs(xpy)<re.panic_threshold) xpy=xpy<0?-re.panic_threshold:re.panic_threshold;
    VECTOR<T,3> P=re.P_From_Strain(F,1,simplex).To_Vector();
    T S12=(P.x+P.y)/xpy, D12=re.base.Ex_Ey_x_y(x,y,z,simplex);
    dP_dF.x2112 = (D12-S12)/2;
    dP_dF.x2121 = (D12+S12)/2;
    if(fabs(xpz)<re.panic_threshold) xpz=xpz<0?-re.panic_threshold:re.panic_threshold;
    T S13=(P.x+P.z)/xpz, D13=re.base.Ex_Ez_x_z(x,y,z,simplex);
    dP_dF.x3113 = (D13-S13)/2;
    dP_dF.x3131 = (D13+S13)/2;
    if(fabs(ypz)<re.panic_threshold) ypz=ypz<0?-re.panic_threshold:re.panic_threshold;
    T S23=(P.y+P.z)/ypz, D23=re.base.Ey_Ez_y_z(x,y,z,simplex);
    dP_dF.x3223 = (D23-S23)/2;
    dP_dF.x3232 = (D23+S23)/2;
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int triangle) const
{
    return Isotropic_Stress_Derivative_Helper(*this,F,dP_dF,triangle);
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> RC2_EXTRAPOLATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int RC2_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> RC2_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const
{
    SYMMETRIC_MATRIX<T,d> strain_rate=scale*(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*constant_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(constant_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
//#####################################################################
// Function Compute_s
//#####################################################################
template<class T> static T
Compute_s(const VECTOR<T,2>& f,const int simplex,T extrapolation_cutoff)
{
    VECTOR<T,2> fm1=f-1;
    QUADRATIC<T> quadratic(fm1.Product(),fm1.Sum(),1-extrapolation_cutoff);
    T a=quadratic(0);
    T b=quadratic(1);
    PHYSBAM_ASSERT(a>0);
    if(b==0) return 1;
    if((a>0)==(b>0)) return -1;
    ITERATIVE_SOLVER<T> iterative_solver;
    return iterative_solver.Bisection_Secant_Root(quadratic,0,1);
}
//#####################################################################
// Function Compute_s
//#####################################################################
template<class T> static T
Compute_s(const VECTOR<T,3>& f,const int simplex,T extrapolation_cutoff)
{
    VECTOR<T,3> fm1=f-1;
    CUBIC<T> cubic(fm1.Product(),DIAGONAL_MATRIX<T,3>(fm1).Cofactor_Matrix().Trace(),fm1.Sum(),1-extrapolation_cutoff);
    T a=cubic(0);
    T b=cubic(1);
    PHYSBAM_ASSERT(a>0);
    if(b==0) return 1;
    if((a>0)==(b>0)) return -1;
    ITERATIVE_SOLVER<T> iterative_solver;
    return iterative_solver.Bisection_Secant_Root(cubic,0,1);
}
//#####################################################################
// Function Compute_E
//#####################################################################
template<class T,int d> bool RC2_EXTRAPOLATED<T,d>::HELPER::
Compute_E(const GENERAL_ENERGY<T>& base,T extrapolation_cutoff,const TV& f,const int simplex)
{
    TV fm1=f-1;
    s=Compute_s(f,simplex,extrapolation_cutoff);
    if(s==-1) return false;
    q=fm1*s+1;
    z=(fm1/q).Sum();
    xi=1/z;
    phi=base.E(q,simplex);
    m=1/fm1.Magnitude();
    u=m*fm1;
    g=base.dE(q,simplex);
    h=TV::Dot_Product(f-q,u);
    H=base.ddE(q,simplex);
    Hu=H*u;
    b=TV::Dot_Product(g,u);
    c=TV::Dot_Product(u,Hu);
    E=phi+h*b+(T)0.5*sqr(h)*c;
    return true;
}
//#####################################################################
// Function Compute_dE
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::HELPER::
Compute_dE(const GENERAL_ENERGY<T>& base,const TV& f,const int simplex)
{
    TV fm1=f-1;
    ds=-s*xi/q;
    dq=s+MATRIX<T,d>::Outer_Product(fm1,ds);
    dz=(T)1/q-dq.Transpose_Times(fm1/sqr(q));
    dxi=-sqr(xi)*dz;
    dphi=dq.Transpose_Times(g);
    dm=-cube(m)*fm1;
    du=MATRIX<T,d>::Outer_Product(fm1,dm)+m;
    dg=H*dq;
    dh=((T)1-dq).Transpose_Times(u)+du.Transpose_Times(f-q);
    db=dg.Transpose_Times(u)+du.Transpose_Times(g);
    base.dddE(q,simplex,&TT(1));
    for(int i=1; i<=d; i++) Tu+=u(i)*TT(i);
    uTu=Tu*u;
    dc=dq.Transpose_Times(uTu)+(T)2*du.Transpose_Times(Hu);
    dE=dphi+dh*(b+h*c)+h*db+(T)0.5*sqr(h)*dc;
}
//#####################################################################
// Function Compute_ddE
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::HELPER::
Compute_ddE(const GENERAL_ENERGY<T>& base,const TV& f,const int simplex)
{
    TV fm1=f-1;
    T m2=sqr(m),m3=m*m2;
    dds.From_Matrix(-xi*MATRIX<T,d>::Outer_Product((T)1/q,ds)-MATRIX<T,d>::Outer_Product(s/q,dxi)+DIAGONAL_MATRIX<T,d>((T)s*xi/sqr(q))*dq);
    for(int i=1; i<=d; i++){
        MATRIX<T,d> t;
        for(int j=1; j<=d; j++){t(i,j)+=ds(j);t(j,i)+=ds(j);}
        ddq(i)=fm1(i)*dds+t.Symmetric_Part();}
    ddz+=SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dq,DIAGONAL_MATRIX<T,d>((T)2*fm1/cube(q)));
    ddz-=(DIAGONAL_MATRIX<T,d>((T)2/sqr(q))*dq).Symmetric_Part();
    for(int i=1; i<=d; i++) ddz-=fm1(i)/sqr(q(i))*ddq(i);
    ddxi=2*cube(xi)*SYMMETRIC_MATRIX<T,d>::Outer_Product(dz)-sqr(xi)*ddz;
    ddphi=SYMMETRIC_MATRIX<T,d>::Transpose_Times_With_Symmetric_Result(dg,dq);
    for(int i=1;i<=d;i++) ddphi+=g(i)*ddq(i);
    ddm=3*m2*m3*SYMMETRIC_MATRIX<T,d>::Outer_Product(fm1)-m3;
    for(int i=1; i<=d; i++){
        MATRIX<T,d> t;
        for(int j=1; j<=d; j++){t(i,j)+=dm(j);t(j,i)+=dm(j);}
        ddu(i)=fm1(i)*ddm+t.Symmetric_Part();}
    for(int i=1; i<=d; i++){
        ddg(i)=SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dq,TT(i));
        for(int j=1; j<=d; j++) ddg(i)+=H(i,j)*ddq(j);}
    ddh=((T)1-dq).Transpose_Times(du*2).Symmetric_Part();
    for(int i=1; i<=d; i++) ddh+=(f(i)-q(i))*ddu(i)-ddq(i)*u(i);
    ddb=(T)2*dg.Transpose_Times(du).Symmetric_Part();
    for(int i=1; i<=d; i++) ddb+=g(i)*ddu(i)+ddg(i)*u(i);
    base.ddddE(q,simplex,&A(1)(1));
    ddc=(T)4*dq.Transpose_Times(Tu*du).Symmetric_Part()+(T)2*SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(du,H);
    for(int i=1; i<=d; i++) for(int j=1; j<=d; j++) ddc+=u(i)*u(j)*SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dq,A(i)(j));
    for(int i=1; i<=d; i++) ddc+=uTu(i)*ddq(i)+(T)2*Hu(i)*ddu(i);
    ddE=ddphi+(b+c*h)*ddh+h*ddb+(T).5*sqr(h)*ddc+MATRIX<T,d>::Outer_Product(dh*2,db+h*dc).Symmetric_Part()+c*SYMMETRIC_MATRIX<T,d>::Outer_Product(dh);
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV>
static void Test_Model_Helper(const char* str,T a0, T a1, TV da0, TV da1, TV df, T e)
{
    T av=TV::Dot_Product(da1+da0,df)/2/e;
    T dif=(a1-a0)/e;
    char buff[1000];
    sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av, dif, fabs(av-dif));
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV,int d>
static void Test_Model_Helper(const char* str,TV a0, TV a1, const MATRIX<T,d>& da0, const MATRIX<T,d>& da1, TV df, T e)
{
    TV av=(da1+da0)*df/2/e;
    TV dif=(a1-a0)/e;
    char buff[1000];
    sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV,int d>
static void Test_Model_Helper(const char* str,TV a0, TV a1, const SYMMETRIC_MATRIX<T,d>& da0, const SYMMETRIC_MATRIX<T,d>& da1, TV df, T e)
{
    Test_Model_Helper(str,a0, a1, MATRIX<T,d>(da0), MATRIX<T,d>(da1), df, e);
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV,int d>
static void Test_Model_Helper(const char* str,const MATRIX<T,d>& a0, const MATRIX<T,d>& a1, const VECTOR<SYMMETRIC_MATRIX<T,d>,d>& da0, const VECTOR<SYMMETRIC_MATRIX<T,d>,d>& da1, TV df, T e)
{
    for(int i=1;i<=TV::m;i++){
        TV av=(da1(i)+da0(i))*df/2/e;
        TV dif=(a1.Transposed().Column(i)-a0.Transposed().Column(i))/e;
        char buff[1000];
        sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());
        LOG::cout<<buff;}
}
//#####################################################################
// Function Test_Model
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::
Test_Model() const
{
    RANDOM_NUMBERS<T> rand;
    for(int i=1;i<=20;i++){
        TV f;
        rand.Fill_Uniform(f,0,2);
        T e=1e-5;
        TV df;
        rand.Fill_Uniform(df,-e,e);
        f=f.Sorted().Reversed();
        if(rand.Get_Uniform_Integer(0,1)==1) f(2)=-f(2);
        LOG::cout<<f<<std::endl;
//        Test(DIAGONAL_MATRIX<T,3>(f),1);

        int simplex=0;
        if(f.Product()>extrapolation_cutoff) continue;
        HELPER h0;
        if(!h0.Compute_E(base,extrapolation_cutoff,f,simplex)) continue;
        h0.Compute_dE(base,f,simplex);
        h0.Compute_ddE(base,f,simplex);
        HELPER h1;
        if(!h1.Compute_E(base,extrapolation_cutoff,f+df,simplex)) continue;
        h1.Compute_dE(base,f+df,simplex);
        h1.Compute_ddE(base,f+df,simplex);
#define XX(k) Test_Model_Helper(#k,h0.k, h1.k, h0.d##k, h1.d##k, df, e);Test_Model_Helper(#k,h0.d##k, h1.d##k, h0.dd##k, h1.dd##k, df, e);
        XX(m);
        XX(h);
        XX(phi);
        XX(E);
        XX(z);
        XX(xi);
        XX(s);
        XX(q);
        XX(u);
        XX(g);
        XX(b);
        XX(c);
    }
}
template class RC2_EXTRAPOLATED<float,2>;
template class RC2_EXTRAPOLATED<float,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RC2_EXTRAPOLATED<double,2>;
template class RC2_EXTRAPOLATED<double,3>;
#endif
