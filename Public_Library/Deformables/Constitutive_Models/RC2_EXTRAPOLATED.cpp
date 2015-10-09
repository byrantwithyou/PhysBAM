//#####################################################################
// Copyright 2011, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Tools/Polynomials/CUBIC.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/RC2_EXTRAPOLATED.h>
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
Update_Lame_Constants(const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient_input)
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
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T J=F.To_Vector().Product();
    if(J<extrapolation_cutoff){
        HELPER helper;
        bool b=helper.Compute_E(base,extrapolation_cutoff,F.To_Vector(),id);
        if(b) return helper.E;}
    return base.E(F.To_Vector(),id);
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> RC2_EXTRAPOLATED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    T J=F.To_Vector().Product();
    if(J<extrapolation_cutoff){
        HELPER helper;
        bool b=helper.Compute_E(base,extrapolation_cutoff,F.To_Vector(),id);
        if(b){
            helper.Compute_dE(base,F.To_Vector(),id);
//            Add_Debug_Particle(F.To_Vector(),VECTOR<T,3>(1,0,0));
//            Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,-helper.dE/youngs_modulus);
            return DIAGONAL_MATRIX<T,d>(helper.dE);}}
//    Add_Debug_Particle(F.To_Vector(),VECTOR<T,3>(0,1,0));
//    Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,-base.dE(F.To_Vector(),id)/youngs_modulus);

    return DIAGONAL_MATRIX<T,d>(base.dE(F.To_Vector(),id));
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const RC2_EXTRAPOLATED<T,2>& re,const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int id)
{
    T J=F.To_Vector().Product();
    T x=F.x.x,y=F.x.y,xpy=x+y;
    if(fabs(xpy)<re.panic_threshold) xpy=xpy<0?-re.panic_threshold:re.panic_threshold;
    if(J<re.extrapolation_cutoff){
        typename RC2_EXTRAPOLATED<T,2>::HELPER helper;
        bool b=helper.Compute_E(re.base,re.extrapolation_cutoff,F.To_Vector(),id);
        if(b){
            helper.Compute_dE(re.base,F.To_Vector(),id);
            helper.Compute_ddE(re.base,F.To_Vector(),id);
            dP_dF.x0000=helper.ddE.x00;
            dP_dF.x1100=helper.ddE.x10;
            dP_dF.x1111=helper.ddE.x11;
            T S=helper.dE.Sum()/xpy,D=helper.dE_it.x;
            dP_dF.x1001=(D-S)/2;
            dP_dF.x1010=(D+S)/2;
            return;}}
    dP_dF.x0000=re.base.Exx(x,y,id);
    dP_dF.x1100=re.base.Exy(x,y,id);
    dP_dF.x1111=re.base.Eyy(x,y,id);
    T S=re.P_From_Strain(F,id).Trace()/xpy,D=re.base.Ex_Ey_x_y(x,y,id);
    dP_dF.x1001=(D-S)/2;
    dP_dF.x1010=(D+S)/2;
}
///#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T> static void
Isotropic_Stress_Derivative_Helper(const RC2_EXTRAPOLATED<T,3>& re,const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int id)
{
    T J=F.To_Vector().Product();
    T x=F.x.x,y=F.x.y,z=F.x.z,xpy=x+y,xpz=x+z,ypz=y+z;
    if(fabs(xpy)<re.panic_threshold) xpy=xpy<0?-re.panic_threshold:re.panic_threshold;
    if(fabs(xpz)<re.panic_threshold) xpz=xpz<0?-re.panic_threshold:re.panic_threshold;
    if(fabs(ypz)<re.panic_threshold) ypz=ypz<0?-re.panic_threshold:re.panic_threshold;
    if(J<re.extrapolation_cutoff){
        typename RC2_EXTRAPOLATED<T,3>::HELPER helper;
        bool b=helper.Compute_E(re.base,re.extrapolation_cutoff,F.To_Vector(),id);
        if(b){
            helper.Compute_dE(re.base,F.To_Vector(),id);
            helper.Compute_ddE(re.base,F.To_Vector(),id);
            dP_dF.x0000=helper.ddE.x00;
            dP_dF.x1111=helper.ddE.x11;
            dP_dF.x2222=helper.ddE.x22;
            dP_dF.x1100=helper.ddE.x10;
            dP_dF.x2200=helper.ddE.x20;
            dP_dF.x2211=helper.ddE.x21;
            T ss1=sqr(F.x.x),ss2=sqr(F.x.y),ss3=sqr(F.x.z);
            T s01=ss1-ss2,s02=ss1-ss3,s12=ss2-ss3;
            if(fabs(s01)<re.panic_threshold) s01=s01<0?-re.panic_threshold:re.panic_threshold;
            if(fabs(s02)<re.panic_threshold) s02=s02<0?-re.panic_threshold:re.panic_threshold;
            if(fabs(s12)<re.panic_threshold) s12=s12<0?-re.panic_threshold:re.panic_threshold;
            T S01=(helper.dE.x+helper.dE.y)/xpy,D01=helper.dE_it.z;
            dP_dF.x1001=(D01-S01)/2;
            dP_dF.x1010=(D01+S01)/2;
            T S02=(helper.dE.x+helper.dE.z)/xpz,D02=helper.dE_it.y;
            dP_dF.x2002=(D02-S02)/2;
            dP_dF.x2020=(D02+S02)/2;
            T S12=(helper.dE.y+helper.dE.z)/ypz,D12=helper.dE_it.x;
            dP_dF.x2112=(D12-S12)/2;
            dP_dF.x2121=(D12+S12)/2;
            return;}}
    dP_dF.x0000=re.base.Exx(x,y,z,id);
    dP_dF.x1111=re.base.Eyy(x,y,z,id);
    dP_dF.x2222=re.base.Ezz(x,y,z,id);
    dP_dF.x1100=re.base.Exy(x,y,z,id);
    dP_dF.x2200=re.base.Exz(x,y,z,id);
    dP_dF.x2211=re.base.Eyz(x,y,z,id);
    VECTOR<T,3> P=re.P_From_Strain(F,id).To_Vector();
    T S01=(P.x+P.y)/xpy,D01=re.base.Ex_Ey_x_y(x,y,z,id);
    dP_dF.x1001=(D01-S01)/2;
    dP_dF.x1010=(D01+S01)/2;
    T S02=(P.x+P.z)/xpz,D02=re.base.Ex_Ez_x_z(x,y,z,id);
    dP_dF.x2002=(D02-S02)/2;
    dP_dF.x2020=(D02+S02)/2;
    T S12=(P.y+P.z)/ypz,D12=re.base.Ey_Ez_y_z(x,y,z,id);
    dP_dF.x2112=(D12-S12)/2;
    dP_dF.x2121=(D12+S12)/2;
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const
{
    return Isotropic_Stress_Derivative_Helper(*this,F,dP_dF,id);
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> RC2_EXTRAPOLATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=(alpha.m?alpha(id):constant_alpha),id_beta=(beta.m?beta(id):constant_beta);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
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
P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=(alpha.m?alpha(id):constant_alpha),id_beta=(beta.m?beta(id):constant_beta);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(id_alpha/TV::dimension+dd*dd)-dd;
    SYMMETRIC_MATRIX<T,d> s=sb*strain_rate+sa*strain_rate.Trace();
    *(MATRIX<T,d>*)aggregate.Get_Array_Pointer()+=s;
}
//#####################################################################
// Function P_From_Strain_Rate_Second_Half
//#####################################################################
template<class T,int d> MATRIX<T,d> RC2_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const int id) const
{
    T id_alpha=(alpha.m?alpha(id):constant_alpha),id_beta=(beta.m?beta(id):constant_beta);
    SYMMETRIC_MATRIX<T,d> strain_rate=(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(id_alpha/TV::dimension+dd*dd)-dd;
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
Compute_s(const VECTOR<T,3>& f,const int id,T extrapolation_cutoff)
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
Compute_E(const GENERAL_ENERGY<T>& base,T extrapolation_cutoff,const TV& f,const int id)
{
    TV fm1=f-1;
    s=Compute_s(f,id,extrapolation_cutoff);
    if(s==-1) return false;
    q=fm1*s+1;
    z=(fm1/q).Sum();
    xi=1/z;
    phi=base.E(q,id);
    m=1/fm1.Magnitude();
    u=m*fm1;
    g=base.dE(q,id);
    h=TV::Dot_Product(f-q,u);
    H=base.ddE(q,id);
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
Compute_dE(const GENERAL_ENERGY<T>& base,const TV& f,const int id)
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
    base.dddE(q,id,&TT(0));
    for(int i=0;i<d;i++) Tu+=u(i)*TT(i);
    uTu=Tu*u;
    dc=dq.Transpose_Times(uTu)+(T)2*du.Transpose_Times(Hu);
    dE=dphi+dh*(b+h*c)+h*db+(T)0.5*sqr(h)*dc;
}
//#####################################################################
// Function id_it
//#####################################################################
template<class T> MATRIX<T,3>
id_it(const VECTOR<T,3>& f)
{
    T a=1/(f.z-f.y),b=1/(f.x-f.z),c=1/(f.y-f.x);
    return MATRIX<T,3>(0,-a,a,b,0,-b,-c,c,0);
}
//#####################################################################
// Function id_it
//#####################################################################
template<class T> MATRIX<T,2,1>
id_it(const VECTOR<T,2>& f)
{
    T a=1/(f.x-f.y);
    MATRIX<T,2,1> r;
    r(1,1)=a;
    r(2,1)=-a;
    return r;
}
//#####################################################################
// Function id_it
//#####################################################################
template<class T> VECTOR<T,3>
qiqi(const VECTOR<T,3>& q)
{
    return (T)1/VECTOR<T,3>(q.z*q.y,q.x*q.z,q.x*q.y);
}
//#####################################################################
// Function id_it
//#####################################################################
template<class T> VECTOR<T,1>
qiqi(const VECTOR<T,2>& q)
{
    return VECTOR<T,1>(1/(q.x*q.y));
}
//#####################################################################
// Function Compute_ddE
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::HELPER::
Compute_ddE(const GENERAL_ENERGY<T>& base,const TV& f,const int id)
{
    TV fm1=f-1;
    T m2=sqr(m),m3=m*m2;
    dds.From_Matrix(-xi*MATRIX<T,d>::Outer_Product((T)1/q,ds)-MATRIX<T,d>::Outer_Product(s/q,dxi)+DIAGONAL_MATRIX<T,d>((T)s*xi/sqr(q))*dq);
    for(int i=0;i<d;i++){
        MATRIX<T,d> t;
        for(int j=0;j<d;j++){t(i,j)+=ds(j);t(j,i)+=ds(j);}
        ddq(i)=fm1(i)*dds+t.Symmetric_Part();}
    ddz+=SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dq,DIAGONAL_MATRIX<T,d>((T)2*fm1/cube(q)));
    ddz-=(DIAGONAL_MATRIX<T,d>((T)2/sqr(q))*dq).Symmetric_Part();
    for(int i=0;i<d;i++) ddz-=fm1(i)/sqr(q(i))*ddq(i);
    ddxi=2*cube(xi)*SYMMETRIC_MATRIX<T,d>::Outer_Product(dz)-sqr(xi)*ddz;
    ddphi=SYMMETRIC_MATRIX<T,d>::Transpose_Times_With_Symmetric_Result(dg,dq);
    for(int i=0;i<d;i++) ddphi+=g(i)*ddq(i);
    ddm=3*m2*m3*SYMMETRIC_MATRIX<T,d>::Outer_Product(fm1)-m3;
    for(int i=0;i<d;i++){
        MATRIX<T,d> t;
        for(int j=0;j<d;j++){t(i,j)+=dm(j);t(j,i)+=dm(j);}
        ddu(i)=fm1(i)*ddm+t.Symmetric_Part();}
    for(int i=0;i<d;i++){
        ddg(i)=SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dq,TT(i));
        for(int j=0;j<d;j++) ddg(i)+=H(i,j)*ddq(j);}
    ddh=((T)1-dq).Transpose_Times(du*2).Symmetric_Part();
    for(int i=0;i<d;i++) ddh+=(f(i)-q(i))*ddu(i)-ddq(i)*u(i);
    ddb=(T)2*dg.Transpose_Times(du).Symmetric_Part();
    for(int i=0;i<d;i++) ddb+=g(i)*ddu(i)+ddg(i)*u(i);
    base.ddddE(q,id,&A(0)(0));
    ddc=(T)4*dq.Transpose_Times(Tu*du).Symmetric_Part()+(T)2*SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(du,H);
    for(int i=0;i<d;i++) for(int j=0;j<d;j++) ddc+=u(i)*u(j)*SYMMETRIC_MATRIX<T,d>::Conjugate_With_Transpose(dq,A(i)(j));
    for(int i=0;i<d;i++) ddc+=uTu(i)*ddq(i)+(T)2*Hu(i)*ddu(i);
    ddE=ddphi+(b+c*h)*ddh+h*ddb+(T).5*sqr(h)*ddc+MATRIX<T,d>::Outer_Product(dh*2,db+h*dc).Symmetric_Part()+c*SYMMETRIC_MATRIX<T,d>::Outer_Product(dh);

    base.Compute_it(q,id,g_it,H_xit,H_iitt,T_xxit,T_xiitt,T_iiittt,T_itit);
    g_it*=s; // Correct for the location of f differing from that of q.
    H_xit*=s;
    H_iitt*=s;
    T_xxit*=s;
    T_xiitt*=s;
    T_iiittt*=s;
    T_itit*=s;
    u_it.Fill(m);
    VECTOR<T,TV::SPIN::m> Hu_it,uTu_it;
    for(int i=0;i<TV::SPIN::m;i++){
        int a=(TV::m==3?(i+1)%3:i),b=(a+1)%3;
        Hu_it(i)=H(a,a)*u_it(i)+H_iitt(i)*u(b)-H(a,b)*u_it(i);
        T u2_it=sqr(m)*(fm1(a)+fm1(b));
        uTu_it(i)=2*T_itit(i)*u(a)*u(b)+T_iiittt(i)*sqr(u(a))+TT(b)(b,b)*u2_it-T_itit(i)*sqr(u(a))-TT(a)(b,b)*u2_it;}
    for(int i=0;i<(TV::m==3)*TV::m;i++){
        int a=(i+1)%3,b=(a+1)%3;
        Hu_it(i)+=H_xit(i)*u(i);
        uTu_it(i)+=T_xxit(i)*sqr(u(i))+2*u(i)*T_xiitt(i)*u(b)+2*u(i)*(TT(i)(a,a)*u_it(i)-TT(i)(a,b)*u_it(i));}

    dE_it=g_it+h*Hu_it+(T).5*sqr(s*h)/m*TV::Dot_Product(uTu,u)*xi*SYMMETRIC_MATRIX<T,d>::Outer_Product((T)1/q).Off_Diagonal_Part()+(T).5*sqr(h)*s*uTu_it;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T>
static void Print_Helper(const char* str,T a0,T a1)
{
    char buff[1000];
    sprintf(buff,"============ test ============ %s %8.5f %8.5f (%8.5f)\n",str,a0,a1,fabs(a0-a1));
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV>
static void Test_Model_Helper(const char* str,T a0,T a1,TV da0,TV da1,TV df,T e)
{
    T av=TV::Dot_Product(da1+da0,df)/2/e;
    T dif=(a1-a0)/e;
    Print_Helper(str,av,dif);
}
//#####################################################################
// Function Print_Helper
//#####################################################################
template<class T,int d>
static void Print_Helper(const char* str,const VECTOR<T,d>& a0,const VECTOR<T,d> &a1)
{
    char buff[1000];
    sprintf(buff,"============ test ============ %s %8.5f %8.5f (%8.5f)\n",str,a0.Magnitude(),a1.Magnitude(),(a0-a1).Magnitude());
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV,int d>
static void Test_Model_Helper(const char* str,TV a0,TV a1,const MATRIX<T,d,d>& da0,const MATRIX<T,d,d>& da1,TV df,T e)
{
    TV av=(da1+da0)*df/2/e;
    TV dif=(a1-a0)/e;
    Print_Helper(str,av,dif);
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV,int d>
static void Test_Model_Helper(const char* str,TV a0,TV a1,const SYMMETRIC_MATRIX<T,d>& da0,const SYMMETRIC_MATRIX<T,d>& da1,TV df,T e)
{
    Test_Model_Helper(str,a0,a1,MATRIX<T,d>(da0),MATRIX<T,d>(da1),df,e);
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class T,class TV,int d>
static void Test_Model_Helper(const char* str,const MATRIX<T,d,d>& a0,const MATRIX<T,d,d>& a1,const VECTOR<SYMMETRIC_MATRIX<T,d>,d>& da0,const VECTOR<SYMMETRIC_MATRIX<T,d>,d>& da1,TV df,T e)
{
    for(int i=0;i<TV::m;i++){
        TV av=(da1(i)+da0(i))*df/2/e;
        TV dif=(a1.Row(i)-a0.Row(i))/e;
        char buff[1000];
        sprintf(buff,"============ test ============ %s %8.5f %8.5f (%8.5f)\n",str,av.Magnitude(),dif.Magnitude(),(av-dif).Magnitude());
        LOG::cout<<buff;}
}
//#####################################################################
// Function Test_it
//#####################################################################
template<class T> static void
Test_it(const char* str,const VECTOR<T,3>& f,const VECTOR<T,3>& w,const VECTOR<T,3>& w_it)
{
    Print_Helper(str,(w.y-w.z)/(f.y-f.z),w_it.x);
    Print_Helper(str,(w.x-w.z)/(f.x-f.z),w_it.y);
    Print_Helper(str,(w.x-w.y)/(f.x-f.y),w_it.z);
}
//#####################################################################
// Function Test_it
//#####################################################################
template<class T> static void
Test_it(const char* str,const VECTOR<T,2>& f,const VECTOR<T,2>& w,const VECTOR<T,1>& w_it)
{
    Print_Helper(str,(w.x-w.y)/(f.x-f.y),w_it.x);
}
//#####################################################################
// Function Test_it
//#####################################################################
template<class T> static void
Test_it(const char* str,const VECTOR<T,3>& f,const MATRIX<T,3>& w,const MATRIX<T,3>& w_it)
{
    Print_Helper(str,(w.Column(1)-w.Column(2))/(f.y-f.z),w_it.Column(0));
    Print_Helper(str,(w.Column(0)-w.Column(2))/(f.x-f.z),w_it.Column(1));
    Print_Helper(str,(w.Column(0)-w.Column(1))/(f.x-f.y),w_it.Column(2));
}
//#####################################################################
// Function Test_it
//#####################################################################
template<class T> static void
Test_it(const char* str,const VECTOR<T,2>& f,const MATRIX<T,2>& w,MATRIX<T,2,1> w_it)
{
    Print_Helper(str,(w.Column(0)-w.Column(1))/(f.x-f.y),w_it.Column(0));
}
//#####################################################################
// Function Test_Model
//#####################################################################
template<class T,int d> void RC2_EXTRAPOLATED<T,d>::
Test_Model() const
{
    RANDOM_NUMBERS<T> rand;
    for(int i=0;i<20;i++){
        TV f;
        rand.Fill_Uniform(f,0,2);
        T e=1e-5;
        TV df;
        rand.Fill_Uniform(df,-e,e);
        f=f.Sorted().Reversed();
        if(rand.Get_Uniform_Integer(0,1)==1) f(2)=-f(2);
        LOG::cout<<f<<std::endl;
        this->Test(DIAGONAL_MATRIX<T,TV::m>(f),1);

        int id=0;
        if(f.Product()>extrapolation_cutoff) continue;
        HELPER h0;
        if(!h0.Compute_E(base,extrapolation_cutoff,f,id)) continue;
        h0.Compute_dE(base,f,id);
        h0.Compute_ddE(base,f,id);
        HELPER h1;
        if(!h1.Compute_E(base,extrapolation_cutoff,f+df,id)) continue;
        h1.Compute_dE(base,f+df,id);
        h1.Compute_ddE(base,f+df,id);
#define XX(k) Test_Model_Helper(#k,h0.k,h1.k,h0.d##k,h1.d##k,df,e);Test_Model_Helper(#k,h0.d##k,h1.d##k,h0.dd##k,h1.dd##k,df,e);
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

#define IT(w) Test_it(#w "_it",f,h0.w,h0.w##_it)
        IT(dE);
        IT(g);
        IT(u);
    }
}
namespace PhysBAM{
template class RC2_EXTRAPOLATED<float,2>;
template class RC2_EXTRAPOLATED<float,3>;
template class RC2_EXTRAPOLATED<double,2>;
template class RC2_EXTRAPOLATED<double,3>;
}
