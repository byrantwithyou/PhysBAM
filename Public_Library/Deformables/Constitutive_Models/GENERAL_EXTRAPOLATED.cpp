 //#####################################################################
// Copyright 2003-2011, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Alexey Stomakhin, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Math_Tools/cube.h>
#include <Tools/Math_Tools/pow.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/GENERAL_EXTRAPOLATED.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> GENERAL_EXTRAPOLATED<T,d>::
GENERAL_EXTRAPOLATED(GENERAL_ENERGY<T>& ge_input,const T youngs_modulus_input,const T poissons_ratio_input,const T Rayleigh_coefficient,const T extrapolation_cutoff_input,
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
template<class T,int d> GENERAL_EXTRAPOLATED<T,d>::
~GENERAL_EXTRAPOLATED()
{
    delete &base;
}
//#####################################################################
// Update Lame Constants
//#####################################################################
template<class T,int d> void GENERAL_EXTRAPOLATED<T,d>::
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
template<class T,int d> T GENERAL_EXTRAPOLATED<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    return Energy_Density_Helper(F,id);
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T GENERAL_EXTRAPOLATED<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,2>& F,const int id) const
{
    T a = extrapolation_cutoff;
    T x = std::max(F.x.x,a);
    T y = std::max(F.x.y,a);
    T dx = F.x.x - extrapolation_cutoff;
    T dy = F.x.y - extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
    T E=base.E(x,y,id);
    if(dx < 0) E+=base.Ex(x,y,id)*dx + k*dx*dx;
    if(dy < 0) E+=base.Ey(x,y,id)*dy + k*dy*dy;
    if((dx < 0) && (dy < 0)) E+=base.Exy(x,y,id)*dx*dy;
    return E;
}
//#####################################################################
// Function Energy_Density_Helper
//#####################################################################
template<class T,int d> T GENERAL_EXTRAPOLATED<T,d>::
Energy_Density_Helper(const DIAGONAL_MATRIX<T,3>& F,const int id) const
{
    T a = extrapolation_cutoff;
    T x = std::max(F.x.x,a);
    T y = std::max(F.x.y,a);
    T z = std::max(F.x.z,a);
    T dx = F.x.x - extrapolation_cutoff;
    T dy = F.x.y - extrapolation_cutoff;
    T dz = F.x.z - extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;
    T E=base.E(x,y,z,id);
    if(dx < 0) E+=base.Ex(x,y,z,id)*dx + k*dx*dx;
    if(dy < 0) E+=base.Ey(x,y,z,id)*dy + k*dy*dy;
    if(dz < 0) E+=base.Ez(x,y,z,id)*dz + k*dz*dz;
    if((dx < 0) && (dy < 0)) E+=base.Exy(x,y,z,id)*dx*dy;
    if((dx < 0) && (dz < 0)) E+=base.Exz(x,y,z,id)*dx*dz;
    if((dy < 0) && (dz < 0)) E+=base.Eyz(x,y,z,id)*dy*dz;
    if((dx < 0) && (dy < 0) && (dz < 0)) E+=base.Exyz(x,y,z,id)*dx*dy*dz;
    return E;
}
//#####################################################################
// Function P_From_Strain
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,d> GENERAL_EXTRAPOLATED<T,d>::
P_From_Strain(const DIAGONAL_MATRIX<T,d>& F,const int id) const
{
    return P_From_Strain_Helper(F,id);
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,2> GENERAL_EXTRAPOLATED<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,2>& F,const int id) const
{
    T a = extrapolation_cutoff;
    T x = std::max(F.x.x,a);
    T y = std::max(F.x.y,a);
    T dx = F.x.x - extrapolation_cutoff;
    T dy = F.x.y - extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;

    DIAGONAL_MATRIX<T,2> result;
    result.x.x = base.Ex(x,y,id);
    result.x.y = base.Ey(x,y,id);
    if(dx < 0)
    {
        result.x.x += 2*k*dx;
        result.x.y += base.Exy(x,y,id)*dx;
    }
    if(dy < 0)
    {
        result.x.x += base.Exy(x,y,id)*dy;
        result.x.y += 2*k*dy;
    }
    return result;
}
//#####################################################################
// Function P_From_Strain_Helper
//#####################################################################
template<class T,int d> DIAGONAL_MATRIX<T,3> GENERAL_EXTRAPOLATED<T,d>::
P_From_Strain_Helper(const DIAGONAL_MATRIX<T,3>& F,const int id) const
{
    T a = extrapolation_cutoff;
    T x = std::max(F.x.x,a);
    T y = std::max(F.x.y,a);
    T z = std::max(F.x.z,a);
    T dx = F.x.x - extrapolation_cutoff;
    T dy = F.x.y - extrapolation_cutoff;
    T dz = F.x.z - extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;

    DIAGONAL_MATRIX<T,3> result;
    result.x.x = base.Ex(x,y,z,id);
    result.x.y = base.Ey(x,y,z,id);
    result.x.z = base.Ez(x,y,z,id);
    if(dx < 0)
    {
        result.x.x += 2*k*dx;
        result.x.y += base.Exy(x,y,z,id)*dx;
        result.x.z += base.Exz(x,y,z,id)*dx;
    }
    if(dy < 0)
    {
        result.x.x += base.Exy(x,y,z,id)*dy;
        result.x.y += 2*k*dy;
        result.x.z += base.Eyz(x,y,z,id)*dy;
    }
    if(dz < 0)
    {
        result.x.x += base.Exz(x,y,z,id)*dz;
        result.x.y += base.Eyz(x,y,z,id)*dz;
        result.x.z += 2*k*dz;
    }
    if((dy < 0) && (dz < 0)) result.x.x += base.Exyz(x,y,z,id)*dy*dz;
    if((dx < 0) && (dz < 0)) result.x.y += base.Exyz(x,y,z,id)*dx*dz;
    if((dx < 0) && (dy < 0)) result.x.z += base.Exyz(x,y,z,id)*dx*dy;
    return result;
}
//#####################################################################
// Function Isotropic_Stress_Derivative
//#####################################################################
template<class T,int d> void GENERAL_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dP_dF,const int id) const
{
    return Isotropic_Stress_Derivative_Helper(F,dP_dF,id);
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void GENERAL_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,2>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dP_dF,const int id) const
{
    T a = extrapolation_cutoff;
    T x = std::max(F.x.x,a);
    T y = std::max(F.x.y,a);
    T dx = F.x.x - extrapolation_cutoff;
    T dy = F.x.y - extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;

    dP_dF.x0000 = 0;
    dP_dF.x1111 = 0;
    dP_dF.x1100 = base.Exy(x,y,id);

    if(dx < 0) dP_dF.x0000 += 2*k;
    else dP_dF.x0000 += base.Exx(x,y,id);

    if(dy < 0) dP_dF.x1111 += 2*k;
    else dP_dF.x1111 += base.Eyy(x,y,id);

    T xpy = F.x.x+F.x.y; if(fabs(xpy)<panic_threshold) xpy=xpy<0?-panic_threshold:panic_threshold;
    T r=((dx < 0) != (dy < 0))?((F.x.x!=F.x.y)?dx/(F.x.x-F.x.y):1):0;
    if(dy<0) r=1-r;

    if((dx < 0) && (dy >= 0)) // Rx
    {
        dP_dF.x1111 += base.Exyy(x,y,id)*dx;
    }
    else if((dx >= 0) && (dy < 0)) // Ry
    {
        dP_dF.x0000 += base.Exxy(x,y,id)*dy;
    }

    T S=P_From_Strain_Helper(F,id).Trace()/xpy, D=(2*k-base.Exy(x,y,id))*r+(1-r)*base.Ex_Ey_x_y(x,y,id);
    dP_dF.x1001 = (D-S)/2;
    dP_dF.x1010 = (D+S)/2;
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function Isotropic_Stress_Derivative_Helper_Part
//#####################################################################
template<class T,int d> void GENERAL_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative_Helper_Part(T fx,T fy,T fz,const int id,T& xxxx,T& yyzz,T& yzyz,T& yzzy) const
{
    T a = extrapolation_cutoff;
    T x = std::max(fx,a);
    T y = std::max(fy,a);
    T z = std::max(fz,a);
    T dx = fx - extrapolation_cutoff;
    T dy = fy - extrapolation_cutoff;
    T dz = fz - extrapolation_cutoff;
    T k = extra_force_coefficient*youngs_modulus;

    yyzz = base.Eyz(x,y,z,id);

    if(dx < 0)
    {
        xxxx = 2*k;
        yyzz += base.Exyz(x,y,z,id)*dx;
    }
    else
    {
        xxxx = base.Exx(x,y,z,id);
        if(dy < 0) xxxx += base.Exxy(x,y,z,id)*dy;
        if(dz < 0) xxxx += base.Exxz(x,y,z,id)*dz;
        if((dy < 0) && (dz < 0)) xxxx += base.Exxyz(x,y,z,id)*dy*dz;
    }

    DIAGONAL_MATRIX<T,3> P=P_From_Strain_Helper(DIAGONAL_MATRIX<T,3>(fx,fy,fz),id);
    T ypz = fy+fz; if(fabs(ypz)<panic_threshold) ypz=ypz<0?-panic_threshold:panic_threshold;
    T ryz=((dy < 0) != (dz < 0))?((fy!=fz)?dy/(fy-fz):1):0;
    if(dz<0) ryz=1-ryz;
    T Cyz=base.Ey_Ez_y_z(x,y,z,id);
    if(dx<0) Cyz+=base.Exy_Exz_y_z(x,y,z,id)*dx;
    T Syz=(P.x.y+P.x.z)/ypz, Dyz=(2*k-base.Eyz(x,y,z,id))*ryz+(1-ryz)*Cyz;
    yzzy = (Dyz-Syz)/2;
    yzyz = (Dyz+Syz)/2;
}
///#####################################################################
// Function Isotropic_Stress_Derivative_Helper
//#####################################################################
template<class T,int d> void GENERAL_EXTRAPOLATED<T,d>::
Isotropic_Stress_Derivative_Helper(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dP_dF,const int id) const
{
    Isotropic_Stress_Derivative_Helper_Part(F.x.x,F.x.y,F.x.z,id,dP_dF.x0000,dP_dF.x2211,dP_dF.x2121,dP_dF.x2112);
    Isotropic_Stress_Derivative_Helper_Part(F.x.y,F.x.z,F.x.x,id,dP_dF.x1111,dP_dF.x2200,dP_dF.x2020,dP_dF.x2002);
    Isotropic_Stress_Derivative_Helper_Part(F.x.z,F.x.x,F.x.y,id,dP_dF.x2222,dP_dF.x1100,dP_dF.x1010,dP_dF.x1001);
    if(enforce_definiteness) dP_dF.Enforce_Definiteness();
}
//#####################################################################
// Function P_From_Strain_Rate
//#####################################################################
template<class T,int d> MATRIX<T,d> GENERAL_EXTRAPOLATED<T,d>::
P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const int id) const
{
    T id_alpha=(alpha.m?alpha(id):constant_alpha),id_beta=(beta.m?beta(id):constant_beta);
    SYMMETRIC_MATRIX<T,d> strain_rate=F_dot.Symmetric_Part(); // use linear damping because of problems with inverting elements...
    return 2*id_beta*strain_rate+id_alpha*strain_rate.Trace();
}
//#####################################################################
// Function P_From_Strain_Rate_Forces_Size
//#####################################################################
template<class T,int d> int GENERAL_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Forces_Size() const
{
    return sizeof(MATRIX<T,d>)/sizeof(T);
}
//#####################################################################
// Function P_From_Strain_Rate_First_Half
//#####################################################################
template<class T,int d> void GENERAL_EXTRAPOLATED<T,d>::
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
template<class T,int d> MATRIX<T,d> GENERAL_EXTRAPOLATED<T,d>::
P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<const T> aggregate,const int id) const
{
    T id_alpha=(alpha.m?alpha(id):constant_alpha),id_beta=(beta.m?beta(id):constant_beta);
    SYMMETRIC_MATRIX<T,d> strain_rate=(*(const MATRIX<T,d>*)aggregate.Get_Array_Pointer()).Symmetric_Part(); // use linear damping because of problems with inverting elements...
    T sb=sqrt(2*id_beta);
    T dd=sb/TV::dimension;
    T sa=sqrt(id_alpha/TV::dimension+dd*dd)-dd;
    return sb*strain_rate+sa*strain_rate.Trace();
}
namespace PhysBAM{
template class GENERAL_EXTRAPOLATED<float,2>;
template class GENERAL_EXTRAPOLATED<float,3>;
template class GENERAL_EXTRAPOLATED<double,2>;
template class GENERAL_EXTRAPOLATED<double,3>;
}
