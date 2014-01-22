//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_2X2.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
#include <Deformables/Constitutive_Models/ISOTROPIC_CONSTITUTIVE_MODEL.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,int d> ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
ISOTROPIC_CONSTITUTIVE_MODEL()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,int d> ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
~ISOTROPIC_CONSTITUTIVE_MODEL()
{
}
//#####################################################################
// Function dP_From_dF
//#####################################################################
template<class T,int d> MATRIX<T,d> ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
dP_From_dF(const MATRIX<T,d>& dF,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,const T scale,const int simplex) const
{
    return scale*dPi_dF.Differential(dF);
}
//#####################################################################
// Function Update_State_Dependent_Auxiliary_Variables
//#####################################################################
template<class T,int d> void ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
Update_State_Dependent_Auxiliary_Variables(const DIAGONAL_MATRIX<T,d>& F,const int simplex)
{
}
//#####################################################################
// Function Energy_Density
//#####################################################################
template<class T,int d> T ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
Energy_Density(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Report_Diagnostics
//#####################################################################
template<class T> static void
Report_Diagnostics(const char* str,T x,T y)
{
    char buff[1000];
    sprintf(buff, "%6s %8.5f %8.5f (%8.5f)", str, x, y, fabs(x-y));
    LOG::cout<<buff<<std::endl;
}
//#####################################################################
// Function Report_Diagnostics
//#####################################################################
template<class T> static void
Report_Diagnostics(const char* str,T g0,T g1,T dg0,T dg1,T e)
{
    T av=(dg0+dg1)/2;
    T dif=(g1-g0)/e;
    Report_Diagnostics(str, av, dif);
}
//#####################################################################
// Function Report_Diagnostics
//#####################################################################
template<class T> static void
Report_Diagnostics(const DIAGONAL_MATRIX<T,2>& F,T E0,T* E,const DIAGONAL_MATRIX<T,2>& P0,const DIAGONAL_MATRIX<T,2>* P,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>& dPi_dF0,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>* dPi_dF,T e)
{
    Report_Diagnostics("Px",E0,E[0],P0.x.x,P[0].x.x,e);
    Report_Diagnostics("Py",E0,E[1],P0.x.y,P[1].x.y,e);

    Report_Diagnostics("Hxx",P0.x.x,P[0].x.x,dPi_dF0.x1111,dPi_dF[0].x1111,e);
    Report_Diagnostics("Hyy",P0.x.y,P[1].x.y,dPi_dF0.x2222,dPi_dF[1].x2222,e);
    Report_Diagnostics("Hxy",P0.x.x,P[1].x.x,dPi_dF0.x2211,dPi_dF[1].x2211,e);
    T ss1=sqr(F.x.x),ss2=sqr(F.x.y);
    T s12=1/(ss1-ss2);
    Report_Diagnostics("x2112",dPi_dF0.x2112,(-P0.x.y*F.x.x+P0.x.x*F.x.y)*s12);
    Report_Diagnostics("x2121",dPi_dF0.x2121,(-P0.x.y*F.x.y+P0.x.x*F.x.x)*s12);
}
//#####################################################################
// Function Report_Diagnostics
//#####################################################################
template<class T> static void
Report_Diagnostics(const DIAGONAL_MATRIX<T,3>& F,T E0,T* E,const DIAGONAL_MATRIX<T,3>& P0,const DIAGONAL_MATRIX<T,3>* P,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF0,
    const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>* dPi_dF,T e)
{
    Report_Diagnostics("Px",E0,E[0],P0.x.x,P[0].x.x,e);
    Report_Diagnostics("Py",E0,E[1],P0.x.y,P[1].x.y,e);
    Report_Diagnostics("Pz",E0,E[2],P0.x.z,P[2].x.z,e);

    Report_Diagnostics("Hxx",P0.x.x,P[0].x.x,dPi_dF0.x1111,dPi_dF[0].x1111,e);
    Report_Diagnostics("Hyy",P0.x.y,P[1].x.y,dPi_dF0.x2222,dPi_dF[1].x2222,e);
    Report_Diagnostics("Hzz",P0.x.z,P[2].x.z,dPi_dF0.x3333,dPi_dF[2].x3333,e);
    Report_Diagnostics("Hxy",P0.x.x,P[1].x.x,dPi_dF0.x2211,dPi_dF[1].x2211,e);
    Report_Diagnostics("Hxz",P0.x.x,P[2].x.x,dPi_dF0.x3311,dPi_dF[2].x3311,e);
    Report_Diagnostics("Hyz",P0.x.y,P[2].x.y,dPi_dF0.x3322,dPi_dF[2].x3322,e);
    T ss1=sqr(F.x.x),ss2=sqr(F.x.y),ss3=sqr(F.x.z);
    T s12=1/(ss1-ss2),s13=1/(ss1-ss3),s23=1/(ss2-ss3);
    Report_Diagnostics("x2112",dPi_dF0.x2112,(-P0.x.y*F.x.x+P0.x.x*F.x.y)*s12);
    Report_Diagnostics("x2121",dPi_dF0.x2121,(-P0.x.y*F.x.y+P0.x.x*F.x.x)*s12);
    Report_Diagnostics("x3113",dPi_dF0.x3113,(-P0.x.z*F.x.x+P0.x.x*F.x.z)*s13);
    Report_Diagnostics("x3131",dPi_dF0.x3131,(-P0.x.z*F.x.z+P0.x.x*F.x.x)*s13);
    Report_Diagnostics("x3223",dPi_dF0.x3223,(-P0.x.z*F.x.y+P0.x.y*F.x.z)*s23);
    Report_Diagnostics("x3232",dPi_dF0.x3232,(-P0.x.z*F.x.z+P0.x.y*F.x.y)*s23);
}
//#####################################################################
// Function Test
//#####################################################################
template<class T,int d> void ISOTROPIC_CONSTITUTIVE_MODEL<T,d>::
Test(const DIAGONAL_MATRIX<T,d>& F,const int simplex) const
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    T e=(T)1e-6;
    DIAGONAL_MATRIX<T,d> dF;

    T E0=Energy_Density(F,simplex),E[d];
    DIAGONAL_MATRIX<T,d> P0=P_From_Strain(F,1,simplex),P[d];
    DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d> dPi_dF0,dPi_dF[d];
    Isotropic_Stress_Derivative(F,dPi_dF0,simplex);
    for(int i=0;i<d;i++){
        DIAGONAL_MATRIX<T,d> Fd=F+DIAGONAL_MATRIX<T,d>(VECTOR<T,d>::Axis_Vector(i+1)*e);
        E[i]=Energy_Density(Fd,simplex);
        P[i]=P_From_Strain(Fd,1,simplex);
        Isotropic_Stress_Derivative(Fd,dPi_dF[i],simplex);}
    Report_Diagnostics(F,E0,E,P0,P,dPi_dF0,dPi_dF,e);
}
namespace PhysBAM{
template class ISOTROPIC_CONSTITUTIVE_MODEL<float,2>;
template class ISOTROPIC_CONSTITUTIVE_MODEL<float,3>;
template class ISOTROPIC_CONSTITUTIVE_MODEL<double,2>;
template class ISOTROPIC_CONSTITUTIVE_MODEL<double,3>;
}
