//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONALIZED_STRESS_DERIVATIVE
//##################################################################### 
#ifndef __DIAGONALIZED_STRESS_DERIVATIVE__
#define __DIAGONALIZED_STRESS_DERIVATIVE__

#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,int d> class DIAGONALIZED_STRESS_DERIVATIVE;

template<class T>
class DIAGONALIZED_STRESS_DERIVATIVE<T,2>
{
    typedef VECTOR<T,2> TV;
public:
    DIAGONAL_MATRIX<T,2> F;
    SYMMETRIC_MATRIX<T,2> S;
    SYMMETRIC_MATRIX<T,3> dSdC;

    MATRIX<T,2> Differential(const MATRIX<T,2>& dF) const
    {SYMMETRIC_MATRIX<T,2> dC((F*dF).Twice_Symmetric_Part());
    VECTOR<T,3> ds=dSdC*VECTOR<T,3>(dC.x00,dC.x11,sqrt((T)2)*dC.x10);
    SYMMETRIC_MATRIX<T,2> dS(ds.x,sqrt((T)1/2)*ds.z,ds.y);
    return dF*S+F*dS;}

    void Enforce_Definiteness()
    {S=S.Positive_Definite_Part();dSdC=dSdC.Positive_Definite_Part();}
};

template<class T>
class DIAGONALIZED_STRESS_DERIVATIVE<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    DIAGONAL_MATRIX<T,3> F;
    SYMMETRIC_MATRIX<T,3> S,dSdC_d,dSdC_s;
    MATRIX<T,3> dSdC_ds;

    MATRIX<T,3> Differential(const MATRIX<T,3>& dF) const
    {SYMMETRIC_MATRIX<T,3> dC((F*dF).Twice_Symmetric_Part());
    VECTOR<T,3> dC_d(dC.x00,dC.x11,dC.x22),dC_s(sqrt((T)2)*dC.x10,sqrt((T)2)*dC.x20,sqrt((T)2)*dC.x21);
    VECTOR<T,3> dS_d(dSdC_d*dC_d+dSdC_ds*dC_s),dS_s(dSdC_ds.Transpose_Times(dC_d)+dSdC_s*dC_s);
    SYMMETRIC_MATRIX<T,3> dS(dS_d.x,sqrt((T)1/2)*dS_s.x,sqrt((T)1/2)*dS_s.y,dS_d.y,sqrt((T)1/2)*dS_s.z,dS_d.z);
    return dF*S+F*dS;}

    void Enforce_Definiteness()
    {S=S.Positive_Definite_Part();
    DIAGONAL_MATRIX<T,3> D_d,D_s;MATRIX<T,3> V_d,V_s;dSdC_d.Fast_Solve_Eigenproblem(D_d,V_d);dSdC_s.Fast_Solve_Eigenproblem(D_s,V_s);MATRIX<T,3> M(V_d.Transpose_Times(dSdC_ds*V_s));
    D_d.x00=max(D_d.x00,abs(M(0,0))+abs(M(0,1))+abs(M(0,2)));D_d.x11=max(D_d.x11,abs(M(1,0))+abs(M(1,1))+abs(M(1,2)));D_d.x22=max(D_d.x22,abs(M(2,0))+abs(M(2,1))+abs(M(2,2)));
    D_s.x00=max(D_s.x00,abs(M(0,0))+abs(M(1,0))+abs(M(2,0)));D_s.x11=max(D_s.x11,abs(M(0,1))+abs(M(1,1))+abs(M(2,1)));D_s.x22=max(D_s.x22,abs(M(0,2))+abs(M(1,2))+abs(M(2,2)));
    dSdC_d=SYMMETRIC_MATRIX<T,3>::Conjugate(V_d,D_d);dSdC_s=SYMMETRIC_MATRIX<T,3>::Conjugate(V_s,D_s);}
};
}
#endif
