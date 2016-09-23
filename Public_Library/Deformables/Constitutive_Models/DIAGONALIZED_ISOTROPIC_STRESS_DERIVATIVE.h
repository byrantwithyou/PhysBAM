//#####################################################################
// Copyright 2004-2007, Geoffrey Irving, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE
//#####################################################################
#ifndef __DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE__
#define __DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE__

#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Core/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,int d> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;

template<class T>
class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,1>
{
public:
    T x0000;

    MATRIX<T,1> Differential(const MATRIX<T,1>& dF) const
    {return MATRIX<T,1>(x0000*dF.x[0]);}

//#####################################################################
    void Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,1>& F,const VECTOR<T,1>& dE_ds,const SYMMETRIC_MATRIX<T,1>& dE_dsds); // Not robust
    void Enforce_Definiteness();

    SYMMETRIC_MATRIX<T,1> Get_Hessian_Block() const
    {return SYMMETRIC_MATRIX<T,1>(x0000);}

    void Set_Hessian_Block(const MATRIX<T,1>& M)
    {x0000=M(0,0);}
//#####################################################################
};

template<class T>
class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,2>
{
public:
    T x0000,x1100,x1111; // 2x2 block
    T x1010,x1001; // 2x2 block

    MATRIX<T,2> Differential(const MATRIX<T,2>& dF) const
    {return MATRIX<T,2>(x0000*dF.x[0]+x1100*dF.x[3],x1010*dF.x[1]+x1001*dF.x[2],
        x1001*dF.x[1]+x1010*dF.x[2],x1100*dF.x[0]+x1111*dF.x[3]);}

//#####################################################################
    void Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,2>& F,const VECTOR<T,2>& dE_ds,const SYMMETRIC_MATRIX<T,2>& dE_dsds); // Not robust
    void Enforce_Definiteness();

    SYMMETRIC_MATRIX<T,2> Get_Hessian_Block() const
    {return SYMMETRIC_MATRIX<T,2>(x0000,x1100,x1111);}

    void Set_Hessian_Block(const SYMMETRIC_MATRIX<T,2>& M)
    {x0000=M(0,0);x1100=M(1,0);x1111=M(1,1);}
//#####################################################################
};

template<class T>
class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>
{
public:
    T x0000,x1100,x2200,x1111,x2211,x2222; // 3x3 block
    T x1010,x1001; // 2x2 block
    T x2020,x2002; // 2x2 block
    T x2112,x2121; // 2x2 block

    MATRIX<T,3> Differential(const MATRIX<T,3>& dF) const
    {return MATRIX<T,3>(x0000*dF.x[0]+x1100*dF.x[4]+x2200*dF.x[8],x1010*dF.x[1]+x1001*dF.x[3],x2020*dF.x[2]+x2002*dF.x[6],
        x1001*dF.x[1]+x1010*dF.x[3],x1100*dF.x[0]+x1111*dF.x[4]+x2211*dF.x[8],x2121*dF.x[5]+x2112*dF.x[7],
        x2002*dF.x[2]+x2020*dF.x[6],x2112*dF.x[5]+x2121*dF.x[7],x2200*dF.x[0]+x2211*dF.x[4]+x2222*dF.x[8]);}

//#####################################################################
    void Enforce_Definiteness(const T eigenvalue_clamp_percentage=(T)0,const T epsilon=(T)1e-4);
    void Compute_From_Singular_Value_Derivatives(const DIAGONAL_MATRIX<T,3>& F,const VECTOR<T,3>& dE_ds,const SYMMETRIC_MATRIX<T,3>& dE_dsds); // Not robust
    SYMMETRIC_MATRIX<T,3> Get_Hessian_Block() const
    {return SYMMETRIC_MATRIX<T,3>(x0000,x1100,x2200,x1111,x2211,x2222);}

    void Set_Hessian_Block(const SYMMETRIC_MATRIX<T,3>& M)
    {x0000=M(0,0);x1100=M(1,0);x2200=M(2,0);x1111=M(1,1);x2211=M(2,1);x2222=M(2,2);}
//#####################################################################
};
}
#endif
