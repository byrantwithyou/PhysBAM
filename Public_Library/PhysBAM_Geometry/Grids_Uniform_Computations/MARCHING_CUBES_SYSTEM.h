//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
// Class MARCHING_CUBES_SYSTEM
//#####################################################################
#ifndef __MARCHING_CUBES_SYSTEM__
#define __MARCHING_CUBES_SYSTEM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>

namespace PhysBAM{
//#####################################################################
// Class MARCHING_CUBES_VECTOR
//#####################################################################
template<class TV>
class MARCHING_CUBES_VECTOR:public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;

public:

    ARRAY<TV> x; 

    MARCHING_CUBES_VECTOR(){}
    virtual ~MARCHING_CUBES_VECTOR(){}

    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE{
        x+=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x;
        return *this;}

    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE{
        x-=debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x;
        return *this;}

    BASE& operator*=(const T a) PHYSBAM_OVERRIDE{
        x*=a; return *this;}

    void Copy(const T c1,const BASE& bv1) PHYSBAM_OVERRIDE{
        x=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).x*c1;}

    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE{
        x=debug_cast<const MARCHING_CUBES_VECTOR&>(bv1).x*c1+debug_cast<const MARCHING_CUBES_VECTOR&>(bv2).x;}

    int Raw_Size() const PHYSBAM_OVERRIDE
    {return x.m*TV::m;}

    T& Raw_Get(int i) PHYSBAM_OVERRIDE{
        return x(i/TV::m)(i%TV::m);}

    KRYLOV_VECTOR_BASE<T>* Clone_Default() const PHYSBAM_OVERRIDE{
        MARCHING_CUBES_VECTOR* V=new MARCHING_CUBES_VECTOR;
        V->x.Resize(x.m); return V;}

    void Resize(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE{
        x.Resize(debug_cast<const MARCHING_CUBES_VECTOR&>(bv).x.m);}
};
//#####################################################################
// Class MARCHING_CUBES_SYSTEM
//#####################################################################
template<class TV>
class MARCHING_CUBES_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
public:

    typedef typename TV::SCALAR T;
    typedef MARCHING_CUBES_VECTOR<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef MATRIX<T,TV::m> TM;

    struct BLOCK
    {
        VECTOR<int,TV::m+1> index;
        VECTOR<VECTOR<TM,TV::m+1>,TV::m+1> matrix;
    };
    
    ARRAY<BLOCK> blocks;

    MARCHING_CUBES_SYSTEM():BASE(false,false){}
    virtual ~MARCHING_CUBES_SYSTEM(){}

//#####################################################################

    void Multiply(const KRYLOV_VECTOR_BASE<T>& bv_input,KRYLOV_VECTOR_BASE<T>& bv_result) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE;

    void Project(KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE{}
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE{}
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bv) const PHYSBAM_OVERRIDE{Project(bv);}
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& br,KRYLOV_VECTOR_BASE<T>& bz) const PHYSBAM_OVERRIDE{}

//#####################################################################

    T Set_Matrix_Block_And_Rhs(const VECTOR<int,TV::m+1> index,const VECTOR<TV,TV::m+1> particles,INDIRECT_ARRAY<ARRAY<TV>,VECTOR<int,TV::m+1>&> rhs);
};
}
#endif
