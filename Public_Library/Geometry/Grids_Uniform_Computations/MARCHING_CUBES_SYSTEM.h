//#####################################################################
// Copyright 2012, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//#####################################################################
// Class MARCHING_CUBES_SYSTEM
//#####################################################################
#ifndef __MARCHING_CUBES_SYSTEM__
#define __MARCHING_CUBES_SYSTEM__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Matrices/MATRIX.h>
#include <Geometry/Grids_Uniform_Computations/MARCHING_CUBES_VECTOR.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>

namespace PhysBAM{
//#####################################################################
// Class MARCHING_CUBES_SYSTEM
//#####################################################################
template<class TV>
class MARCHING_CUBES_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
public:
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::SCALAR T;
    typedef MARCHING_CUBES_VECTOR<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef MATRIX<T,TV::m> TM;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;

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

    static void Compute_Active_List(ARRAY<VECTOR<int,TV::m+1> >& active_list,const HASHTABLE<VECTOR<int,2>,T_SURFACE*>& interface,const ARRAY<int>& reverse_index_map);
    T Set_Matrix_And_Rhs(MARCHING_CUBES_VECTOR<TV>& rhs,const ARRAY<VECTOR<int,TV::m+1> >& active_list,const ARRAY<int>& index_map,const ARRAY<int>& reverse_index_map,ARRAY_VIEW<const TV> X);
    T Set_Matrix_Block_And_Rhs(const VECTOR<int,TV::m+1> index,const VECTOR<TV,TV::m+1> particles,INDIRECT_ARRAY<ARRAY<TV>,VECTOR<int,TV::m+1>&>* rhs);
    static T Set_Rhs(MARCHING_CUBES_VECTOR<TV>& rhs,const ARRAY<VECTOR<int,TV::m+1> >& active_list,const ARRAY<int>& index_map,ARRAY_VIEW<const TV> X);
    static T Compute_Energy(const ARRAY<VECTOR<int,TV::m+1> >& active_list,ARRAY_VIEW<const TV> X);
    static void Test_System(const ARRAY<VECTOR<int,TV::m+1> >& active_list,const ARRAY<int>& index_map,const ARRAY<int>& reverse_index_map);

};
}
#endif
