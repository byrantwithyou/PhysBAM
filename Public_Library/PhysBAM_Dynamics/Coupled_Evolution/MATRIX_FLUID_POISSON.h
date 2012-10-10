//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_POISSON
//#####################################################################
#ifndef __MATRIX_FLUID_POISSON__
#define __MATRIX_FLUID_POISSON__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class TV> class GRID;

template<class TV>
class MATRIX_FLUID_POISSON:public NONCOPYABLE
{
    enum WORKAROUND {d=TV::dimension};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;

    const COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    const T_ARRAYS_SCALAR& one_over_rho_c_squared;
public:
    SPARSE_MATRIX_FLAT_NXN<T> poisson;
    ARRAY<int> map;

    MATRIX_FLUID_POISSON(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input,
        const T_ARRAYS_SCALAR& one_over_rho_c_squared_input);

//#####################################################################
    void Compute(const SPARSE_MATRIX_FLAT_MXN<T>& gradient,const ARRAY<T>& one_over_fluid_mass,const T dt,const bool use_preconditioner);
    void Compute_Preconditioner();
    void Apply_Preconditioner(ARRAY<T>& pressure) const;
    void Times_Add(const ARRAY<T>& pressure_in,ARRAY<T>& pressure_out) const;
    void Times(const ARRAY<T>& pressure_in,ARRAY<T>& pressure_out) const;
    void Test_Matrix(const bool print_matrix=false) const;
    void Print_Each_Matrix(int n) const;
//#####################################################################
};
}
#endif
