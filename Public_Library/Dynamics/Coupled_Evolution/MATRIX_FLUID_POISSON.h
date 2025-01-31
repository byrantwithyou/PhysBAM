//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_POISSON
//#####################################################################
#ifndef __MATRIX_FLUID_POISSON__
#define __MATRIX_FLUID_POISSON__
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>

namespace PhysBAM{
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class TV> class GRID;

template<class TV>
class MATRIX_FLUID_POISSON
{
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;

    const COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    const ARRAY<T,TV_INT>& one_over_rho_c_squared;
public:
    SPARSE_MATRIX_FLAT_MXN<T> poisson;
    ARRAY<int> map;

    MATRIX_FLUID_POISSON(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input,
        const ARRAY<T,TV_INT>& one_over_rho_c_squared_input);
    MATRIX_FLUID_POISSON(const MATRIX_FLUID_POISSON&) = delete;
    void operator=(const MATRIX_FLUID_POISSON&) = delete;

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
