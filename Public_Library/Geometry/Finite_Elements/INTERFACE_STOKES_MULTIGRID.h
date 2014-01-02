//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERFACE_STOKES_MULTIGRID
//#####################################################################
#ifndef __INTERFACE_STOKES_MULTIGRID__
#define __INTERFACE_STOKES_MULTIGRID__
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Geometry/Finite_Elements/INTERFACE_STOKES_SYSTEM_COLOR.h>

namespace PhysBAM{

template<class TV>
class INTERFACE_STOKES_MULTIGRID:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    typedef INTERFACE_STOKES_SYSTEM_VECTOR_COLOR<TV> T_VECTOR;

    struct LEVEL
    {
        INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss;
        ARRAY<SPARSE_MATRIX_FLAT_MXN<T> > pressure_poisson; // per color
        T_VECTOR tmp0,tmp1,tmp2;

        void Interior_Smoother(T_VECTOR& z,const T_VECTOR& x) const; // z should be initial guess
        void Boundary_Smoother(T_VECTOR& z,const T_VECTOR& x) const; // z should be initial guess

        LEVEL()
            :iss(0)
        {
        }
    };

    ARRAY<LEVEL> levels;

    INTERFACE_STOKES_MULTIGRID(int num_levels,INTERFACE_STOKES_SYSTEM_COLOR<TV>* iss);
    ~INTERFACE_STOKES_MULTIGRID();

    void Construct_Level(int l);

    void Coarsen(T_VECTOR& z,const T_VECTOR& x,int fine_level) const;
    void Prolongation(T_VECTOR& z,const T_VECTOR& x,int fine_level) const;
    void Exact_Solve(T_VECTOR& z,const T_VECTOR& x) const;

    void Apply_Preconditioner(T_VECTOR& z,const T_VECTOR& x,bool initial_guess);
    void Update();

//#####################################################################
};
}
#endif
