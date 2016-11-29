//#####################################################################
// Copyright 2010, Craig Schroeder, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED
//#####################################################################
#ifndef __MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED__
#define __MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED__
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>
#include <Dynamics/Coupled_Evolution/MATRIX_FLUID_INTERPOLATION_BASE.h>

namespace PhysBAM{

template<class TV>
class MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED:public MATRIX_FLUID_INTERPOLATION_BASE<TV>
{
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;typedef VECTOR<int,d> TV_INT;
    using MATRIX_FLUID_INTERPOLATION_BASE<TV>::index_map;

    struct STENCIL
    {
        VECTOR<VECTOR<PAIR<int,T>,4>,d> s;
    };

    ARRAY<STENCIL> stencils;

public:
    struct ENTRY
    {
        FACE_INDEX<d> face_index;
        int side;
        TV X;
    };

    ARRAY<ENTRY> entries;

    MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input);
    virtual ~MATRIX_FLUID_INTERPOLATION_EXTRAPOLATED();

//#####################################################################
    COUPLING_CONSTRAINT_ID Number_Of_Constraints() const override;
    void Compute(int ghost_cells) override;
    void Times_Add(const ARRAY<T>& faces,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const override;
    void Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,ARRAY<T>& faces) const override;
    void Print() const override;
    void Print_Each_Matrix(int n) const override;
    void Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_FLUID_MASS<TV>& fluid_mass) const override;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const override;
//#####################################################################
};
}
#endif
