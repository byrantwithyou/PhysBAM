//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_INTERPOLATION_BASE
//#####################################################################
#ifndef __MATRIX_FLUID_INTERPOLATION_BASE__
#define __MATRIX_FLUID_INTERPOLATION_BASE__
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>

namespace PhysBAM{

template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;
template<class TV> class GENERALIZED_FLUID_MASS;

template<class TV>
class MATRIX_FLUID_INTERPOLATION_BASE:public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;
protected:
    COLLISION_AWARE_INDEX_MAP<TV>& index_map;

public:
    MATRIX_FLUID_INTERPOLATION_BASE(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input);
    MATRIX_FLUID_INTERPOLATION_BASE(const MATRIX_FLUID_INTERPOLATION_BASE&) = delete;
    void operator=(const MATRIX_FLUID_INTERPOLATION_BASE&) = delete;
    virtual ~MATRIX_FLUID_INTERPOLATION_BASE();

//#####################################################################
    virtual COUPLING_CONSTRAINT_ID Number_Of_Constraints() const=0;
    virtual void Compute(int ghost_cells)=0;
    virtual void Times_Add(const ARRAY<T>& faces,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const=0;
    void Times(const ARRAY<T>& faces,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const;
    virtual void Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,ARRAY<T>& faces) const=0;
    void Transpose_Times(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,ARRAY<T>& faces) const;
    virtual void Print() const=0;
    void Test_Matrix() const;
    virtual void Print_Each_Matrix(int n) const=0;
    virtual void Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_FLUID_MASS<TV>& fluid_mass) const=0;
//#####################################################################
};
}
#endif
