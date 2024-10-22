//#####################################################################
// Copyright 2010, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_VISCOUS_FORCES
//#####################################################################
#ifndef __MATRIX_VISCOUS_FORCES__
#define __MATRIX_VISCOUS_FORCES__
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/PAIR.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Dynamics/Coupled_Evolution/VISCOUS_FORCE_ID.h>

namespace PhysBAM{
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;

template<class TV>
class MATRIX_VISCOUS_FORCES:public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    enum WORKAROUND{d=TV::m};
private:
    const GRID<TV>& grid;
    const COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    struct ENTRY
    {
        T weight;
        int face_index;
        VISCOUS_FORCE_ID viscous_id;
        
        ENTRY()
        {}
        
        ENTRY(T weight_input,int face_index_input,VISCOUS_FORCE_ID id)
            :weight(weight_input),face_index(face_index_input),viscous_id(id)
        {}
    };

    ARRAY<ENTRY> entries;
    VISCOUS_FORCE_ID last_id;

public:
    MATRIX_VISCOUS_FORCES(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input);
    MATRIX_VISCOUS_FORCES(const MATRIX_VISCOUS_FORCES&) = delete;
    void operator=(const MATRIX_VISCOUS_FORCES&) = delete;
    virtual ~MATRIX_VISCOUS_FORCES();

//#####################################################################
    void Compute(const T dt,const ARRAY<bool,FACE_INDEX<d> >& psi_N,T mu);
    void Times_Add(const ARRAY<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const;
    void Times(const ARRAY<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const;
    void Transpose_Times_Add(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,ARRAY<T>& velocities) const;
    void Transpose_Times(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,ARRAY<T>& velocities) const;
    VISCOUS_FORCE_ID Viscous_Forces_Size() const;
    void Test_Matrix() const;
    void Print_Each_Matrix(int n) const;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const override;
//#####################################################################
};
}
#endif
