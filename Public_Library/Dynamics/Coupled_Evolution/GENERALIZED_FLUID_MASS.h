//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERALIZED_FLUID_MASS
//#####################################################################
#ifndef __GENERALIZED_FLUID_MASS__
#define __GENERALIZED_FLUID_MASS__
#include <Core/Arrays/ARRAY.h>
#include <Core/Matrices/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;

template<class TV>
class GENERALIZED_FLUID_MASS:public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;

    const COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    const ARRAY<T,FACE_INDEX<TV::m> >& beta; // 1/rho
    const ARRAY<T>& constrained_beta;

public:
    ARRAY<T> one_over_fluid_mass_at_faces;

    GENERALIZED_FLUID_MASS(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input,const ARRAY<T,FACE_INDEX<TV::m> >& beta_input,const ARRAY<T>& constrained_beta_input);
    GENERALIZED_FLUID_MASS(const GENERALIZED_FLUID_MASS&) = delete;
    void operator=(const GENERALIZED_FLUID_MASS&) = delete;
    virtual ~GENERALIZED_FLUID_MASS();

    void Compute();
    void Inverse_Times(const ARRAY<T>& faces_in,ARRAY<T>& faces_out) const;
    void Inverse_Times_Add(const ARRAY<T>& faces_in,ARRAY<T>& faces_out) const;
    void Second_Order_Mass_Correction(const ARRAY<T,TV_INT>& phi);    
    void Compute_For_Two_Phase_Pressure_Jump(const ARRAY<T,TV_INT>& phi,const ARRAY<T,TV_INT>& density);

    void Print_Each_Matrix(int n) const;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const override;
};
}
#endif
