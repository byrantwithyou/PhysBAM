//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KANG_POISSON_VISCOSITY
//#####################################################################
#ifndef __KANG_POISSON_VISCOSITY__
#define __KANG_POISSON_VISCOSITY__

#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class FLUIDS_PARAMETERS_UNIFORM;

template<class TV>
class KANG_POISSON_VISCOSITY
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:

    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters;
    const ARRAY<T,TV_INT>& old_phi;

    bool print_matrix;
    bool test_system;

    KANG_POISSON_VISCOSITY(FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_input,const ARRAY<T,TV_INT>& old_phi_input);
    ~KANG_POISSON_VISCOSITY();

//#####################################################################
    void Project_Fluid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt,T time) const;
    void Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt,T time) const;
    int Cell_Index(const TV_INT& cell) const;
    T Pressure_Jump(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const TV_INT& cell,T dt) const;
//#####################################################################
};
}
#endif
