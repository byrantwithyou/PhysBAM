//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KANG_POISSON_VISCOSITY
//#####################################################################
#ifndef __KANG_POISSON_VISCOSITY__
#define __KANG_POISSON_VISCOSITY__

#include <Grid_Tools/Grids/FACE_INDEX.h>
namespace PhysBAM{

template<class TV> class FLUIDS_PARAMETERS_UNIFORM;

template<class TV>
class KANG_POISSON_VISCOSITY
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:

    FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters;
    const ARRAY<T,TV_INT>& old_phi;

    bool print_matrix;
    bool test_system;
    mutable ARRAY<T,FACE_INDEX<TV::m> > face_velocities_ghost;
    ARRAY<bool,TV_INT> psi_D;
    ARRAY<bool,FACE_INDEX<TV::dimension> > psi_N;
    ARRAY<T,TV_INT> psi_D_value;
    ARRAY<T,FACE_INDEX<TV::dimension> > psi_N_value;

    KANG_POISSON_VISCOSITY(FLUIDS_PARAMETERS_UNIFORM<TV>& fluids_parameters_input,const ARRAY<T,TV_INT>& old_phi_input);
    ~KANG_POISSON_VISCOSITY();

//#####################################################################
    void Project_Fluid(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt) const;
    void Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,T dt,bool implicit) const;
    void Apply_Viscosity(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,int axis,T dt,bool implicit) const;
    T Pressure_Jump(const TV_INT& cell,T dt) const;
    MATRIX<T,TV::m> Viscosity_Jump(const TV_INT& cell) const;
    MATRIX<T,TV::m> Viscosity_Jump(const FACE_INDEX<TV::m>& cell) const;
    T Face_Phi(const FACE_INDEX<TV::m>& face) const;
//#####################################################################
};
}
#endif
