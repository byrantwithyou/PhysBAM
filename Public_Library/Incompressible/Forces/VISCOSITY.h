//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISCOSITY
//#####################################################################
#ifndef __VISCOSITY__
#define __VISCOSITY__

#include <Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{
template<class TV> class LAPLACE_UNIFORM;

template<class TV>
class VISCOSITY:public INCOMPRESSIBLE_FLUIDS_FORCES<TV>
{
    
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
    LAPLACE_UNIFORM<TV>& elliptic_solver;
public:
    T density;
    T viscosity;
    bool implicit_viscosity;
    bool use_explicit_part_of_implicit_viscosity;
    const ARRAY<T,TV_INT>& variable_viscosity;
    bool use_variable_viscosity;
    int maximum_implicit_viscosity_iterations;
    bool use_psi_R;

    VISCOSITY(LAPLACE_UNIFORM<TV>& elliptic_solver_input,const ARRAY<T,TV_INT>& variable_viscosity_input,const T density_input,const T viscosity_input,bool implicit_viscosity_input,
        bool use_explicit_part_of_implicit_viscosity_input,bool use_variable_viscosity_input,int maximum_implicit_viscosity_iterations_input,bool use_psi_R_input);
    virtual ~VISCOSITY();

//#####################################################################
    void Add_Explicit_Forces(const GRID<TV>& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Add_Implicit_Forces_Before_Projection(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Add_Implicit_Forces_Projection(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Grids(const GRID<TV>& grid) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
