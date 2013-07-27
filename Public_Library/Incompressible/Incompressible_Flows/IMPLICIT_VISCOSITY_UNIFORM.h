//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_VISCOSITY_UNIFORM
//#####################################################################
#ifndef __IMPLICIT_VISCOSITY_UNIFORM__
#define __IMPLICIT_VISCOSITY_UNIFORM__

#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class TV> class LAPLACE_UNIFORM;
template<class T> class MPI_UNIFORM_GRID;

template<class TV>
class IMPLICIT_VISCOSITY_UNIFORM:public NONCOPYABLE
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef typename ARRAY<T,TV_INT>::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef AVERAGING_UNIFORM<TV> T_AVERAGING;
    typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;

protected:
    LAPLACE_UNIFORM<TV>& elliptic_solver;
    const ARRAY<T,TV_INT>& variable_viscosity;
    T density;
    T viscosity;
    T_MPI_GRID* mpi_grid;
    const int axis;
    GRID<TV> face_grid; // a mac grid with the faces of the axis as cells
    LAPLACE_UNIFORM<TV>* heat_solver;
    ARRAY<T,TV_INT> u;
    bool use_variable_viscosity;
    bool use_psi_R;
public:

    IMPLICIT_VISCOSITY_UNIFORM(LAPLACE_UNIFORM<TV>& elliptic_solver_input,const ARRAY<T,TV_INT>& variable_viscosity_input,const T density_input,const T viscosity_input,T_MPI_GRID* mpi_grid_input,
        const int axis_input,bool use_variable_viscosity_input,bool use_psi_R_input);
    virtual ~IMPLICIT_VISCOSITY_UNIFORM();

//#####################################################################
    void Viscous_Update(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T dt,const T time,const int maximum_implicit_viscosity_iterations);
    static void Variable_Viscosity_Explicit_Part(const T density,const ARRAY<T,TV_INT>& variable_viscosity,const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T dt,const T time);
protected:
    virtual void Allocate_Heat_Solver();
    virtual void Setup_Viscosity(const T dt);
    virtual void Setup_Boundary_Conditions(const T_FACE_ARRAYS_SCALAR& face_velocities);
//#####################################################################
};
}
#endif
