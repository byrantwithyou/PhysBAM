//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_VISCOSITY_UNIFORM
//#####################################################################
#ifndef __IMPLICIT_VISCOSITY_UNIFORM__
#define __IMPLICIT_VISCOSITY_UNIFORM__

#include <Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class LAPLACE_UNIFORM;

template<class T_GRID>
class IMPLICIT_VISCOSITY_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef AVERAGING_UNIFORM<T_GRID> T_AVERAGING;
    typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;

protected:
    LAPLACE_UNIFORM<T_GRID>& elliptic_solver;
    const T_ARRAYS_SCALAR& variable_viscosity;
    T density;
    T viscosity;
    T_MPI_GRID* mpi_grid;
    const int axis;
    T_GRID face_grid; // a mac grid with the faces of the axis as cells
    LAPLACE_UNIFORM<T_GRID>* heat_solver;
    T_ARRAYS_SCALAR u;
    bool use_variable_viscosity;
    bool use_psi_R;
public:

    IMPLICIT_VISCOSITY_UNIFORM(LAPLACE_UNIFORM<T_GRID>& elliptic_solver_input,const T_ARRAYS_SCALAR& variable_viscosity_input,const T density_input,const T viscosity_input,T_MPI_GRID* mpi_grid_input,
        const int axis_input,bool use_variable_viscosity_input,bool use_psi_R_input);
    virtual ~IMPLICIT_VISCOSITY_UNIFORM();

//#####################################################################
    void Viscous_Update(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T dt,const T time,const int maximum_implicit_viscosity_iterations);
    static void Variable_Viscosity_Explicit_Part(const T density,const T_ARRAYS_SCALAR& variable_viscosity,const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T dt,const T time);
protected:
    virtual void Allocate_Heat_Solver();
    virtual void Setup_Viscosity(const T dt);
    virtual void Setup_Boundary_Conditions(const T_FACE_ARRAYS_SCALAR& face_velocities);
//#####################################################################
};
}
#endif
