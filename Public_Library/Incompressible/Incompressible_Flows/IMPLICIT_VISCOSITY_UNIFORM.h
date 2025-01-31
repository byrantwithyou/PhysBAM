//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_VISCOSITY_UNIFORM
//#####################################################################
#ifndef __IMPLICIT_VISCOSITY_UNIFORM__
#define __IMPLICIT_VISCOSITY_UNIFORM__

#include <Grid_PDE/Interpolation/AVERAGING_UNIFORM.h>
namespace PhysBAM{

template<class TV> class LAPLACE_UNIFORM;
template<class T> class MPI_UNIFORM_GRID;
template<class TV> class BOUNDARY_CONDITION_DOUBLE_FINE;

template<class TV>
class IMPLICIT_VISCOSITY_UNIFORM
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
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
    BOUNDARY_CONDITION_DOUBLE_FINE<TV>* bc_fine=0;
    
    IMPLICIT_VISCOSITY_UNIFORM(LAPLACE_UNIFORM<TV>& elliptic_solver_input,const ARRAY<T,TV_INT>& variable_viscosity_input,const T density_input,const T viscosity_input,T_MPI_GRID* mpi_grid_input,
        const int axis_input,bool use_variable_viscosity_input,bool use_psi_R_input);
    IMPLICIT_VISCOSITY_UNIFORM(const IMPLICIT_VISCOSITY_UNIFORM&) = delete;
    void operator=(const IMPLICIT_VISCOSITY_UNIFORM&) = delete;
    virtual ~IMPLICIT_VISCOSITY_UNIFORM();

//#####################################################################
    void Viscous_Update(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const T dt,const T time,const int maximum_implicit_viscosity_iterations);
    static void Variable_Viscosity_Explicit_Part(const T density,const ARRAY<T,TV_INT>& variable_viscosity,const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,const T dt,const T time);
protected:
    virtual void Allocate_Heat_Solver();
    virtual void Setup_Viscosity(const T dt);
    virtual void Setup_Boundary_Conditions(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);
//#####################################################################
};
}
#endif
