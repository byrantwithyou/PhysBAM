//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM
//#####################################################################
#ifndef __IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM__
#define __IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
#include <string>
namespace PhysBAM{
template<class T> class MPI_UNIFORM_GRID;

template<class TV> class POISSON_COLLIDABLE_UNIFORM;
template<class T_LAPLACE> class HEAT_LAPLACE;
template<class TV> class PROJECTION_DYNAMICS_UNIFORM;
class VIEWER_DIR;

template<class TV>
class IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM:public IMPLICIT_VISCOSITY_UNIFORM<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    typedef ARRAY<int,TV_INT> T_ARRAYS_INT;
    typedef AVERAGING_UNIFORM<TV> T_AVERAGING;
    typedef MPI_UNIFORM_GRID<TV> T_MPI_GRID;
public:
    typedef IMPLICIT_VISCOSITY_UNIFORM<TV> BASE;
    using BASE::face_grid;using BASE::heat_solver;using BASE::u;using BASE::axis;
    using BASE::mpi_grid;using BASE::use_variable_viscosity;

    PROJECTION_DYNAMICS_UNIFORM<TV>& projection;
    ARRAY<T> densities;
    ARRAY<T> viscosities;

    IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM(PROJECTION_DYNAMICS_UNIFORM<TV>& projection_input,const ARRAY<T,TV_INT>& variable_viscosity_input,const ARRAY<T>& densities_input,const ARRAY<T>& viscosities_input,T_MPI_GRID* mpi_grid_input,const int axis_input,bool use_variable_viscosity_input);
    virtual ~IMPLICIT_VISCOSITY_MULTIPHASE_UNIFORM();

//#####################################################################
private:
    void Allocate_Heat_Solver() override;
    void Setup_Viscosity(const T dt) override;
    void Setup_Boundary_Conditions(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities) override;
    void Calculate_Velocity_Jump();
    void Debug_Write(const VIEWER_DIR& viewer_dir);
//#####################################################################
};
}
#endif
