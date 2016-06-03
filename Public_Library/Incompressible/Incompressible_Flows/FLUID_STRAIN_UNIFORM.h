//#####################################################################
// Copyright 2002-2006, Doug Enright, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andy Lutomirksi, Jonathan Su, Paul-James White.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_STRAIN_UNIFORM  
//#####################################################################
#ifndef __FLUID_STRAIN_UNIFORM__
#define __FLUID_STRAIN_UNIFORM__

#include <Tools/Advection/ADVECTION_FORWARD.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_UNIFORM_FORWARD.h>
#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Incompressible/Boundaries/BOUNDARY_FORWARD.h>
#include <Incompressible/Incompressible_Flows/FLUID_STRAIN.h>
namespace PhysBAM{

template<class T> class EXTERNAL_STRAIN_ADJUSTMENT;
template<class TV> class LEVELSET_MULTIPLE;
template<class TV,class T2> class BOUNDARY;

template<class TV>
class FLUID_STRAIN_UNIFORM:public FLUID_STRAIN<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<SYMMETRIC_MATRIX<T,TV::m> ,TV_INT> T_ARRAYS_SYMMETRIC_MATRIX;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,T> T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM_BETA<TV,SYMMETRIC_MATRIX<T,TV::m> > T_ADVECTION_SEMI_LAGRANGIAN_SYMMETRIC_MATRIX;
public:
    using FLUID_STRAIN<T>::viscosity_index;using FLUID_STRAIN<T>::strainrate_time;using FLUID_STRAIN<T>::elastic_modulus;
    using FLUID_STRAIN<T>::plasticity_alpha;using FLUID_STRAIN<T>::plasticity_gamma;

    GRID<TV> grid;
    T_ARRAYS_SYMMETRIC_MATRIX e; // strain tensor
    BOUNDARY<TV,SYMMETRIC_MATRIX<T,TV::m> >* e_boundary;
    ADVECTION<TV,SYMMETRIC_MATRIX<T,TV::m> >* e_advection;
    EXTERNAL_STRAIN_ADJUSTMENT<T>* external_strain_adjustment;
private:               
    BOUNDARY<TV,SYMMETRIC_MATRIX<T,TV::m> >& e_boundary_default;
    T_ADVECTION_SEMI_LAGRANGIAN_SYMMETRIC_MATRIX& e_advection_default;
    mutable bool cfl_called;
public:

    FLUID_STRAIN_UNIFORM(const GRID<TV>& grid_input);
    ~FLUID_STRAIN_UNIFORM();

    void Initialize_Grid(const GRID<TV>& grid_input)
    {assert(grid_input.Is_MAC_Grid());grid=grid_input;e.Resize(grid.Domain_Indices());}

    void Set_Custom_Boundary(BOUNDARY<TV,SYMMETRIC_MATRIX<T,TV::m> >& e_boundary_input)
    {e_boundary=&e_boundary_input;}

    void Set_Custom_Advection(ADVECTION<TV,SYMMETRIC_MATRIX<T,TV::m> >& e_advection_input)
    {e_advection=&e_advection_input;}

    void Set_External_Strain_Adjustment(EXTERNAL_STRAIN_ADJUSTMENT<T>& external_strain_adjustment_input)
    {external_strain_adjustment=&external_strain_adjustment_input;}

//#####################################################################
    void Update_Strain_Equation(const T dt,const T time,const T density,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,
        const ARRAY<T,TV_INT>& phi_ghost,const int number_of_ghost_cells);
    void Update_Strain_Equation_Multiphase(const T dt,const T time,const T density,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,
        const LEVELSET_MULTIPLE<TV>& levelset,const int region,const int number_of_ghost_cells);
private:
    void Update_Strain_Equation_Helper_Cell_Centered(const T dt,const T time,const T density,const T heaviside_bandwidth,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,
        ARRAY<TV,TV_INT>& V,const ARRAY<T,TV_INT>& phi_ghost,const int number_of_ghost_cells);
public:
    void Extrapolate_Strain_Across_Interface(ARRAY<T,TV_INT>& phi_ghost,const T band_width=3);
    T CFL(const T density) const;
//#####################################################################
};

// not implemented in one dimension
template<class T>
class FLUID_STRAIN_UNIFORM<VECTOR<T,1> >:public FLUID_STRAIN<T>
{
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
public:
//#####################################################################
    typedef ARRAY<SYMMETRIC_MATRIX<T,TV::m> ,TV_INT> T_ARRAYS_SYMMETRIC_MATRIX;

    T_ARRAYS_SYMMETRIC_MATRIX e; // strain tensor
    BOUNDARY<TV,SYMMETRIC_MATRIX<T,TV::m> >* e_boundary;

    FLUID_STRAIN_UNIFORM(const GRID<TV>& grid_input){PHYSBAM_NOT_IMPLEMENTED();}
    void Initialize_Grid(const GRID<TV>& grid_input){PHYSBAM_NOT_IMPLEMENTED();}
    void Set_Custom_Boundary(BOUNDARY<TV,SYMMETRIC_MATRIX<T,TV::m> >& e_boundary_input){e_boundary=&e_boundary_input;}
    void Update_Strain_Equation(const T dt,const T time,const T density,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,
        const ARRAY<T,VECTOR<int,1> >& phi_ghost,const int number_of_ghost_cells){PHYSBAM_NOT_IMPLEMENTED();}
    void Update_Strain_Equation_Multiphase(const T dt,const T time,const T density,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,
        const LEVELSET_MULTIPLE<TV>& levelset,const int region,const int number_of_ghost_cells){PHYSBAM_NOT_IMPLEMENTED();}
    void Extrapolate_Strain_Across_Interface(ARRAY<T,VECTOR<int,1> >& phi_ghost,const T band_width=3){PHYSBAM_NOT_IMPLEMENTED();}
    T CFL(const T density) const{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
};

}
#endif
