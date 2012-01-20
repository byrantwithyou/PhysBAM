//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_CONSERVATIVE_UNIFORM
//#####################################################################
#ifndef __ADVECTION_CONSERVATIVE_UNIFORM__
#define __ADVECTION_CONSERVATIVE_UNIFORM__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/CUBIC_MN_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> //  T_AVERAGING=AVERAGING_UNIFORM<T_GRID>, T_INTERPOLATION=LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2>
class ADVECTION_CONSERVATIVE_UNIFORM:public ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
public:
    using ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>::Update_Advection_Equation_Cell;
    template<class T3> struct REBIND{typedef ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T3,T_AVERAGING,typename T_INTERPOLATION::template REBIND<T3>::TYPE> TYPE;};
    template<class T_INTERPOLATION_2> struct REBIND_INTERPOLATION{typedef ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION_2> TYPE;};

    bool use_second_order;
    bool clamp_weights;
    bool use_collision_object;
    T2 total_mass_lost;
    T2 total_mass_gained;
    T2 total_potential_lost;
    T2 total_potential_gained;
    RANGE<TV> ghost_box;
    SPHERE<TV> collision_object;
    VECTOR<VECTOR<bool,2>,TV::dimension> solid_walls_hack_axis;
    ARRAY<T,TV_INT> momentum_lost;
    bool find_closest_point;
    int num_cells;
    int number_of_ghost_cells;
    int cfl;
    T density;

    ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >,FACE_INDEX<TV::dimension> > weights_to; //w ij
    ARRAY<ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >,FACE_INDEX<TV::dimension> > weights_from; //w ki
    ARRAY<T,FACE_INDEX<TV::dimension> > sum_jc; //wbar_jc
    ARRAY<T,TV_INT> sum_jc_cell; //wbar_jc
    int num_iterations,num_diffusion_iterations;
    PARTICLE_LEVELSET_UNIFORM<T_GRID>* pls;
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid;
    VECTOR<VECTOR<bool,2>,TV::dimension> mpi_boundary;
    int evenodd,evenodd_cell;
    ARRAY<T,TV_INT>* smoke_density;

//#####################################################################
    ADVECTION_CONSERVATIVE_UNIFORM();
    void Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost=0,const T_ARRAYS_T2* Z_max_ghost=0,T_ARRAYS_T2* Z_min=0,T_ARRAYS_T2* Z_max=0);
    void Update_Advection_Equation_Cell(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost=0,const T_ARRAYS_T2* Z_max_ghost=0,T_ARRAYS_T2* Z_min=0,T_ARRAYS_T2* Z_max=0);
    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max);
     void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max);
    void Update_Advection_Equation_Face_Lookup(T_GRID& grid,T_ARRAYS_SCALAR& phi,T_ARRAYS_SCALAR& phi2,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max);
    void Clamp_Weights_To_Grid(const RANGE<TV_INT>& inside_domain,ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& X);
    void Clamp_Weights_To_Objects(const GRID<TV>& grid,ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& weights);
    void Clamp_Weights_To_Objects(const GRID<TV>& grid,ARRAY<PAIR<TV_INT,T> >& weights);
    void Clamp_X_To_Grid(const GRID<TV>& grid,TV& X);
    bool Is_Outside(const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face);
    bool Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face);
    bool Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const TV_INT& index);
    void Clean_Weights(ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& weights);
    void Cell_Diffusion_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z);
    void Face_Diffusion_Helper(const GRID<TV>& grid,FACE_INDEX<TV::dimension>& first_face_index,FACE_INDEX<TV::dimension>& second_face_index,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside);
    void Face_Diffusion_Helper(FACE_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside=0);
    void Face_Diffusion_Helper(CELL_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside=0);
    void Face_Diffusion(const T_GRID& grid,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,T_BOUNDARY& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside=0);
    void Cell_Diffusion(const T_GRID& grid,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,T_BOUNDARY_T2& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum);
//#####################################################################
};
}
#endif
