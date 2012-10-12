//#####################################################################
// Copyright 2002-2006, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET
//#####################################################################
#ifndef __LEVELSET__
#define __LEVELSET__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <cassert>
#include <cfloat>
namespace PhysBAM{

template<class T_GRID> struct COLLISION_GEOMETRY_COLLECTION_POLICY;
template<class T_GRID> struct INTERPOLATION_POLICY;
template<class T_GRID> class LEVELSET_CALLBACKS;
template<class T_GRID> struct BOUNDARY_POLICY;
template<class T_GRID> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID;
template<class T_GRID,class T2> class BOUNDARY_UNIFORM;
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class T,class T_GRID=GRID<VECTOR<T,1> > >
class LEVELSET:public NONCOPYABLE
{
    STATIC_ASSERT((IS_SCALAR<T>::value));
    typedef typename T_GRID::VECTOR_T TV;typedef VECTOR<int,TV::m> TV_INT;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
    typedef typename REBIND<T_INTERPOLATION_SCALAR,TV>::TYPE T_INTERPOLATION_VECTOR;
    typedef typename REBIND<T_LINEAR_INTERPOLATION_SCALAR,TV>::TYPE T_LINEAR_INTERPOLATION_VECTOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef ARRAY<bool,FACE_INDEX<TV::m> > T_FACE_ARRAYS_BOOL;
public:
    T small_number;
    T_INTERPOLATION_SCALAR *interpolation,*curvature_interpolation,*secondary_interpolation;
    T_INTERPOLATION_VECTOR* normal_interpolation;
    T max_time_step;
    bool curvature_motion;
    T sigma;
public:
    bool refine_fmm_initialization_with_iterative_solver;
    int fmm_initialization_iterations;
    T fmm_initialization_iterative_tolerance;
    T fmm_initialization_iterative_drift_fraction;

    BOUNDARY_UNIFORM<T_GRID,T>* boundary;
    LEVELSET_CALLBACKS<T_GRID>* levelset_callbacks;
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >* collision_body_list;
    const T_FACE_ARRAYS_BOOL* face_velocities_valid_mask_current;
    bool collision_aware_signed_distance,clamp_phi_with_collision_bodies;
//protected:
    BOUNDARY_UNIFORM<T_GRID,T>& boundary_default;
    static T_LINEAR_INTERPOLATION_SCALAR interpolation_default;
    T_INTERPOLATION_SCALAR *collision_aware_interpolation_plus,*collision_aware_interpolation_minus,*collision_unaware_interpolation;
    static T_LINEAR_INTERPOLATION_VECTOR normal_interpolation_default;
    T collidable_phi_replacement_value;
public:
    ARRAY<bool,TV_INT> valid_mask_current;
    ARRAY<bool,TV_INT> valid_mask_next;
protected:

    LEVELSET();
    virtual ~LEVELSET();

public:
    void Set_Small_Number(const T small_number_input=1e-8)
    {small_number=small_number_input;}

    void Set_Max_Time_Step(const T max_time_step_input=1e8)
    {max_time_step=max_time_step_input;}

    void Set_Curvature_Motion(const T sigma_input=1) // times by the grid spacing
    {curvature_motion=true;sigma=sigma_input;assert(sigma >= 0);}

    void Initialize_FMM_Initialization_Iterative_Solver(const bool refine_fmm_initialization_with_iterative_solver_input=true,const int fmm_initialization_iterations_input=10,
        const T fmm_initialization_iterative_tolerance_input=1e-2,const T fmm_initialization_iterative_drift_fraction_input=.1)
    {refine_fmm_initialization_with_iterative_solver=refine_fmm_initialization_with_iterative_solver_input;fmm_initialization_iterations=fmm_initialization_iterations_input;
    fmm_initialization_iterative_tolerance=fmm_initialization_iterative_tolerance_input;fmm_initialization_iterative_drift_fraction=fmm_initialization_iterative_drift_fraction_input;}

    virtual void Initialize_Levelset_Grid_Values()=0; // should call version below

protected:
    void Initialize_Levelset_Grid_Values(const T_GRID& grid_input)
    {}
public:

    void Initialize_Valid_Masks(const T_GRID& grid_input)
    {
        valid_mask_current.Resize(grid_input.Cell_Indices(3),true,true,true);
        valid_mask_next.Resize(grid_input.Cell_Indices(3),false);
    }

    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,T>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Custom_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {interpolation=&interpolation_input;}

    void Set_Custom_Secondary_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {secondary_interpolation=&interpolation_input;}

    void Set_Custom_Normal_Interpolation(T_INTERPOLATION_VECTOR& normal_interpolation_input)
    {normal_interpolation=&normal_interpolation_input;}

    void Set_Custom_Curvature_Interpolation(T_INTERPOLATION_SCALAR& curvature_interpolation_input)
    {curvature_interpolation=&curvature_interpolation_input;}

    void Enable_Collision_Aware_Interpolation(const int sign)
    {PHYSBAM_ASSERT(!collision_unaware_interpolation);collision_unaware_interpolation=interpolation;
    interpolation=sign==1?collision_aware_interpolation_plus:collision_aware_interpolation_minus;}

    void Disable_Collision_Aware_Interpolation()
    {interpolation=collision_unaware_interpolation;collision_unaware_interpolation=0;}

    void Set_Levelset_Callbacks(LEVELSET_CALLBACKS<T_GRID>& levelset_callbacks_input)
    {levelset_callbacks=&levelset_callbacks_input;}

    void Set_Collision_Body_List(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_body_list_input,const bool set_secondary_interpolation=false);

    void Set_Face_Velocities_Valid_Mask(const T_FACE_ARRAYS_BOOL* face_velocities_valid_mask_current_input)
    {
        face_velocities_valid_mask_current=face_velocities_valid_mask_current_input;
    }
};
}
#endif
