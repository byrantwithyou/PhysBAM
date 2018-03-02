//#####################################################################
// Copyright 2002-2006, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET
//#####################################################################
#ifndef __LEVELSET__
#define __LEVELSET__

#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/INTERVAL.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Vectors/SCALAR_POLICY.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <cassert>
#include <cfloat>
namespace PhysBAM{

template<class TV> struct INTERPOLATION_POLICY;
template<class TV> class LEVELSET_CALLBACKS; // TODO: invalid dependency
template<class TV> struct BOUNDARY_POLICY;
template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID;
template<class TV,class T2> class BOUNDARY;
template<class TV> class LEVELSET;

template<class TV>
class LEVELSET
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef INTERPOLATION_UNIFORM<TV,T> T_INTERPOLATION_SCALAR;
    typedef INTERPOLATION_UNIFORM<TV,TV> T_INTERPOLATION_VECTOR;
    typedef LINEAR_INTERPOLATION_UNIFORM<TV,TV> T_LINEAR_INTERPOLATION_VECTOR;
    typedef FACE_LOOKUP_UNIFORM<TV> T_FACE_LOOKUP;
public:
    typedef int HAS_UNTYPED_READ_WRITE;
    T small_number;
    T_INTERPOLATION_SCALAR *interpolation,*curvature_interpolation,*secondary_interpolation;
    T_INTERPOLATION_VECTOR* normal_interpolation;
    T max_time_step;
    T sigma;
    T half_band_width;

    BOUNDARY<TV,T>* boundary;
    LEVELSET_CALLBACKS<TV>* levelset_callbacks;
    const ARRAY<bool,FACE_INDEX<TV::m> >* face_velocities_valid_mask_current;
//protected:
    BOUNDARY<TV,T>& boundary_default;
    static LINEAR_INTERPOLATION_UNIFORM<TV,T> interpolation_default;
    static T_LINEAR_INTERPOLATION_VECTOR normal_interpolation_default;
    ARRAY<bool,TV_INT> valid_mask_current;
    ARRAY<bool,TV_INT> valid_mask_next;
    GRID<TV>& grid;
    ARRAY<T,TV_INT>& phi;
    ARRAY<TV,TV_INT>* normals;
    ARRAY<T,TV_INT> *curvature;
    ARRAY<INTERVAL<T>,TV_INT> *cell_range;
    int number_of_ghost_cells;

    LEVELSET(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,const int number_of_ghost_cells_input=3);
    LEVELSET(const LEVELSET&) = delete;
    void operator=(const LEVELSET&) = delete;
    ~LEVELSET();

    void Set_Small_Number(const T small_number_input=1e-8)
    {small_number=small_number_input;}

    void Set_Max_Time_Step(const T max_time_step_input=1e8)
    {max_time_step=max_time_step_input;}

protected:
    void Initialize_Levelset_Grid_Values(const GRID<TV>& grid_input)
    {}
public:

    void Set_Band_Width(const T number_of_cells=6)
    {half_band_width=number_of_cells*grid.dX.Max()/2;}

    void Initialize_Valid_Masks(const GRID<TV>& grid_input)
    {
        valid_mask_current.Resize(grid_input.Cell_Indices(3),use_init,true);
        valid_mask_next.Resize(grid_input.Cell_Indices(3),no_init);
    }

    void Set_Custom_Boundary(BOUNDARY<TV,T>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Custom_Secondary_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {secondary_interpolation=&interpolation_input;}

    void Set_Custom_Normal_Interpolation(T_INTERPOLATION_VECTOR& normal_interpolation_input)
    {normal_interpolation=&normal_interpolation_input;}

    void Set_Custom_Curvature_Interpolation(T_INTERPOLATION_SCALAR& curvature_interpolation_input)
    {curvature_interpolation=&curvature_interpolation_input;}

    void Set_Levelset_Callbacks(LEVELSET_CALLBACKS<TV>& levelset_callbacks_input)
    {levelset_callbacks=&levelset_callbacks_input;}

    void Set_Face_Velocities_Valid_Mask(const ARRAY<bool,FACE_INDEX<TV::m> >* face_velocities_valid_mask_current_input)
    {
        face_velocities_valid_mask_current=face_velocities_valid_mask_current_input;
    }

    void Initialize_Levelset_Grid_Values()
    {Initialize_Levelset_Grid_Values(grid);}

    T Phi(const TV& location) const
    {return interpolation->Clamped_To_Array(grid,phi,location);}

    T Phi_Secondary(const TV& location) const
    {return secondary_interpolation->Clamped_To_Array(grid,phi,location);}

    T Extended_Phi(const TV& location) const
    {TV clamped_location(grid.Clamp(location));
    T phi_value=interpolation->Clamped_To_Array(grid,phi,clamped_location);
    T magnitude_squared=(clamped_location-location).Magnitude_Squared();
    if(magnitude_squared) phi_value=sqrt(magnitude_squared+sqr(max((T)0,phi_value)));
    return phi_value;}

    T Curvature(const TV& location) const // later add finite difference for curvature, like in normal above
    {assert(curvature);return curvature_interpolation->Clamped_To_Array(grid,*curvature,location);}

    bool Lazy_Inside(const TV& clamped_X,const T contour_value=0) const
    {assert(cell_range);TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    const INTERVAL<T>& range=(*cell_range)(index);
    if(range.min_corner>contour_value) return false;
    else if(range.max_corner<=contour_value) return true;
    return interpolation->From_Base_Node(grid,phi,clamped_X,index)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& clamped_X,T& phi_value,const T contour_value=0) const
    {assert(cell_range);TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    if((*cell_range)(index).min_corner>contour_value) return false;
    phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    return phi_value<=contour_value;}

    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const
    {if(grid.domain.Lazy_Inside_Half_Open(unclamped_X)) return Lazy_Inside(unclamped_X,contour_value);
    if(contour_value<=0) return false;
    T sqr_contour_value=sqr(contour_value);
    TV clamped_X(grid.Clamp(unclamped_X));
    T magnitude_squared=(unclamped_X-clamped_X).Magnitude_Squared();
    if(magnitude_squared>sqr_contour_value) return false;
    TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    T phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    return magnitude_squared+sqr(max((T)0,phi_value))<=sqr_contour_value;}

    // phi_value only guaranteed to be set if function returns true
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const
    {assert(cell_range);TV clamped_X(grid.Clamp(unclamped_X));
    T magnitude_squared=(unclamped_X-clamped_X).Magnitude_Squared();
    if(contour_value <= 0){if(magnitude_squared > 0) return false;}
    else if(magnitude_squared > sqr(contour_value)) return false;
    TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    if((*cell_range)(index).min_corner>contour_value) return false;
    phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    if(magnitude_squared>0) phi_value=sqrt(magnitude_squared+sqr(max((T)0,phi_value)));
    return phi_value<=contour_value;}

    bool Lazy_Outside(const TV& clamped_X,const T contour_value=0) const
    {return !Lazy_Inside(clamped_X,contour_value);}

    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const
    {return !Lazy_Inside_Extended_Levelset(unclamped_X,contour_value);}

    // phi_value only guaranteed to be set if function returns true
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const
    {assert(cell_range);TV clamped_X(grid.Clamp(unclamped_X));
    T magnitude_squared=(unclamped_X-clamped_X).Magnitude_Squared();
    TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    if(magnitude_squared==0 && (*cell_range)(index).max_corner<=contour_value) return false;
    phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    if(magnitude_squared>0) phi_value=sqrt(magnitude_squared+sqr(max((T)0,phi_value)));
    return phi_value>contour_value;}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,grid,phi);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,grid,phi);}

    static TV Normal_At_Node(const GRID<TV>& grid,const ARRAY<T,TV_INT>& phi,const TV_INT& index)
    {TV N;for(int d=0;d<TV::m;d++){TV_INT a(index),b(index);a(d)--;b(d)++;N(d)=(phi(b)-phi(a))*grid.one_over_dX(d);}return N.Normalized();}

    T Compute_Curvature(const TV_INT& index) const
    {return Compute_Curvature(phi,index);}

//#####################################################################
    T CFL(const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities) const;
    T CFL(const ARRAY<TV,TV_INT>& velocity) const;
    TV Iterative_Find_Interface(TV left,TV right,const int iterations=3) const;
    void Compute_Gradient(ARRAY<TV,TV_INT>& gradient,const T time=0) const;
    void Compute_Normals(const T time=0);
    TV Gradient(const TV& location) const;
    TV Normal(const TV& location) const;
    TV Extended_Normal(const TV& location) const;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const ARRAY<T,TV_INT>& phi_input,const TV_INT& index) const;
    TV Gradient(const ARRAY<T,TV_INT>& phi_input,const TV_INT& index) const;
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true);
    T Compute_Curvature(const ARRAY<T,TV_INT>& phi_input,const TV_INT& index) const;
    void Compute_Curvature(const T time=0);
    T Compute_Curvature(const TV& location) const;
    void Fast_Marching_Method(const T time=0,const T stopping_distance=0,const ARRAY<TV_INT>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    void Get_Signed_Distance_Using_FMM(ARRAY<T,TV_INT>& signed_distance,const T time=0,const T stopping_distance=0,const ARRAY<TV_INT>* seed_indices=0,
        const bool add_seed_indices_for_ghost_cells=false,int process_sign=0);
    T Approximate_Surface_Size(const T interface_thickness=3,const T time=0) const;
    VECTOR<T,TV::m-1> Principal_Curvatures(const TV& X) const;
//#####################################################################
};
}
#endif
