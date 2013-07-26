//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM__
#define __LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM__

#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
namespace PhysBAM{

template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV,class T2>
class LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list;
    const ARRAY<bool,TV_INT>* cell_valid_value_mask;
    T2 default_cell_replacement_value;
    bool extrapolate_to_invalid_cell_values;
    LINEAR_INTERPOLATION_UNIFORM<TV,T2> interpolation;

    LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input,const ARRAY<bool,TV_INT>* cell_valid_value_mask_input,
        const T2& default_cell_replacement_value_input,const bool extrapolate_to_invalid_cell_values_input=true);
    virtual ~LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM();

    void Set_Default_Replacement_Value(const T2& default_cell_replacement_value_input)
    {default_cell_replacement_value=default_cell_replacement_value_input;}

    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAY<T2,TV_INT>& u,const TV& X,bool& interpolation_valid) const
    {return From_Base_Node(grid,u,X,INTERPOLATION_UNIFORM<TV,T2>::Clamped_Index_End_Minus_One(grid,u,X),interpolation_valid);}
    
//#####################################################################
    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    // TODO: this only works for cells, which would ideally be enforced at compile time
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& cell) const PHYSBAM_OVERRIDE;
    // TODO: this only works for cells, which would ideally be enforced at compile time
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& cell,bool& interpolation_valid,VECTOR<T2,2>* extrema=0) const;
    VECTOR<T2,2> Extrema_Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const PHYSBAM_OVERRIDE;
    VECTOR<T2,2> Extrema_From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X,const TV_INT& cell) const PHYSBAM_OVERRIDE;
protected:
    T2 Trace_And_Get_Value(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& cell,TV_INT* cells,int& valid_mask,const int mask_value=0) const;
    virtual T2 Invalid_Value_Replacement(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& cell,const TV& intersection_point,
        const COLLISION_GEOMETRY_ID body_id,const int aggregate_id,bool& valid,const T ray_t_max=0) const;
    void Extrapolate_To_Invalid_Values(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& cell,T2* values,int& valid_mask) const;
//#####################################################################
};
}
#endif
