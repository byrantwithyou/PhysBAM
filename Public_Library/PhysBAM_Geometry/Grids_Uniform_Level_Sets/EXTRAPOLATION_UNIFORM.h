//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Eran Guendelman, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_UNIFORM  
//##################################################################### 
//
// Extrapolates the values of u in the normal direction from phi <= 0 to phi > 0,
// overwriting u where phi > 0, or possibly some real cells when using the Isobaric Fix.
//
//#####################################################################
#ifndef __EXTRAPOLATION_UNIFORM__
#define __EXTRAPOLATION_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class EXTRAPOLATION_UNIFORM:public EXTRAPOLATION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef ARRAYS_ND_BASE<T,TV_INT> T_ARRAYS_BASE;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:
    template<class T3> struct REBIND{typedef EXTRAPOLATION_UNIFORM<T_GRID,T3> TYPE;};
    typedef  EXTRAPOLATION<T_GRID,T2> BASE;
    using BASE::band_width;using BASE::boundary;using BASE::isobaric_fix_width;using BASE::Add_To_Heap;using BASE::Remove_Root_From_Heap;
   
    ARRAYS_ND_BASE<T2,TV_INT>& u; // variable to be extrapolated
    const T_ARRAYS_BASE& phi;
    ARRAY<TV_INT>* seed_indices;
    ARRAYS_ND_BASE<bool,TV_INT>* seed_done;
    bool collision_aware_extrapolation;
    const ARRAYS_ND_BASE<VECTOR<bool,TV::m>,TV_INT>* neighbors_visible; // used for collision aware extrapolation
protected:
    T_GRID node_grid;
    T optimization_scale[3]; // dz^2/dy^2, dz^2/dx^2, dy^2/dx^2, optimizations, indexed by missing axis
    TV_INT dimension_start,dimension_end;
    int ghost_cells;
    T small_number;
public:

    EXTRAPOLATION_UNIFORM(const T_GRID& grid,const T_ARRAYS_BASE& phi_input,ARRAYS_ND_BASE<T2,TV_INT>& u_input,const int ghost_cells_input);
    ~EXTRAPOLATION_UNIFORM();

    void Set_Band_Width(const T number_of_cells=(T)5)
    {band_width=number_of_cells*node_grid.dX.Max();}
    
    void Set_Band_Width_To_Entire_Domain()
    {band_width=2*node_grid.domain.Edge_Lengths().Magnitude();}
    
    void Set_Isobaric_Fix_Width(const int number_of_cells=0)
    {isobaric_fix_width=number_of_cells*node_grid.dX.Max();}

    void Set_Custom_Seed_Indices(ARRAY<TV_INT>* seed_indices_input=0)
    {seed_indices=seed_indices_input;}

    void Set_Custom_Seed_Done(ARRAYS_ND_BASE<bool,TV_INT>* seed_done_input=0)
    {seed_done=seed_done_input;}

    void Set_Collision_Aware_Extrapolation(const ARRAYS_ND_BASE<VECTOR<bool,TV::m>,TV_INT>& neighbors_visible_input)
    {collision_aware_extrapolation=true;neighbors_visible=&neighbors_visible_input;}

    void Set_Small_Number(const T small_number_input=1e-8)
    {small_number=small_number_input;}

    bool Neighbor_Visible(const int neighbor_number,const TV_INT& current_index) // neighbor_number between 1 and 3 -- right, top, back
    {return !neighbors_visible->Valid_Index(current_index) || (*neighbors_visible)(current_index)(neighbor_number);}

//#####################################################################
    void Extrapolate(const T time=0,const bool fill_ghost_cells=true);
private:
    void Initialize(const T_ARRAYS_BASE& phi,ARRAYS_ND_BASE<bool,TV_INT>& done,ARRAYS_ND_BASE<bool,TV_INT>& close,ARRAY<TV_INT>& heap,int& heap_length);
    void Update_Close_Point(ARRAYS_ND_BASE<T2,TV_INT>& u,const T_ARRAYS_BASE& phi,const ARRAYS_ND_BASE<bool,TV_INT>& done,const TV_INT& index);
//#####################################################################
}; 
}  
#endif
