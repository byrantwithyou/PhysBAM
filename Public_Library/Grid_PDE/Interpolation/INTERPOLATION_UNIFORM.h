//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERPOLATION_UNIFORM__
#define __INTERPOLATION_UNIFORM__

#include <Core/Arrays/ARRAY.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Grid_Tools/Grids/BLOCK_UNIFORM.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM_FORWARD.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV>
class INTERPOLATION_UNIFORM
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    STATIC_ASSERT((is_same<typename GRID<TV>::GRID_TAG,UNIFORM_TAG<TV> >::value));
public:
    template<class T3> struct REBIND{typedef INTERPOLATION_UNIFORM<TV,T3,T_FACE_LOOKUP> TYPE;};

    INTERPOLATION_UNIFORM() = default;
    INTERPOLATION_UNIFORM(const INTERPOLATION_UNIFORM&) = delete;
    void operator=(const INTERPOLATION_UNIFORM&) = delete;
    virtual ~INTERPOLATION_UNIFORM()
    {}

    static TV_INT Clamped_Index(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& location,int start_offset=0,int end_offset=0) 
    {return clamp(grid.Index(location),u.domain.min_corner+start_offset,u.domain.max_corner-(end_offset+1));}

    static TV_INT Clamped_Index_End_Minus_One(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& location)
    {return Clamped_Index(grid,u,location,0,1);}

    static TV_INT Clamped_Index_Interior(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& location)
    {return Clamped_Index(grid,u,location,1,1);}

    static TV_INT Clamped_Index_Interior_End_Minus_One(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& location)
    {return Clamped_Index(grid,u,location,1,2);}

    void Populate_New_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const GRID<TV>& grid_new,ARRAYS_ND_BASE<T2,TV_INT>& u_new)
    {for(CELL_ITERATOR<TV> iterator(grid_new,u_new.domain);iterator.Valid();iterator.Next()){ // CELL ITERATOR works for nodal
        u_new(iterator.Cell_Index())=Clamped_To_Array(grid,u,grid_new.X(iterator.Cell_Index()));}}

    void Populate_New_Array(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& u,const GRID<TV>& grid_new,ARRAY<T,FACE_INDEX<TV::m> >& u_new)
    {FACE_LOOKUP_UNIFORM<TV> lookup(u);
    for(FACE_ITERATOR<TV> iterator(grid_new);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        u_new.Component(axis)(iterator.Face_Index())=Clamped_To_Array_Face_Component(axis,grid,lookup,iterator.Location());}}

    T2 Clamped_To_Array_Cell(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
    {return Clamped_To_Array(grid,u,X);}

    T2 Clamped_To_Array_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
    {return Clamped_To_Array(grid,u,X);}

    TV Clamped_To_Array_Face(const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(grid.Is_MAC_Grid());return From_Block_Face(grid,BLOCK_UNIFORM<TV>(grid,X,u.Number_Of_Ghost_Cells()),u,X);}

    T Clamped_To_Array_Face_Component(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(grid.Is_MAC_Grid());return From_Block_Face_Component(axis,grid,BLOCK_UNIFORM<TV>(grid,X,u.Number_Of_Ghost_Cells()),u,X);}
    
    ARRAY<PAIR<FACE_INDEX<TV::m>,T> > Clamped_To_Array_Face_Component_Weights(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(grid.Is_MAC_Grid());return From_Block_Face_Component_Weights(axis,grid,BLOCK_UNIFORM<TV>(grid,X,u.Number_Of_Ghost_Cells()),u,X);}

    VECTOR<T2,2> Extrema_Clamped_To_Array_Cell(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const
    {return Extrema_Clamped_To_Array(grid,u_min,u_max,X);}

    VECTOR<T2,2> Extrema_Clamped_To_Array_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const
    {return Extrema_Clamped_To_Array(grid,u_min,u_max,X);}

    VECTOR<TV,2> Extrema_Clamped_To_Array_Face(const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {assert(grid.Is_MAC_Grid());return Extrema_From_Block_Face(grid,BLOCK_UNIFORM<TV>(grid,X,u_min.Number_Of_Ghost_Cells()),u_min,u_max,X);}

    VECTOR<T,2> Extrema_Clamped_To_Array_Face_Component(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {assert(grid.Is_MAC_Grid());return Extrema_From_Block_Face_Component(axis,grid,BLOCK_UNIFORM<TV>(grid,X,u_min.Number_Of_Ghost_Cells()),u_min,u_max,X);}

//#####################################################################
    virtual T2 Clamped_To_Array_No_Extrema(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<T2,TV::m> Clamped_To_Array_Gradient(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual SYMMETRIC_MATRIX<T2,TV::m> Clamped_To_Array_Hessian(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TV From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual ARRAY<PAIR<FACE_INDEX<TV::m>,T> > From_Block_Face_Component_Weights(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<T2,2> Extrema_Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const 
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<T2,2> Extrema_From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X,const TV_INT& index) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<TV,2> Extrema_From_Block_Face(const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<T,2> Extrema_From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,
        const typename T_FACE_LOOKUP::LOOKUP& u_min,const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
