//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC_INTERPOLATION_UNIFORM 
//#####################################################################
#ifndef __QUADRATIC_INTERPOLATION_UNIFORM__
#define __QUADRATIC_INTERPOLATION_UNIFORM__

#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> >
class QUADRATIC_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_Interior_End_Minus_One;

public:
    int ghost_cells;

    QUADRATIC_INTERPOLATION_UNIFORM(): ghost_cells(0) {}
    virtual ~QUADRATIC_INTERPOLATION_UNIFORM() {}

    // T2 Clamped_To_Array_No_Extrema(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override
    // {return Clamped_To_Array(grid,u,X);}

    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const override;
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const override;
    TV_INT Base_Index(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const;
    TV_INT Base_Index_Face(const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,int axis,const TV& X) const;

    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
//    VECTOR<T2,2> Extrema_Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const  override;
//    VECTOR<T2,2> Extrema_From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X,const TV_INT& index) const override;
    T From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override;
//    ARRAY<PAIR<FACE_INDEX<TV::m>,T> > From_Block_Face_Component_Weights(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override;
//#####################################################################
};
}
#endif
