//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC_INTERPOLATION_UNIFORM 
//#####################################################################
#ifndef __QUADRATIC_INTERPOLATION_UNIFORM__
#define __QUADRATIC_INTERPOLATION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID> >
class QUADRATIC_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    //typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef ARRAYS_ND_BASE<T2,VECTOR<int,TV::dimension> > T_ARRAYS_T2;
public:
    typedef INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_Interior_End_Minus_One;

public:
    int ghost_cells;

    QUADRATIC_INTERPOLATION_UNIFORM();
    ~QUADRATIC_INTERPOLATION_UNIFORM();

    // T2 Clamped_To_Array_No_Extrema(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE
    // {return Clamped_To_Array(grid,u,X);}

    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE
    {return From_Base_Node_Weights_Helper(grid,u,X,index);}

    T2 From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE
    {return From_Base_Node_Helper(grid,u,X,index);}

    T2 From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    T2 From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    T2 From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    ARRAY<PAIR<TV_INT,T> > From_Base_Node_Weights_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
    TV_INT Base_Index(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const;
    TV_INT Base_Index_Face(const T_GRID& grid,const typename T_FACE_LOOKUP::LOOKUP& u,int axis,const TV& X) const;

    T2 Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const PHYSBAM_OVERRIDE;
//    VECTOR<T2,2> Extrema_Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X) const  PHYSBAM_OVERRIDE;
//    VECTOR<T2,2> Extrema_From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;
    T From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE;
//    ARRAY<PAIR<FACE_INDEX<T_GRID::VECTOR_T::dimension>,T> > From_Block_Face_Component_Weights(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const;
    T2 From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const;
    T2 From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const;
//#####################################################################
};
}
#endif
