//#####################################################################
// Copyright 2015
//#####################################################################
// Class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM 
//#####################################################################
#ifndef __CUBIC_MONOTONIC_INTERPOLATION_UNIFORM__
#define __CUBIC_MONOTONIC_INTERPOLATION_UNIFORM__

#include <Tools/Interpolation/CUBIC_MN_INTERPOLATION.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP>
class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_Interior_End_Minus_One;

    CUBIC_MN_INTERPOLATION<T,T2> cubic_mn_interpolation;

    CUBIC_MONOTONIC_INTERPOLATION_UNIFORM();
    ~CUBIC_MONOTONIC_INTERPOLATION_UNIFORM();

    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const override;
    T2 Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const;
    TV_INT Base_Index(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const;
    TV_INT Base_Index_Face(const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,int axis,const TV& X) const;
    T From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const override;
//#####################################################################
};
}
#endif
