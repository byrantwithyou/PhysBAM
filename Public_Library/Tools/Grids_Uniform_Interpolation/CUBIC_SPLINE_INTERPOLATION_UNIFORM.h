//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_SPLINE_INTERPOLATION_UNIFORM 
//#####################################################################
#ifndef __CUBIC_SPLINE_INTERPOLATION_UNIFORM__
#define __CUBIC_SPLINE_INTERPOLATION_UNIFORM__

#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <Tools/Interpolation/CUBIC_MN_INTERPOLATION.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV>
class CUBIC_SPLINE_INTERPOLATION_UNIFORM:public INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP> BASE;
    using BASE::Clamped_Index_Interior_End_Minus_One;

    CUBIC_MN_INTERPOLATION<T,T2> cubic_mn_interpolation;

    CUBIC_SPLINE_INTERPOLATION_UNIFORM();
    ~CUBIC_SPLINE_INTERPOLATION_UNIFORM();

    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const;
    VECTOR<T2,TV::m> Clamped_To_Array_Gradient(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    SYMMETRIC_MATRIX<T2,TV::m> Clamped_To_Array_Hessian(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const PHYSBAM_OVERRIDE;
    T2 From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
