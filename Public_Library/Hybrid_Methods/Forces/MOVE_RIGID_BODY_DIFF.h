//#####################################################################
// Copyright 2017.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MOVE_RIGID_BODY_DIFF
//#####################################################################
#ifndef __MOVE_RIGID_BODY_DIFF__
#define __MOVE_RIGID_BODY_DIFF__

#include <Core/Matrices/FRAME.h>
#include <Tools/Tensors/TENSOR.h>
namespace PhysBAM{
template<class TV> class TWIST;

template<class TV>
struct MOVE_RIGID_BODY_DIFF
{
    typedef typename TV::SCALAR T;
    FRAME<TV> frame;
    TENSOR<T,TV::m,TV::m,TV::SPIN::m> tensor;
    MATRIX<T,TV::m> R;

    void Compute(const FRAME<TV>& frame0,const TWIST<TV>& dt_twist);
    TV Frame_Times(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m>& dZdL,
        MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const;
    TV Frame_Inverse_Times(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m>& dZdL,
        MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const;
    TV Rotate(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const;
    TV Inverse_Rotate(TV v,MATRIX<T,TV::m>& dZdv,MATRIX<T,TV::m,TV::SPIN::m>& dZdA) const;
};
}
#endif
