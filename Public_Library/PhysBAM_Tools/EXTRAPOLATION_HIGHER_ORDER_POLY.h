//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_HIGHER_ORDER_POLY
//#####################################################################
#ifndef __EXTRAPOLATION_HIGHER_ORDER_POLY__
#define __EXTRAPOLATION_HIGHER_ORDER_POLY__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

// TODO: Limit to near interface
template<class TV,class T2>
class EXTRAPOLATION_HIGHER_ORDER_POLY:public NONCOPYABLE
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;

public:
    EXTRAPOLATION_HIGHER_ORDER_POLY();
    ~EXTRAPOLATION_HIGHER_ORDER_POLY();
    struct MASK
    {
        virtual bool Inside(const TV_INT& index)=0;
    };
    struct MASK_FACE
    {
        virtual bool Inside(const FACE_INDEX<TV::m>& index)=0;
    };

    static void Extrapolate_Node(const GRID<TV>& grid,MASK& inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& x,int order,int fill_width,T order_reduction_penalty=3);
    static void Extrapolate_Cell(const GRID<TV>& grid,MASK& inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& x,int order,int fill_width,T order_reduction_penalty=3);
    static void Extrapolate_Face(const GRID<TV>& grid,MASK_FACE& inside_mask,int ghost,ARRAY<T2,FACE_INDEX<TV::m> >& x,int order,int fill_width,T order_reduction_penalty=3);
//#####################################################################
};
}
#endif
