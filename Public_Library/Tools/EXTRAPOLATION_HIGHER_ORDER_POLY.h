//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_HIGHER_ORDER_POLY
//#####################################################################
#ifndef __EXTRAPOLATION_HIGHER_ORDER_POLY__
#define __EXTRAPOLATION_HIGHER_ORDER_POLY__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR.h>
#include <boost/function.hpp>
namespace PhysBAM{

// TODO: Limit to near interface
template<class TV,class T2>
class EXTRAPOLATION_HIGHER_ORDER_POLY:public NONCOPYABLE
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;

public:
    EXTRAPOLATION_HIGHER_ORDER_POLY();
    ~EXTRAPOLATION_HIGHER_ORDER_POLY();
    static void Extrapolate_Node(const GRID<TV>& grid,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& x,int order,int fill_width,int order_reduction_penalty=3);
    static void Extrapolate_Cell(const GRID<TV>& grid,boost::function<bool(const TV_INT& index)> inside_mask,int ghost,ARRAYS_ND_BASE<T2,TV_INT>& x,int order,int fill_width,int order_reduction_penalty=3);
    static void Extrapolate_Face(const GRID<TV>& grid,boost::function<bool(const FACE_INDEX<TV::m>& index)> inside_mask,int ghost,ARRAY<T2,FACE_INDEX<TV::m> >& x,int order,int fill_width,int order_reduction_penalty=3);
//#####################################################################
};
}
#endif
