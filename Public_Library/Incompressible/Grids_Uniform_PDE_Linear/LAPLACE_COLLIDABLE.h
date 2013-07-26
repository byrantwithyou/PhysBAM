//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_COLLIDABLE
//#####################################################################
#ifndef __LAPLACE_COLLIDABLE__
#define __LAPLACE_COLLIDABLE__

#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class TV>
class LAPLACE_COLLIDABLE
{
    typedef typename TV::SCALAR T;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;

public:
    bool second_order_cut_cell_method;
    T second_order_cut_cell_threshold;
    LEVELSET<TV>* levelset; // used in second order accurate cut cell method
    T_FACE_ARRAYS_SCALAR u_interface; // interface boundary condition - 2nd order method

protected:
    LAPLACE_COLLIDABLE()
        :second_order_cut_cell_method(false),second_order_cut_cell_threshold((T)1e-3),levelset(0)
    {}
public:
    virtual ~LAPLACE_COLLIDABLE()
    {}

    virtual void Use_External_Level_Set(LEVELSET<TV>& cell_centered_levelset)
    {levelset=&cell_centered_levelset;}

    virtual void Set_Up_Second_Order_Cut_Cell_Method(const bool use_second_order_cut_cell_method_input=true)
    {second_order_cut_cell_method=use_second_order_cut_cell_method_input;}

//#####################################################################
};
}
#endif

