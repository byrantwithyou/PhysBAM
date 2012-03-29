//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_STENCIL_UNIFORM
//#####################################################################
#ifndef __BASIS_STENCIL_UNIFORM__
#define __BASIS_STENCIL_UNIFORM__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class T,int rank,int d> struct STATIC_POLYNOMIAL;

template<class TV,int d>
struct BASIS_STENCIL_UNIFORM
{
    typedef VECTOR<int,TV::m> TV_INT;typedef typename TV::SCALAR T;
    TV_INT center_offset; // In units of dx/2
    TV dX;

    struct ENTRY
    {
        RANGE<TV_INT> region; // In units of dx/2
        STATIC_POLYNOMIAL<T,TV::m,d> polynomial; // In units of dx
    };

    ARRAY<ENTRY> stencils;

    struct DICED
    {
        TV_INT index_offset;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,d> polynomial;
    };

    ARRAY<DICED> diced;

    BASIS_STENCIL_UNIFORM(TV dx);
    ~BASIS_STENCIL_UNIFORM();

    void Set_Center()
    {center_offset=TV_INT();}

    void Set_Node()
    {center_offset=TV_INT()-1;}

    void Set_Face(int axis)
    {center_offset=-TV_INT::Axis_Vector(axis);}

    void Add_Symmetric_Entry(const ENTRY& e, int mask=-1); // 1=x, 2=y, 4=z
    void Exchange_Axes(int m,int n);
    void Scale_Axes(TV scale);
    void Dice_Stencil();
    void Set_Constant_Stencil();
    void Set_Multilinear_Stencil();
    void Differentiate(int v);

    void Print() const;
};
}
#endif
