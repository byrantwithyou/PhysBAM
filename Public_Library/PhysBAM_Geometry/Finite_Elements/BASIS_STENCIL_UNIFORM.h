//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BASIS_STENCIL_UNIFORM
//#####################################################################
#ifndef __BASIS_STENCIL_UNIFORM__
#define __BASIS_STENCIL_UNIFORM__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Symbolics/MULTIVARIATE_POLYNOMIAL.h>
namespace PhysBAM{

template<class TV>
struct BASIS_STENCIL_UNIFORM
{
    typedef VECTOR<int,TV::m> TV_INT;
    TV_INT center_offset; // In units of dx/2

    struct ENTRY
    {
        RANGE<TV_INT> region; // In units of dx/2
        MULTIVARIATE_POLYNOMIAL<TV> polynomial; // In units of dx
    };

    ARRAY<ENTRY> stencils;

    struct DICED
    {
        TV_INT index_offset;
        RANGE<TV_INT> range; // Subset of [-1,1)
        MULTIVARIATE_POLYNOMIAL<TV> polynomial;
    };

    ARRAY<DICED> diced;

    BASIS_STENCIL_UNIFORM();
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
};
}
#endif
