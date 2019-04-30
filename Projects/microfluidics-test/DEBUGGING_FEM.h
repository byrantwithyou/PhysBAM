//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DEBUGGING_FEM__
#define __DEBUGGING_FEM__
#include "BLOCK_MATRIX.h"
#include "BLOCK_VECTOR.h"
#include "COMPONENT_LAYOUT_FEM.h"

namespace PhysBAM{

template<class T> struct CACHED_ELIMINATION_MATRIX;

template<class T>
struct DEBUGGING_FEM
{
    typedef VECTOR<T,2> TV2;
    typedef VECTOR<int,2> IV2;
    typedef VECTOR<T,3> TV3;
    COMPONENT_LAYOUT_FEM<T>& cl;

    DEBUGGING_FEM(COMPONENT_LAYOUT_FEM<T>& cl);
    
    void Visualize_Block_State(BLOCK_ID b) const;
    void Visualize_Block_Dofs(BLOCK_ID b) const;
    void Visualize_Solution(const BLOCK_VECTOR<TV2>& U,BLOCK_ID b,bool remap_dofs) const;
    void Visualize_Solution(const BLOCK_VECTOR<TV3>& U,BLOCK_ID b,bool remap_dofs) const;
    void Visualize_Flat_Dofs() const;
    void Visualize_Ticks(BLOCK_ID b,bool reference_ticks) const;
    void Check_Analytic_Solution() const;
    void Dump_World_Space_System() const;
    void Dump_World_Space_Vector(const char* name) const;
    void Visualize_Tetrahedron(BLOCK_ID b) const;
};

}
#endif
