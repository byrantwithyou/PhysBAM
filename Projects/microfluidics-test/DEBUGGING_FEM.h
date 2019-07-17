//#####################################################################
// Copyright 2019.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __DEBUGGING_FEM__
#define __DEBUGGING_FEM__
#include "BLOCK_MATRIX.h"
#include "BLOCK_VECTOR.h"
#include "COMPONENT_LAYOUT_FEM.h"
#include <string>

namespace PhysBAM{

template<class T> struct CACHED_ELIMINATION_MATRIX;
template<class TV> struct MATRIX_CONSTRUCTION_FEM;

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
    void Visualize_Tetrahedron(BLOCK_ID b) const;
    void Visualize_Tetrahedron_Dofs(const MATRIX_CONSTRUCTION_FEM<TV3>& mc) const;
    void Visualize_Tetrahedron_Dofs(const MATRIX_CONSTRUCTION_FEM<TV2>& mc) const {}
    template<int d>
    void Highlight_Dof(BLOCK_ID b,int vep,int r,int dim) const;
    void Visualize_Domain(const std::string& name,bool fill,const IV2& dim,const RANGE<TV2>& anno,const RANGE<TV2>& bound) const;
    void Visualize_Meshing(const std::string& name,const IV2& dim,const RANGE<TV2>& domain) const;
};

}
#endif
