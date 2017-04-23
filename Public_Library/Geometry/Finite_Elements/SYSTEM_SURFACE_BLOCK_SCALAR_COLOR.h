//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_SURFACE_BLOCK_SCALAR_COLOR 
//#####################################################################

#ifndef __SYSTEM_SURFACE_BLOCK_SCALAR_COLOR__
#define __SYSTEM_SURFACE_BLOCK_SCALAR_COLOR__

#include <Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>
#include <functional>

namespace PhysBAM{

template<class TV> class CELL_DOMAIN_INTERFACE_COLOR;

template<class TV,int static_degree>
class SYSTEM_SURFACE_BLOCK_SCALAR_COLOR
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

    SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV> *helper;

public:

    struct OVERLAP_POLYNOMIAL
    {
        int flat_index_offset;
        int flat_index_diff_ref;
        int subcell; // flags indicating fine cells
        STATIC_POLYNOMIAL<T,TV::m,static_degree> polynomial;
    };
    
    T scale;
    ARRAY<ARRAY<T> >* rhs;
    bool use_discontinuous_scalar_field;
    std::function<T(const TV& X,int color0,int color1)> u_jump;
    std::function<T(const TV& X,int color0,int color1)> j_surface;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;

    SYSTEM_SURFACE_BLOCK_SCALAR_COLOR() = default;
    SYSTEM_SURFACE_BLOCK_SCALAR_COLOR(const SYSTEM_SURFACE_BLOCK_SCALAR_COLOR&) = delete;
    void operator=(const SYSTEM_SURFACE_BLOCK_SCALAR_COLOR&) = delete;

    template<int d>
    void Initialize(SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,
        bool use_discontinuous_scalar_field_input,std::function<T(const TV& X,int color0,int color1)> u_jump_input,
        std::function<T(const TV& X,int color0,int color1)> j_surface_input,ARRAY<ARRAY<T> >& rhs_input,T scale_input);

    void Add_Entry(int constraint_index,int flat_index_diff_ref,int color,T value)
    {helper->data(color)(constraint_index,flat_index_diff_ref)+=value*scale;}

    void Add_Constraint_Rhs_Entry(int constraint_index,int color,T value)
    {helper->rhs_data(color)(constraint_index)+=value*scale;}

    int Flat_Diff(int i)
    {return helper->flat_diff(i);}

    void Resize()
    {helper->Resize();}
};
}
#endif
