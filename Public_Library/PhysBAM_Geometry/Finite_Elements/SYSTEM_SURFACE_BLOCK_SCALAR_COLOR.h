//#####################################################################
// Copyright 2012, Craig Schroeder, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYSTEM_SURFACE_BLOCK_SCALAR_COLOR 
//#####################################################################

#ifndef __SYSTEM_SURFACE_BLOCK_SCALAR_COLOR__
#define __SYSTEM_SURFACE_BLOCK_SCALAR_COLOR__

#include <PhysBAM_Tools/Symbolics/STATIC_POLYNOMIAL.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Finite_Elements/SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR.h>
#include <PhysBAM_Geometry/Finite_Elements/BOUNDARY_CONDITIONS_SCALAR_COLOR.h>

namespace PhysBAM{

template<class TV> class CELL_DOMAIN_INTERFACE_COLOR;

template<class TV,int static_degree>
class SYSTEM_SURFACE_BLOCK_SCALAR_COLOR:public NONCOPYABLE
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
    ARRAY<VECTOR_ND<T> >* f_surface;
    BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>* abc;
    ARRAY<OVERLAP_POLYNOMIAL> overlap_polynomials;

    template<int d>
    void Initialize(SYSTEM_SURFACE_BLOCK_SCALAR_HELPER_COLOR<TV>& helper_input,const BASIS_STENCIL_UNIFORM<TV,d>& s,
        BOUNDARY_CONDITIONS_SCALAR_COLOR<TV>* abc_input,ARRAY<VECTOR_ND<T> >& f_surface_input,T scale_input);

    void Add_Entry(int constraint_index,int orientation,int flat_index_diff_ref,int color,T value)
    {helper->data(color)(constraint_index,flat_index_diff_ref)+=value*scale;}

    int Flat_Diff(int i)
    {return helper->flat_diff(i);}

    void Resize()
    {helper->Resize();}
};
}
#endif
