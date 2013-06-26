//#####################################################################
// Copyright 2009, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS
//#####################################################################
#ifndef __INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS__
#define __INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS__
#include <Tools/Grids_Uniform/FACE_INDEX.h>
namespace PhysBAM{
template<class TV> class GRID;

template<class TV>
class INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS
{
    typedef typename TV::SCALAR T;typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;

public:
    INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS()
    {}
    virtual ~INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS()
    {}

//#####################################################################
    virtual void Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const GRID<TV>& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U,
        const ARRAY<bool,TV_INT>& euler_psi,const T_ARRAYS_SCALAR& p_cell,T_FACE_ARRAYS_SCALAR& p_face) const=0;
    virtual void Fill_Incompressible_Beta_Face(const GRID<TV>& grid,T_FACE_ARRAYS_SCALAR& beta_face) const=0;
//#####################################################################
};  
}   
#endif

