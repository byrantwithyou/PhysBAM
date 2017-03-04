//#####################################################################
// Copyright 2006, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_MACCORMACK_UNIFORM
//#####################################################################
#ifndef __ADVECTION_MACCORMACK_UNIFORM__
#define __ADVECTION_MACCORMACK_UNIFORM__

#include <Grid_PDE/Advection/ADVECTION.h>
namespace PhysBAM{

template<class TV,class T2,class T_NESTED_ADVECTION>
class ADVECTION_MACCORMACK_UNIFORM:public ADVECTION<TV,T2>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
private:
    // false to use 1st order advection
    const ARRAY<bool,TV_INT>* node_mask;
    const ARRAY<bool,TV_INT>* cell_mask;
    const ARRAY<bool,FACE_INDEX<TV::m> >* face_mask;
    T_NESTED_ADVECTION& nested_advection;
public:
    bool clamp_extrema; // otherwise fall back to 1st order
    bool ensure_second_order; // sacrifice stability for accuracy.

    ADVECTION_MACCORMACK_UNIFORM(T_NESTED_ADVECTION& nested_advection_input,const ARRAY<bool,TV_INT>* node_mask_input,const ARRAY<bool,TV_INT>* cell_mask_input,const ARRAY<bool,FACE_INDEX<TV::m> >* face_mask_input);

    void Update_Advection_Equation_Node(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost_input=0,const ARRAY<T2,TV_INT>* Z_max_ghost_input=0,ARRAY<T2,TV_INT>* Z_min_input=0,ARRAY<T2,TV_INT>* Z_max_input=0);
    void Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,const FACE_LOOKUP_UNIFORM<TV>& face_velocities,
        BOUNDARY<TV,T2>& boundary,const T dt,const T time,const ARRAY<T2,TV_INT>* Z_min_ghost_input,const ARRAY<T2,TV_INT>* Z_max_ghost_input,ARRAY<T2,TV_INT>* Z_min_input,
        ARRAY<T2,TV_INT>* Z_max_input);
    void Update_Advection_Equation_Face_Lookup(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& Z,const FACE_LOOKUP_UNIFORM<TV>& Z_ghost,
        const FACE_LOOKUP_UNIFORM<TV>& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,const FACE_LOOKUP_UNIFORM<TV>* Z_min_ghost_input,
        const FACE_LOOKUP_UNIFORM<TV>* Z_max_ghost_input,ARRAY<T,FACE_INDEX<TV::m> >* Z_min_input,ARRAY<T,FACE_INDEX<TV::m> >* Z_max_input);
//#####################################################################
};
}
#endif
