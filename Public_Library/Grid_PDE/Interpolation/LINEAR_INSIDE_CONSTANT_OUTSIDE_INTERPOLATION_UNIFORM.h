//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Eran Guendelman, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM
//#####################################################################
#ifndef __LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM__
#define __LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM__

#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV>
class LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM:public LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    template<class T3> struct REBIND{typedef LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<TV,T3,T_FACE_LOOKUP> TYPE;};
    using LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_To_Array;

    LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM()
    {}

    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAY<T2,TV_INT>& u,const TV& X) const
    {RANGE<TV_INT> domain=u.Domain_Indices();TV X_new(clamp(X,grid.X(domain.min_corner),grid.X(domain.max_corner)));
    return LINEAR_INTERPOLATION_UNIFORM<TV,T2>::Clamped_To_Array(grid,u,X_new);}

//#####################################################################
};
}
#endif
