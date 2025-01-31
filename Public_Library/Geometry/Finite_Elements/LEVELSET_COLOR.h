//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_COLOR
//#####################################################################
#ifndef __LEVELSET_COLOR__
#define __LEVELSET_COLOR__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Grid_Tools/Grids/GRID.h>

namespace PhysBAM{

template<class TV>
class LEVELSET_COLOR
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;

public:
    GRID<TV>& grid;
    ARRAY<T,TV_INT>& phi;
    ARRAY<int,TV_INT>& color;

    LEVELSET_COLOR(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,ARRAY<int,TV_INT>& color_input)
        :grid(grid_input),phi(phi_input),color(color_input)
    {}
    LEVELSET_COLOR(const LEVELSET_COLOR&) = delete;
    void operator=(const LEVELSET_COLOR&) = delete;

    T Phi(const TV_INT& index) const
    {return phi(index);}

    int Color(const TV_INT& index) const
    {return color(index);}

    T Phi_And_Color(const TV_INT& index,int& c) const
    {c=color(index);return phi(index);}

    T Phi(const TV& X) const
    {int c;return Phi_And_Color(X,c);}

    int Color(const TV& X) const
    {int c;Phi_And_Color(X,c);return c;}

    T Phi_And_Color(const TV& X,int& c) const;

    void Get_Raw_Levelset_For_Color(ARRAY<T,TV_INT>& color_phi,int c,int ghost) const;
    void Get_Levelset_For_Color(ARRAY<T,TV_INT>& color_phi,int c,int ghost) const;
};
}
#endif
