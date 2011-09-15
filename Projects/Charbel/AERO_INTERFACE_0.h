//#####################################################################
// Copyright 2007, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AERO_INTERFACE_0
//##################################################################### 
#ifndef __AERO_INTERFACE_0__
#define __AERO_INTERFACE_0__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{

template<class T>
class AERO_INTERFACE_0{
  typedef GRID<TV> T_GRID;
  typedef typename T_GRID::VECTOR_T TV;
 private:
    bool is_initialized;
    T_GRID grid;
    LEVELSET_UNIFORM<T_GRID > *levelset;

public:

    AERO_INTERFACE_0(T_GRID domain)
        :is_initialized(false),grid(domain)
    {}

    ~AERO_INTERFACE_0()
    {delete levelset;}
    
    const T Phi(const TV& position) const
    {assert(is_initialized);return levelset->Phi(position);}

    const TV Gradient(const TV& position) const
    {assert(is_initialized);TV gradient;levelset->Compute_Gradient(&gradient);return gradient;}

    void Compute_Level_Set(TRIANGULATED_SURFACE<T>& triangulated_surface)
    {LEVELSET_MAKER_UNIFORM<T> levelset_maker;
    ARRAY<T,VECTOR<int,3> > phi;
    is_initialized=true;
    levelset_maker.Compute_Level_Set(triangulated_surface,grid,phi);
    levelset=new LEVELSET_UNIFORM<T_GRID>(grid,phi);}
};
}
#endif
