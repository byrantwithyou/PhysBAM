//#####################################################################
// Copyright 2003-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_PARTITION
//#####################################################################
#ifndef __PARTICLE_PARTITION__
#define __PARTICLE_PARTITION__

#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_PARTITION
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<int,TV::m> TV_INT;typedef ARRAY<T,TV_INT> T_ARRAYS_T;
public:
    GRID<TV> grid;
    ARRAY<ARRAY<int>,TV_INT> partition;
private:
    bool use_radius;
    T_ARRAYS_T radius;
public:

    PARTICLE_PARTITION(const RANGE<TV>& box,const TV_INT& counts,const GEOMETRY_PARTICLES<TV>& particles,
        const bool use_radius=true,const bool is_mac_grid=true);
    void Add_To_Partition(const TV& location,const int particle_id);
    RANGE<TV_INT> Range(const RANGE<TV>& box) const;
    void Intersection_List(const IMPLICIT_OBJECT<TV>& test_surface,const MATRIX<T,d>& rotation,
        const TV& translation,ARRAY<TV_INT>& intersection_list,const T contour_value=0) const;
    void Proximity_List(const TV& location,const T proximity,ARRAY<int>& proximity_list);
//#####################################################################
};
}
#endif
