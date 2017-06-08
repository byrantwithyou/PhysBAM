//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO
//#####################################################################
#ifndef __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO__
#define __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO__
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAYS_FORWARD.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>

namespace PhysBAM{

template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;
template<class TV> class GRID;
template<class TV> class COLLISION_GEOMETRY;
template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;
template<class TV> class IMPLICIT_BOUNDARY_CONDITION_COLLECTION;

template<class TV>
struct COLLISION_FACE_INFO
{
    FACE_INDEX<TV::m> face;
    int side;
    ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> > simplices;

    bool operator<(const COLLISION_FACE_INFO& cfi) const
    {
        if(face.axis!=cfi.face.axis)
            return face.axis<cfi.face.axis;
        for(int i=0;i<TV::m;i++)
            if(face.index(i)!=cfi.face.index(i))
                return face.index(i)<cfi.face.index(i);
        return false;
    }
};

template<class TV>
class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO
{
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;

    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collision_bodies_affecting_fluid;

public:
    ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID> coupling_bodies;
    ARRAY<COLLISION_FACE_INFO<TV> > collision_face_info;
    GRID<TV>& grid;
    const ARRAY<bool,TV_INT>* outside_fluid;
    const bool& use_collision_face_neighbors;
    T iterator_rasterization_thickness;

    explicit UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& collection);
    UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO&) = delete;
    void operator=(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO&) = delete;

//#####################################################################
    void Initialize_Collision_Aware_Face_Iterator(const ARRAY<bool,TV_INT>& outside_fluid_input,const ARRAY<bool,FACE_INDEX<d> >& kinematic_faces,int ghost_cells,const bool disable_thinshell);
    void Register_Neighbors_As_Collision_Faces();
//#####################################################################
};
}
#endif
