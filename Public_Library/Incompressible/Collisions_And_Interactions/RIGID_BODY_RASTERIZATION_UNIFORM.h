//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#ifndef __RIGID_BODY_RASTERIZATION_UNIFORM__
#define __RIGID_BODY_RASTERIZATION_UNIFORM__
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>

namespace PhysBAM{
template<class TV> class COLLISION_GEOMETRY;
template<class TV, class ID> class OBJECTS_IN_CELL;
template<class TV> class RANGE;

namespace RASTERIZATION{
//#####################################################################
template<class TV> void Rasterize_Object(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID& id);
template<class T,class TV> void Compute_Occupied_Blocks(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,
    ARRAY<bool,VECTOR<int,TV::m> >& occupied,const bool with_body_motion,const T& extra_thickness,const T& body_thickness_factor);
template<class TV> void Rasterize_Box_Onto_Blocks(const GRID<TV>& grid,ARRAY<bool,VECTOR<int,TV::m> >& occupied,const RANGE<TV>& box);
template<class TV> void Rasterize_Box(const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects_in_cell,const RANGE<TV>& box,const COLLISION_GEOMETRY_ID id);
//#####################################################################
};
};
#endif
