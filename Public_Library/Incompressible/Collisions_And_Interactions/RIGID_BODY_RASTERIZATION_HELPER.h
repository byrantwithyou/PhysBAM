//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace RASTERIZATION
//##################################################################### 
#ifndef __RIGID_BODY_RASTERIZATION_HELPER__
#define __RIGID_BODY_RASTERIZATION_HELPER__
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/RANGE.h>
#include <Tools/Vectors/VECTOR.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Incompressible/Collisions_And_Interactions/RIGID_BODY_RASTERIZATION_UNIFORM.h>
namespace PhysBAM{
template<class TV> struct GRID_ARRAYS_POLICY;
template<class TV,class COLLISION_GEOMETRY_ID> class OBJECTS_IN_CELL;

namespace RASTERIZATION{
//#####################################################################
// Function Prepare_For_World_Space_Simplex_Bounding_Box
//#####################################################################
template<class T> void Prepare_For_World_Space_Simplex_Bounding_Box(const RIGID_BODY<VECTOR<T,1> >& rigid_body)
{}
//#####################################################################
// Function Prepare_For_World_Space_Simplex_Bounding_Box
//#####################################################################
template<class T> void Prepare_For_World_Space_Simplex_Bounding_Box(const RIGID_BODY<VECTOR<T,2> >& rigid_body)
{
    if(!rigid_body.simplicial_object->segment_list) const_cast<RIGID_BODY<VECTOR<T,2> >&>(rigid_body).simplicial_object->Update_Segment_List();
}
//#####################################################################
// Function Prepare_For_World_Space_Simplex_Bounding_Box
//#####################################################################
template<class T> void Prepare_For_World_Space_Simplex_Bounding_Box(const RIGID_BODY<VECTOR<T,3> >& rigid_body)
{
    if(!rigid_body.simplicial_object->triangle_list) const_cast<RIGID_BODY<VECTOR<T,3> >&>(rigid_body).simplicial_object->Update_Triangle_List();
}
//#####################################################################
// Function Rasterize_Object_Generic
//#####################################################################
template<class TV> void Rasterize_Object_Generic(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>& objects,
    const COLLISION_GEOMETRY_ID id)
{
    const RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry;
    rigid_collision_geometry=dynamic_cast<const RIGID_COLLISION_GEOMETRY<TV>*>(&collision_geometry);
    if(rigid_collision_geometry) Prepare_For_World_Space_Simplex_Bounding_Box(rigid_collision_geometry->rigid_body);
    if(collision_geometry.Number_Of_Simplices())
        for(int t=0;t<collision_geometry.Number_Of_Simplices();t++){
            typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;
            T_SIMPLEX simplex=collision_geometry.World_Space_Simplex(t);
            VECTOR<TV,TV::dimension> pts;
            for(int d=0;d<TV::dimension;d++) pts(d)=simplex.X(d);
            RANGE<TV> box=RANGE<TV>::Bounding_Box(pts).Thickened(collision_geometry.collision_thickness);
            Rasterize_Box(grid,objects,box,id);}
    else
        for(COLLISION_GEOMETRY_ID i(0);i<collision_geometry.collision_geometries_for_rasterization->m;i++)
            if((*collision_geometry.collision_geometries_for_rasterization)(i) && (*collision_geometry.collision_geometries_for_rasterization)(i)->active)
                Rasterize_Object(*(*collision_geometry.collision_geometries_for_rasterization)(i),grid,objects,id);
}
//#####################################################################
// Macro SPECIALIZE_RASTERIZE_OBJECT
//#####################################################################
#define SPECIALIZE_RASTERIZE_OBJECT(TV) \
    void Rasterize_Object_Generic(const COLLISION_GEOMETRY<TV>&,const GRID<TV>&,OBJECTS_IN_CELL<TV,COLLISION_GEOMETRY_ID>&,const COLLISION_GEOMETRY_ID&) \
    {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Compute_Occupied_Blocks_Generic
//#####################################################################
// NOTE: Since Compute_Occupied_Blocks is a virtual function we can't directly make it a template member
template<class T,class TV> void
Compute_Occupied_Blocks_Generic(const COLLISION_GEOMETRY<TV>& collision_geometry,const GRID<TV>& grid,ARRAY<bool,VECTOR<int,TV::m> >& occupied,
    const bool with_body_motion,const T extra_thickness,const T body_thickness_factor)
{
    const RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry;
    rigid_collision_geometry=dynamic_cast<const RIGID_COLLISION_GEOMETRY<TV>*>(&collision_geometry);
    if(rigid_collision_geometry) Prepare_For_World_Space_Simplex_Bounding_Box(rigid_collision_geometry->rigid_body);
    for(int t=0;t<collision_geometry.Number_Of_Simplices();t++){
        typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;
        T_SIMPLEX simplex=collision_geometry.World_Space_Simplex(t);
        VECTOR<TV,TV::dimension> pts;
        for(int d=0;d<TV::dimension;d++) pts(d)=simplex.X(d);
        RANGE<TV> box=RANGE<TV>::Bounding_Box(pts);
        if(with_body_motion){
            T_SIMPLEX simplex_saved=collision_geometry.World_Space_Simplex(t,1);
            VECTOR<TV,TV::dimension> pts_saved;
            for(int d=0;d<TV::dimension;d++) pts_saved(d)=simplex_saved.X(d);
            box.Enlarge_To_Include_Box(RANGE<TV>::Bounding_Box(pts_saved));}
        box.Change_Size(extra_thickness+body_thickness_factor*collision_geometry.collision_thickness);
        Rasterize_Box_Onto_Blocks(grid,occupied,box);}
}
//#####################################################################
// Macro SPECIALIZE_COMPUTE_OCCUPIED
//#####################################################################
#define SPECIALIZE_COMPUTE_OCCUPIED(TV) \
    void Compute_Occupied_Blocks_Generic(const COLLISION_GEOMETRY<TV>&,const GRID<TV>&,ARRAY<bool>&,const bool,const TV::SCALAR,const TV::SCALAR) \
    {PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
};
};
#endif
