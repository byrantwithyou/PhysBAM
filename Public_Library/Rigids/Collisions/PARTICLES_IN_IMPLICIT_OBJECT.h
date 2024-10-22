//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace PARTICLES_IN_IMPLICIT_OBJECT
//##################################################################### 
#ifndef __PARTICLES_IN_IMPLICIT_OBJECT__
#define __PARTICLES_IN_IMPLICIT_OBJECT__

#include <Core/Matrices/MATRIX_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> struct RIGID_BODY_PARTICLE_INTERSECTION;
template<class T,class ID> class ARRAY_VIEW;
template<class T,class ID> class ARRAY;
template<class TV> class IMPLICIT_OBJECT;
template<class T,int d> class VECTOR;
template<class K,class T> class HASHTABLE;
template<class TV> class GEOMETRY_PARTICLES;
template<class TV,class T_ARRAY> class PARTICLE_HIERARCHY;

template<class TV>
struct PARTICLES_IN_IMPLICIT_OBJECT
{
    static void Append_All_Intersections(RIGID_BODY<TV>& body0, RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy,const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test);
    static void Append_All_Intersections_Points(RIGID_BODY<TV>& body0, RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value);
    static void Append_All_Intersections_Triangles(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy_center_phi_test);
    static void Append_All_Intersections_Edges(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy_center_phi_test);
    static void Get_Interfering_Simplices(const RIGID_BODY<TV>& body0,const RIGID_BODY<TV>& body1,ARRAY<int>& simplex_list,MATRIX<typename TV::SCALAR,TV::m,TV::m>& rotation,TV& translation,const bool use_triangle_hierarchy_center_phi_test);
    static const typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::HIERARCHY& Simplex_Hierarchy(const RIGID_BODY<TV>& rigid_body);
    static void Intersections_Using_Hierarchy(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<int>& simplex_list,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early,MATRIX<typename TV::SCALAR,TV::m,TV::m>& rotation,TV& translation);
    static void Intersections_Using_Hierarchy_And_Edges(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,ARRAY<int>& simplex_list1,ARRAY<int>& simplex_list2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early,MATRIX<typename TV::SCALAR,TV::m,TV::m>& rotation,TV& translation);
    static void Particles_In_Implicit_Object(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early);
    static void Particles_In_Implicit_Object_Hierarchy(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,HASHTABLE<const GEOMETRY_PARTICLES<TV>*,const PARTICLE_HIERARCHY<TV>*>& particle_hierarchies);
    static void Particles_In_Implicit_Object_Partition(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_particle_partition_center_phi_test,const VECTOR<int,TV::m>& particle_partition_size,const bool exit_early);
};
}
#endif
