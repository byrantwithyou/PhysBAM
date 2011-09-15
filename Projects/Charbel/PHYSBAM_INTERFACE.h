//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Interface for AERO-F
//##################################################################### 
#ifndef __PHYSBAM_INTERFACE__
#define __PHYSBAM_INTERFACE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_HIERARCHY;
template<class TV> class GEOMETRY_PARTICLES;
class TRIANGLE_MESH;

template<class T>
struct IntersectionResult {
  int triangleID; // -> -1 if no intersection
  T alpha; // Intersection is at alpha*edgeNode1 + (1-alpha)*edgeNode2
  T zeta[3];  // Intersection is at zeta[0]*triNode1+zeta[1]*triNode2+zeta[2]*triNode3
};

template<class T>
class PhysBAMInterface {
    typedef VECTOR<T,3> TV;

public:
    TRIANGLE_MESH& triangle_mesh;
    GEOMETRY_PARTICLES<TV>& particles;
    TRIANGLE_HIERARCHY<T>* triangle_hierarchy;

    TRIANGULATED_SURFACE<T>* triangulated_surface; // Deprecated

    PhysBAMInterface(TRIANGULATED_SURFACE<T> &triangulated_surface); // Deprecated
    PhysBAMInterface(TRIANGLE_MESH& triangle_mesh,GEOMETRY_PARTICLES<TV>& particles);
    ~PhysBAMInterface();

    // When triangulated_surface.particles have moved, then call the following
    // method to update the triangle hierarchy
    // pass rebuild_hierarchy=true when topology changes
    void Update(const bool rebuild_hierarchy=false);

    int GetClosestTriangle(const TV position,const TV min_corner,const TV max_corner,T* distance=0) const;

    // edges_and_results(i).x is prepopulated with edge information. edges_and_results(i).y stores the result of the intersection.
    void Intersect(const ARRAY<TV>& node_positions,ARRAY<PAIR<VECTOR<int,2>,IntersectionResult<T> > >& edges_and_results,const T thickness=(T)1e-4) const;
    IntersectionResult<T> Intersect(const TV& start,const TV& end,const T thickness) const;
    void Intersect(const ARRAY<TV>& node_positions,ARRAY<bool>& occluded_node,ARRAY<TRIPLE<VECTOR<int,3>,IntersectionResult<T>,IntersectionResult<T> > >& edges_and_results,const T thickness) const;
    bool HasCloseTriangle(const TV min_corner,const TV max_corner,bool* is_occluded) const;
    bool HasCloseTriangle(const TV position,const TV min_corner,const TV max_corner,const T thickness,bool* is_occluded) const;
    void computeSweptNodes(const ARRAY<TV>& node_positions,ARRAY<bool>& swept_node,const T dt,const T thickness);
};

}
#endif
