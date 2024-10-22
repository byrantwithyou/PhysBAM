//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE_INTERACTION_GEOMETRY
//#####################################################################
#ifndef __STRUCTURE_INTERACTION_GEOMETRY__
#define __STRUCTURE_INTERACTION_GEOMETRY__

#include <Core/Arrays/INDIRECT_ARRAY.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/choice.h>
#include <Tools/Parallel_Computation/PARTITION_ID.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV>
class STRUCTURE_INTERACTION_GEOMETRY
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    struct UNUSABLE{};
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT T_SEGMENTED_CURVE;
public:
    GEOMETRY_PARTICLES<TV>& full_particles;
    ARRAY<int> active_indices;
    TRIANGULATED_SURFACE<T>* triangulated_surface;
    T_SEGMENTED_CURVE* segmented_curve;
    POINT_SIMPLICES_1D<T>* point_simplices;
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > subset;
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > particle_hierarchy;
    ARRAY<char> triangulated_surface_processor_masks,segmented_curve_processor_masks,point_processor_masks;
    ARRAY<bool> triangulated_surface_modified,segmented_curve_modified,point_modified; // whether boxes below you have been modified
private:
    bool need_destroy_segmented_curve;
    UNUSABLE unusable;
public:

    STRUCTURE_INTERACTION_GEOMETRY(GEOMETRY_PARTICLES<TV>& full_particles_input)
        :full_particles(full_particles_input),triangulated_surface(0),segmented_curve(0),point_simplices(0),
        subset(full_particles.X,active_indices),particle_hierarchy(subset,false,0),need_destroy_segmented_curve(false)
    {}

    ~STRUCTURE_INTERACTION_GEOMETRY()
    {Clean_Memory();}

    void Clean_Memory()
    {if(need_destroy_segmented_curve) delete segmented_curve;need_destroy_segmented_curve=false;
    active_indices.Clean_Memory();triangulated_surface=0;segmented_curve=0;point_simplices=0;particle_hierarchy.Clean_Memory();}

    void Build_Topological_Structure_Of_Hierarchies()
{
    if(triangulated_surface) triangulated_surface->Initialize_Hierarchy(false);
    if(segmented_curve) segmented_curve->Initialize_Hierarchy(false);
    ARRAY_VIEW<TV> tmp_view(full_particles.X.Get_Array_Pointer(),full_particles.X.Size());
    subset.array.Exchange(tmp_view);
    particle_hierarchy.Initialize_Hierarchy_Using_KD_Tree();}

    const typename conditional<d==2,ARRAY<int>,typename conditional<d==3,ARRAY<VECTOR<int,2> >,UNUSABLE>::type>::type& Edges() const
    {return choice<d-1>(unusable,active_indices,segmented_curve->mesh.elements);}

    const typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d-1>::OBJECT* Face_Mesh_Object() const
    {return choice<d-1>(point_simplices,segmented_curve,triangulated_surface);}

    const ARRAY<bool>& Edge_Modified() const
    {return d==2?point_modified:segmented_curve_modified;}

    const ARRAY<bool>& Face_Modified() const
    {return d==2?segmented_curve_modified:triangulated_surface_modified;}

    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d-1>::HIERARCHY& Face_Hierarchy() const
    {return *Face_Mesh_Object()->hierarchy;}

    bool Has_Edges() const
    {return d==2 || segmented_curve!=0;}

    const typename conditional<(d<2),UNUSABLE,typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d<2?0:d-2>::HIERARCHY>::type& Edge_Hierarchy() const
    {return choice<d-1>(unusable,particle_hierarchy,*segmented_curve->hierarchy);}

    const ARRAY<char>& Edge_Processor_Masks() const
    {return d==2?point_processor_masks:segmented_curve_processor_masks;}

    const ARRAY<char>& Face_Processor_Masks() const
    {return d==2?segmented_curve_processor_masks:triangulated_surface_processor_masks;}

//#####################################################################
    static TRIANGULATED_SURFACE<T>* Triangulated_Surface(STRUCTURE<TV>* structure);
    static T_SEGMENTED_CURVE* Segmented_Curve(STRUCTURE<TV>* structure);
    bool Build_Collision_Geometry(STRUCTURE<TV>& structure);
    void Update_Faces_And_Hierarchies_With_Collision_Free_Positions(ARRAY_VIEW<const T> node_thickness,const T node_thickness_multiplier,ARRAY_VIEW<const TV> X_old_full);
    void Update_Processor_Masks(const PARTITION_ID processor,const ARRAY<PARTITION_ID>& partition_id_from_particle_index);
private:
    template<class T_OBJECT,class T_HIERARCHY>
    void Update_Processor_Masks_Helper(T_OBJECT& object,T_HIERARCHY& hierarchy,const PARTITION_ID processor,const ARRAY<PARTITION_ID>& partition_id_from_particle_index,ARRAY<char>& mask);
//#####################################################################
};
}
#endif
