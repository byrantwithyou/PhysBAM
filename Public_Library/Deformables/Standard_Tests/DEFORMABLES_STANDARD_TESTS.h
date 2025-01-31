//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_STANDARD_TESTS
//#####################################################################
#ifndef __DEFORMABLES_STANDARD_TESTS__
#define __DEFORMABLES_STANDARD_TESTS__

#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/OPENSUBDIV_SURFACE.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>
#include <functional>
namespace PhysBAM{

template<class TV> class SOFT_BINDINGS;
template<class TV> class DEFORMABLE_PARTICLES;
template<class TV> class GEOMETRY_PARTICLES;
template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T> class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE;
template<class TV> class BINDING_LIST;
template<class TV,class T_PARTICLES> class PARTICLES_SUBSET;
template<class TV,int d> class EMBEDDED_MATERIAL_SURFACE;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class GRID;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
class VIEWER_DIR;

template<class TV>
class DEFORMABLES_STANDARD_TESTS
{
protected:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE T_SEGMENTED_CURVE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
public:
    STREAM_TYPE stream_type;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;

    DEFORMABLES_STANDARD_TESTS(STREAM_TYPE stream_type,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input);
    virtual ~DEFORMABLES_STANDARD_TESTS();

    TRIANGULATED_AREA<T>& Create_Mattress(const GRID<TV>& mattress_grid,const bool use_constant_mass,const RIGID_BODY_STATE<TV>& initial_state)
    {return Create_Mattress(mattress_grid,use_constant_mass,&initial_state,1000);}

    template<class T_OBJECT> void
    Substitute_Soft_Bindings_For_Embedded_Nodes(T_OBJECT& object,SOFT_BINDINGS<TV>& soft_bindings,HASHTABLE<int,int>* persistent_soft_bindings=0)
    {Substitute_Soft_Bindings_For_Nodes(object,soft_bindings,persistent_soft_bindings,true);}

    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_BODY_STATE<TV>& initial_state,
        std::function<void (TRIANGULATED_SURFACE<T>&)> clipping_function,ARRAY<int>* particle_indices=0)
    {return Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,&initial_state,clipping_function,particle_indices);}

    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_BODY_STATE<TV>* initial_state,ARRAY<int>* particle_indices=0)
    {return Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,initial_state,0,particle_indices);}

    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_BODY_STATE<TV>& initial_state,ARRAY<int>* particle_indices=0)
    {return Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,&initial_state,particle_indices);}

//#####################################################################
    virtual void Add_Gravity();
    template<class T_STRUCTURE>
    T_STRUCTURE& Copy_And_Add_Structure(T_STRUCTURE& structure,ARRAY<int>* particle_indices=0,bool delete_structure=true);
    void Set_Initial_Particle_Configuration(GEOMETRY_PARTICLES<TV>& particles,const RIGID_BODY_STATE<TV>& state,const bool relative_to_box_center);
    T_SEGMENTED_CURVE& Create_Segmented_Curve(const GRID<VECTOR<T,1> >& square_grid,const RIGID_BODY_STATE<TV>& initial_state,const T density);
    T_SEGMENTED_CURVE& Create_Segmented_Curve(const int m,const RIGID_BODY_STATE<TV>& initial_state,const T initial_radius,const T density);
    TETRAHEDRALIZED_VOLUME<T>& Create_Tetrahedralized_Volume(const std::string& filename,const RIGID_BODY_STATE<TV>& initial_state,const bool relative_to_box_center,
        const bool use_constant_mass,const T density,const T scale=1);
    T_TRIANGULATED_OBJECT& Create_Triangulated_Object(const GRID<TV>& square_grid,const RIGID_BODY_STATE<TV>& initial_state,const T density);
    T_TRIANGULATED_OBJECT& Create_Triangulated_Object(const std::string& filename,const RIGID_BODY_STATE<TV>& initial_state,const bool relative_to_box_center,
        const bool use_constant_mass,const T scale=1);
    T_SEGMENTED_CURVE& Create_Segmented_Curve(const std::string& filename,const RIGID_BODY_STATE<TV>& initial_state,const bool relative_to_box_center,
        const bool use_constant_mass);
    TRIANGULATED_AREA<T>& Create_Mattress(const GRID<VECTOR<T,2> >& mattress_grid,const bool use_constant_mass,const RIGID_BODY_STATE<TV>* initial_state,const T density,const bool reverse_triangles=false);
    TETRAHEDRALIZED_VOLUME<T>& Create_Mattress(const GRID<VECTOR<T,3> >& mattress_grid,const bool use_constant_mass,const RIGID_BODY_STATE<TV>* initial_state,const T density);
    TETRAHEDRALIZED_VOLUME<T>& Create_Cylinder(const CYLINDER<T>& cylinder,int num_elements_height,int num_elements_radius,const bool use_constant_mass,const RIGID_BODY_STATE<TV>* initial_state,const T density);
    template<class T_SHAPE>
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& Create_Embedded_Tetrahedralized_Volume(const T_SHAPE& shape,const RIGID_BODY_STATE<TV>& initial_state,const bool relative_to_box_center);
    template<class T_OBJECT> void
    Substitute_Soft_Bindings_For_Nodes(T_OBJECT& object,SOFT_BINDINGS<TV>& soft_bindings,HASHTABLE<int,int>* persistent_soft_bindings=0,const bool embedded_only=false,
        const bool use_impulses_for_collisions=true);
    LEVELSET_IMPLICIT_OBJECT<TV>* Read_Or_Initialize_Implicit_Surface(const std::string& levelset_filename,const VIEWER_DIR& viewer_dir,TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface) const;
    void Initialize_Tetrahedron_Collisions(const int id_number,VIEWER_DIR& viewer_dir,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,
        TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters,TRIANGULATED_SURFACE<T>* triangulated_surface=0);
    TRIANGULATED_SURFACE<T>& Create_Drifted_Surface(const TRIANGULATED_SURFACE<T>& triangulated_surface,SOFT_BINDINGS<TV>& soft_bindings,const bool use_impulses_for_collisions=false) const;
    template <class T_OBJECT> static void Set_Mass_Of_Particles(const T_OBJECT& volume,const T density,const bool use_constant_mass=false);
    template <int gauss_order> static void Set_Mass_Of_Particles(const OPENSUBDIV_SURFACE<TV,gauss_order>& volume,const T density,const bool use_constant_mass=false);
    void PD_Curl(const T scale,const TV shift,const ROTATION<TV> orient,const T k_p,const int number_of_joints,const bool parent_static=true,const T friction=.5);
    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_BODY_STATE<TV>* initial_state,
        std::function<void (TRIANGULATED_SURFACE<T>&)> clipping_function,ARRAY<int>* particle_indices);
    void Embed_Particles_In_Tetrahedralized_Volume(BINDING_LIST<VECTOR<T,3> >& binding_list,const PARTICLES_SUBSET<VECTOR<T,3>,
        DEFORMABLE_PARTICLES<VECTOR<T,3> > >& particles_to_embed,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const T thickness_over_two);
    void Mark_Hard_Bindings_With_Free_Particles();
    void Find_Intersected_Segments_Triangles(SEGMENTED_CURVE<TV>& segments,TRIANGULATED_SURFACE<T>& surface,ARRAY<bool>* segments_intersected,ARRAY<bool>* triangles_intersected,
        T thickness_over_two,ARRAY<T>* segment_weights,ARRAY<VECTOR<T,3> >* triangle_weights);
    void Embed_Surface_In_Tetrahedralized_Volume(BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,TRIANGULATED_SURFACE<T>& surface,
        TETRAHEDRALIZED_VOLUME<T>& volume,const T thickness_over_two,ARRAY<int>& surface_particle_map,ARRAY<int>& volume_particle_map,
        TRIANGULATED_SURFACE<T>** new_surface,TETRAHEDRALIZED_VOLUME<T>** new_volume,bool bind_edges);
    void Create_Regular_Embedded_Surface(BINDING_LIST<TV>& binding_list,SOFT_BINDINGS<TV>& soft_bindings,TRIANGULATED_SURFACE<T>& surface,T density,int approx_volume,
        const T thickness_over_two,ARRAY<int>& surface_particle_map,TRIANGULATED_SURFACE<T>** new_surface,TETRAHEDRALIZED_VOLUME<T>** new_volume,bool bind_edges);
//#####################################################################
};
}
#endif
