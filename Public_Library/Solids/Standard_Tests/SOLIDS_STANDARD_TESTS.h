//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_STANDARD_TESTS
//#####################################################################
#ifndef __SOLIDS_STANDARD_TESTS__
#define __SOLIDS_STANDARD_TESTS__

#include <Core/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <Core/Log/DEBUG_UTILITIES.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>
#include <Rigids/Standard_Tests/RIGIDS_STANDARD_TESTS.h>
#include <Deformables/Standard_Tests/DEFORMABLES_STANDARD_TESTS.h>
namespace PhysBAM{

template<class TV> class SOLID_BODY_COLLECTION;
template<class TV>
class SOLIDS_STANDARD_TESTS:public DEFORMABLES_STANDARD_TESTS<TV>,public RIGIDS_STANDARD_TESTS<TV>
{
    typedef DEFORMABLES_STANDARD_TESTS<TV> BASE;
    typedef typename BASE::T T;
    typedef typename BASE::T_TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
    typedef typename BASE::T_SEGMENTED_CURVE T_SEGMENTED_CURVE;
public:
    using BASE::Create_Mattress;using BASE::Create_Cloth_Panel;using BASE::Substitute_Soft_Bindings_For_Nodes;
    using RIGIDS_STANDARD_TESTS<TV>::Add_Rigid_Body;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;

    SOLIDS_STANDARD_TESTS(STREAM_TYPE stream_type,const std::string& data_directory,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input);
    virtual ~SOLIDS_STANDARD_TESTS(){}

    TRIANGULATED_AREA<T>& Create_Mattress(const GRID<TV>& mattress_grid,const bool use_constant_mass,const RIGID_BODY_STATE<TV>& initial_state,const T density)
    {return Create_Mattress(mattress_grid,use_constant_mass,&initial_state,density);}

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
    void Bind_Particles_In_Rigid_Body(RIGID_BODY<TV>& rigid_body);
    template<class T_ARRAY> void Bind_Unbound_Particles_In_Rigid_Body(RIGID_BODY<TV>& rigid_body,const T_ARRAY& particle_array);
    template<class T_ARRAY> void Bind_Particles_In_Rigid_Body(RIGID_BODY<TV>& rigid_body,const T_ARRAY& particle_array);
    void PD_Curl(const T scale,const FRAME<TV>& frame,const T k_p,const int number_of_joints,const bool parent_static=true,const T friction=.5);
    RIGID_BODY<TV>* Create_Rigid_Body_From_Tetrahedralized_Volume(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,
        const T density,const T cell_size,const int subdivision_loops=0,const bool perform_manifold_check=false,const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&)=0,
        const bool use_implicit_surface_maker=true);
    RIGID_BODY<TV>* Create_Rigid_Body_From_Fracture_Tetrahedralized_Volume(EMBEDDED_MATERIAL_SURFACE<TV,3>& tetrahedralized_volume,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T density,const T cell_size,const int subdivision_loops=0,const bool perform_manifold_check=false,
        const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&)=0,const bool use_implicit_surface_maker=true);
    RIGID_BODY<TV>* Create_Rigid_Body_From_Triangulated_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T density);
    RIGID_BODY<TV>* Create_Rigid_Body_From_Triangulated_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T density,int levelset_resolution);
    static RIGID_BODY<TV>* Create_Rigid_Body_From_Triangulated_Area(TRIANGULATED_AREA<T>& triangulated_area,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T density,
        const T cell_size,const bool move_body_to_center_of_mass=false,const bool move_only_mesh_particles=false);
//#####################################################################
};
}
#endif
