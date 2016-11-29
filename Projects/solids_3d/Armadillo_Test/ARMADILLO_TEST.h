//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// 1 one mesh, separate fragments
//#####################################################################
#ifndef __ARMADILLO_TEST__
#define __ARMADILLO_TEST__
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET_MAKER.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Forces/BINDING_SPRINGS.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>

namespace PhysBAM{

template<class T_input>
class ARMADILLO_TEST:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    HASHTABLE<int> constrained_nodes;
    bool use_bindings;
    bool use_soft_bindings;
    bool use_impulses;
    bool use_edge_springs;
    bool use_binding_springs;

    ARMADILLO_TEST(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),use_bindings(false),use_soft_bindings(false),
        use_impulses(false),use_edge_springs(false),use_binding_springs(false)
    {
        parse_args.Parse();
        tests.data_directory=data_directory;
    }

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    //RIGID_BODY_PARTICLES<TV>& rigid_body_particles=deformable_body_collection.rigid_body_particles;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    
    LEVELSET_MAKER<T> levelset_maker;
    TRIANGULATED_SURFACE<T>& surface=tests.Create_Triangulated_Object(LOG::sprintf("%s/Triangulated_Surfaces/armadillo_original_400k.tri",data_directory.c_str()),
        RIGID_BODY_STATE<TV>(FRAME<TV>(VECTOR<T,3>(0,(T).6,0))),false,false,(T).01);
    surface.Update_Bounding_Box();RANGE<TV> box=*surface.bounding_box;int size_x=20;
    TV lengths=box.Edge_Lengths();TV aspect=lengths/lengths.x;
    VECTOR<int,3> M((T)size_x*aspect);GRID<TV> grid(M,box);GRID<TV> implicit_grid((VECTOR<int,3>((T)50*aspect)),box);

    // make level set 
    ARRAY<T,VECTOR<int,3> > phi;
    levelset_maker.Compute_Level_Set(surface,implicit_grid,phi);
    LEVELSET_IMPLICIT_OBJECT<TV> implicit(implicit_grid,phi);phi-=(T).001;

    // make a tet volume surrounding
    TETRAHEDRALIZED_VOLUME<T>& temp_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    temp_volume.Initialize_Octahedron_Mesh_And_Particles(grid);
    temp_volume.Discard_Tetrahedrons_Outside_Implicit_Surface(implicit);
    temp_volume.Discard_Valence_Zero_Particles_And_Renumber();
    LOG::cout<<"number of tets = "<<temp_volume.mesh.elements.m<<std::endl;

    TETRAHEDRALIZED_VOLUME<T>& volume=*(TETRAHEDRALIZED_VOLUME<T>*)temp_volume.Append_Particles_And_Create_Copy(particles);
    delete &temp_volume;
    deformable_body_collection.Add_Structure(&volume);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    if(use_bindings){
        // compute the embedding
        volume.Initialize_Hierarchy();
        ARRAY<int> tri_nodes;
        Get_Unique(tri_nodes,surface.mesh.elements.Flattened());
        ARRAY<int> list;
        for(int i=0;i<tri_nodes.m;i++){int p=tri_nodes(i);
            list.Remove_All();
            volume.hierarchy->Intersection_List(particles.X(p),list);
            for(int j=0;j<list.m;j++){int t=list(j);const VECTOR<int,4>& tet_nodes=volume.mesh.elements(t);
                TETRAHEDRON<T> tet(particles.X.Subset(tet_nodes));
                if(!tet.Outside(particles.X(p),(T)1e-3)){
                    //tets_with_embeddings.Set(t);
                    binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,p,tet_nodes,tet.Barycentric_Coordinates(particles.X(p))));
                    continue;}}}
        if(use_soft_bindings){
            // soft bindings
            particles.Preallocate(tri_nodes.m);
            HASHTABLE<int,int> hard_to_soft;
            for(int i=0;i<tri_nodes.m;i++){int hard=tri_nodes(i);
            int soft=particles.Add_Element();
            hard_to_soft.Insert(hard,soft);
            particles.X(soft)=particles.X(hard);
            soft_bindings.Add_Binding(VECTOR<int,2>(soft,hard),use_impulses);}
            for(int i=0;i<surface.mesh.elements.m;i++){
                VECTOR<int,3>& nodes=surface.mesh.elements(i);
                for(int k=0;k<nodes.m;k++) nodes[k]=hard_to_soft.Get(nodes[k]);}}}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();
        
    // set constrained nodes
    ARRAY<int> tet_nodes;
    Get_Unique(tet_nodes,volume.mesh.elements.Flattened());
    for(int i=0;i<tet_nodes.m;i++) if(particles.X(tet_nodes(i)).y<(T).17) constrained_nodes.Insert(tet_nodes(i));

    // tet mass
    T density=TV::m==1?1:TV::m==2?100:1000;
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(volume,density);

    // collision geom
    tests.Add_Ground();

    // add to collision structures
    deformable_body_collection.collisions.collision_structures.Append(&surface);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct mass
    //binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
void After_Initialization() override {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    frame_rate=24;
    last_frame=(int)(64*frame_rate);
    solids_parameters.verbose_dt=true;
    solids_parameters.cfl=(T)4;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    last_frame=20;
    switch(test_number){
        case 1: break;
        case 2: use_bindings=true;break;
        case 3: use_bindings=true;use_soft_bindings=true;use_impulses=true;break;
        case 4: use_bindings=true;use_soft_bindings=true;use_impulses=false;use_binding_springs=true;break;
    }
    LOG::printf("Parameters: %P, %P, %P, %P, %P\n",use_bindings,use_soft_bindings,use_impulses,use_edge_springs,use_binding_springs);
    
    // geometry
    Get_Initial_Data();

    // make forces
    for(int i=0;i<deformable_body_collection.structures.m;i++){
        if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(i))){
            ARRAY<int>& referenced_nodes=*new ARRAY<int>; // hey craig, look a memory leak.  cool andy, why don't you fix it?
            Get_Unique(referenced_nodes,volume->mesh.elements.Flattened());
            for(int i=referenced_nodes.m;i>=1;i--) if(constrained_nodes.Contains(referenced_nodes(i))) referenced_nodes.Remove_Index_Lazy(i);
            solid_body_collection.Add_Force(Create_Edge_Springs(*volume,(T)100,(T)3));
            solid_body_collection.Add_Force(Create_Altitude_Springs(*volume,(T)100,(T)3));
            solid_body_collection.Add_Force(new GRAVITY<TV>(particles,solid_body_collection.rigid_body_collection,&referenced_nodes,0));}}
    if(use_soft_bindings){
        soft_bindings.Initialize_Binding_Mesh();
        if(use_edge_springs){
            LINEAR_SPRINGS<TV>* linear_springs=Create_Edge_Springs(particles,*soft_bindings.binding_mesh,(T)1,(T)3);
            solid_body_collection.Add_Force(linear_springs);
            linear_springs->Clamp_Restlength((T)1e-2);}
        if(use_binding_springs){
            solid_body_collection.Add_Force(Create_Edge_Binding_Springs(particles,*soft_bindings.binding_mesh,(T)1000,(T)3));}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    for(HASHTABLE<int>::ITERATOR i(constrained_nodes);i.Valid();i.Next()) V(i.Key())=TV();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) override
{
    for(HASHTABLE<int>::ITERATOR i(constrained_nodes);i.Valid();i.Next()) V(i.Key())=TV();
}
//#####################################################################
};
}
#endif
