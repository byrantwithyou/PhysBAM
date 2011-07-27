//#####################################################################
// Copyright 2006, Kevin Der, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_EXAMPLE
//#####################################################################
#ifndef __CUTTING_EXAMPLE__
#define __CUTTING_EXAMPLE__

#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
namespace PhysBAM{

template<class T,class RW>
class CUTTING_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;
    ARRAY<int> attached_nodes;

    CUTTING_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),tests(*this,solids_parameters)
    {
        output_directory=STRING_UTILITIES::string_sprintf("Cutting/output");
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-8;
        last_frame=240;
    }

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    //Initialize_Bodies_Without_Embedding();
    Initialize_Bodies_With_Embedding();
}
//#####################################################################
// Function Initialize_Bodies_Without_Embedding
//#####################################################################
void Initialize_Bodies_Without_Embedding()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;

    TETRAHEDRALIZED_VOLUME<T>& tet_volume=tests.Create_Tetrahedralized_Volume("../cutting_better_3d/Output_Orig/embedding_volume.tet",RIGID_BODY_STATE<TV>(FRAME_3D<T>(TV(0,0.5,0))),false,true);
    tet_volume.Update_Bounding_Box();
    //for(int p=1;p<=particles.array_collection->Size();p++) if(particles.X(p).y<tet_volume.bounding_box->ymin+(T).01) attached_nodes.Append(p);
    tests.Add_Ground();

    tests.Create_Triangulated_Object("../cutting_better_3d/Output_Orig/cutting_surface.tri",RIGID_BODY_STATE<TV>(FRAME_3D<T>(TV(0,0.5,0))),false,true);

    // add structures and rigid bodies to collisions
    deformable_object.collisions.collision_structures.Append_Elements(deformable_object.structures);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    solids_parameters.gravity_direction=TV(0,-1,0);

    // correct number nodes
    for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    solid_body_collection.Add_Force(new GRAVITY<TV>(tet_volume));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tet_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));

    solid_body_collection.Update_Fragments();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Initialize_Bodies_With_Embedding
//#####################################################################
void Initialize_Bodies_With_Embedding()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;

    TETRAHEDRALIZED_VOLUME<T>& tet_volume=tests.Create_Tetrahedralized_Volume("../cutting_3d/Output_Dup/embedding_volume.tet",RIGID_BODY_STATE<TV>(FRAME_3D<T>(TV(0,4.5,0))),false,true);
    tet_volume.Update_Bounding_Box();
    tet_volume.Set_Mass_Of_Particles(true);

    TRIANGULATED_SURFACE<T>& boundary_surface=*TRIANGULATED_SURFACE<T>::Create(tet_volume.particles);
    TRIANGLE_MESH boundary_mesh;
    FILE_UTILITIES::Read_From_File<RW>("../cutting_3d/Output_Dup/boundary_mesh.tri",boundary_mesh);
    boundary_surface.mesh.Initialize_Mesh(boundary_mesh);
    boundary_surface.Set_Mass_Of_Particles(true);
    boundary_surface.Update_Bounding_Box();
    solid_body_collection.deformable_object.Add_Structure(&boundary_surface);
    //for(int p=1;p<=particles.array_collection->Size();p++) if(particles.X(p).y>tet_volume.bounding_box->ymax-(T).01) attached_nodes.Append(p);

    ARRAY<ARRAY<int> > embedding_map;
    FILE_UTILITIES::Read_From_File<RW>("../cutting_3d/Output_Dup/embedding_map",embedding_map);
    ARRAY<ARRAY<T> > parent_weights;
    FILE_UTILITIES::Read_From_File<RW>("../cutting_3d/Output_Dup/embedding_weights",parent_weights);
    if(embedding_map.m!=parent_weights.m) PHYSBAM_FATAL_ERROR();
    for(int p=1;p<=embedding_map.m;p++){
        ARRAY<int>& parents=embedding_map(p);
        if(parents.m==2){
            VECTOR<T,2> weights(parent_weights(p)(1),parent_weights(p)(2));
            if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
            for(int i=1;i<=weights.m;i++) weights(i)=max((T)0,min(weights(i),(T)1));
            solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,p,VECTOR<int,2>(parents(1),parents(2)),weights));}
        if(parents.m==3){
            VECTOR<T,3> weights(parent_weights(p)(1),parent_weights(p)(2),parent_weights(p)(3));
            if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
            for(int i=1;i<=weights.m;i++) weights(i)=max((T)0,min(weights(i),(T)1));
            solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,p,VECTOR<int,3>(parents(1),parents(2),parents(3)),weights));}
        if(parents.m==4){
            VECTOR<T,3> weights(parent_weights(p)(1),parent_weights(p)(2),parent_weights(p)(3));
            if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
            for(int i=1;i<=weights.m;i++) weights(i)=max((T)0,min(weights(i),(T)1));
            solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,p,VECTOR<int,4>(parents(1),parents(2),parents(3),parents(4)),weights));}}
    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(boundary_surface,solid_body_collection.deformable_body_collection.binding_list,solid_body_collection.deformable_body_collection.soft_bindings);
    tests.Add_Ground();
    tests.Add_Rigid_Body("sphere",(T)2,(T)0.3);

    // add structures and rigid bodies to collisions
    deformable_object.collisions.collision_structures.Append(&boundary_surface);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    // correct number nodes
    for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    solid_body_collection.Add_Force(new GRAVITY<TV>(tet_volume));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tet_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));

    solid_body_collection.Update_Fragments();

    //for(int i=1;i<=deformable_object.particles_of_fragment.m;i++){
    //    TV center_of_mass=ARRAY<TV>::Average(deformable_object.particles.X.Subset(deformable_object.particles_of_fragment(i)));
    //    deformable_object.particles.X.Subset(deformable_object.particles_of_fragment(i))+=center_of_mass*10;
   // }

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=1;i<=attached_nodes.m;i++) V(attached_nodes(i))=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=1;i<=attached_nodes.m;i++) V(attached_nodes(i))=TV();
}
//#####################################################################
};
}
#endif
