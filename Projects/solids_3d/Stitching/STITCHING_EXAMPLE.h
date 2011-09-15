//#####################################################################
// Copyright 2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STITCHING_EXAMPLE
//#####################################################################
#ifndef __STITCHING_EXAMPLE__
#define __STITCHING_EXAMPLE__

#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/IMPLICIT_ZERO_LENGTH_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TRIANGLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
namespace PhysBAM{

template<class T,class RW>
class STITCHING_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;
    int number_of_objects;

    STITCHING_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),tests(*this,solids_parameters)
    {
        output_directory="Stitching/output";
        last_frame=192;
        frame_rate=24;

        number_of_objects=8;

        solids_parameters.perform_self_collision=false;
    }

    ~STITCHING_EXAMPLE()
    {}

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    // deformable bodies
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<T,TV>& particles=deformable_object.particles;
    BINDING_LIST<T,TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<T,TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    for(int i=1;i<=number_of_objects;i++)
        tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME_3D<T>(TV((T)1.7*(i-4)-(T).85,(T)6,0))),true,true);
    for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    ARRAY<TETRAHEDRALIZED_VOLUME<T>*> tetrahedralized_volumes(deformable_object.structures.m);
    for(int i=1;i<=deformable_object.structures.m;i++){
        tetrahedralized_volumes(i)=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);
        tetrahedralized_volumes(i)->Update_Tetrahedron_List();tetrahedralized_volumes(i)->Initialize_Tetrahedron_Hierarchy();}

    particles.Set_Array_Buffer_Size(1000);
    ARRAY<int> map_to_hard_bound_particles(particles.array_collection->Size());
    for(int volume=1;volume<=tetrahedralized_volumes.m;volume++){
        int particles_embedded=0;
        ARRAY<int> nodes_of_structure;tetrahedralized_volumes(volume)->mesh.elements.Flattened().Get_Unique(nodes_of_structure);
        for(int i=1;i<=nodes_of_structure.m;i++){
            TV& X=deformable_object.particles.X(nodes_of_structure(i));
            ARRAY<VECTOR<int,2> > candidate_embeddings;
            for(int other_volume=1;other_volume<=deformable_object.structures.m;other_volume++) if(other_volume!=volume){
                ARRAY<int> intersection_list;(*tetrahedralized_volumes(other_volume)->tetrahedron_hierarchy).Intersection_List(X,intersection_list);
                for(int t=1;t<=intersection_list.m;t++) candidate_embeddings.Append(VECTOR<int,2>(other_volume,intersection_list(t)));}
            for(int e=candidate_embeddings.m;e>=1;e--){int other_volume,tet;candidate_embeddings(e).Get(other_volume,tet);
                if(!(*tetrahedralized_volumes(other_volume)->tetrahedron_list)(tet).Inside(X)) candidate_embeddings.Remove_Index_Lazy(e);}
            if(candidate_embeddings.m>1) PHYSBAM_FATAL_ERROR();
            if(candidate_embeddings.m==1){int other_volume,tet;candidate_embeddings(1).Get(other_volume,tet);
                VECTOR<T,3> weights;weights=(*tetrahedralized_volumes(other_volume)->tetrahedron_list)(tet).Barycentric_Coordinates(X);
                int embedded_particle=particles.array_collection->Add_Element();
                binding_list.Add_Binding(new LINEAR_BINDING<T,TV,4>(particles,embedded_particle,tetrahedralized_volumes(other_volume)->mesh.elements(tet),weights));
                map_to_hard_bound_particles(nodes_of_structure(i))=embedded_particle;
                particles_embedded++;}}
        LOG::cout<<particles_embedded<<" particles from structure "<<volume<<" were embedded"<<std::endl;}
    SEGMENTED_CURVE<T,TV>& segmented_curve=*SEGMENTED_CURVE<T,TV>::Create(particles);
    for(int i=1;i<=map_to_hard_bound_particles.m;i++) if(map_to_hard_bound_particles(i)) segmented_curve.mesh.elements.Append(VECTOR<int,2>(i,map_to_hard_bound_particles(i)));
    deformable_object.Add_Structure(&segmented_curve);
    for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    tests.Add_Ground();
    tests.Add_Rigid_Body("sphere",(T)3,(T).5);

    // add structures and rigid bodies to collisions
    deformable_object.collisions.collision_structures.Append_Elements(deformable_object.structures);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    // correct number nodes
    for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    binding_list.Clear_Hard_Bound_Particles(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<T,TV>& particles=deformable_object.particles;

    Get_Initial_Data();

    for(int i=1;i<=number_of_objects;i++){
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(i);
        solid_body_collection.Add_Force(new GRAVITY<T,TV>(tetrahedralized_volume));
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)3e5,(T).45,(T).01,(T).25),true,(T).1));}

    SEGMENTED_CURVE<T,TV>& segmented_curve=deformable_object.template Find_Structure<SEGMENTED_CURVE<T,TV>&>();
    IMPLICIT_ZERO_LENGTH_SPRINGS<T,TV>& implicit_zero_restlength_springs=*Create_Edge_Zero_Length_Springs(segmented_curve.mesh,particles,(T)1e6,(T)1);
    solid_body_collection.Add_Force(&implicit_zero_restlength_springs);

    solid_body_collection.Update_Fragments();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
};
}
#endif
