//#####################################################################
// Copyright 2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_EXAMPLE
//#####################################################################
#ifndef __CUTTING_EXAMPLE__
#define __CUTTING_EXAMPLE__

#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_LINEAR_FVM_2D.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT_EVOLUTION_NEWMARK_INCOMPRESSIBLE.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_2D.h>
namespace PhysBAM{

template<class T,class RW>
class CUTTING_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
    typedef VECTOR<T,2> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
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
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<T,TV>& particles=deformable_object.particles;
    
    TRIANGULATED_AREA<T>& triangulated_area=tests.Create_Triangulated_Object("../cutting_better_2d/Output_Dup/embedding_area.tri2d",RIGID_BODY_STATE<TV>(FRAME_2D<T>(TV(0,10))),false,true);
    triangulated_area.Update_Bounding_Box();
    triangulated_area.Set_Mass_Of_Particles(true);
    
    SEGMENTED_CURVE_2D<T>* segmented_curve=SEGMENTED_CURVE_2D<T>::Create(triangulated_area.particles);
    FILE_UTILITIES::Read_From_File<RW>("../cutting_better_2d/Output_Dup/final_boundary_mesh.mesh",segmented_curve->mesh);
    
    solid_body_collection.deformable_object.Add_Structure(segmented_curve);

    ARRAY<ARRAY<int> > embedding_map;
    FILE_UTILITIES::Read_From_File<RW>("../cutting_better_2d/Output_Dup/embedding_map",embedding_map);
    for(int p=1;p<=embedding_map.m;p++){
        ARRAY<int>& parents=embedding_map(p);
        if(parents.m==2){
            VECTOR<T,2> weights=SEGMENT_2D<T>::Barycentric_Coordinates(particles.X(p),particles.X(parents(1)),particles.X(parents(2)));
            if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
            solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<T,TV,2>(particles,p,VECTOR<int,2>(parents(1),parents(2)),weights));}
        if(parents.m==3){
            VECTOR<T,3> weights=TRIANGLE_2D<T>::Barycentric_Coordinates(particles.X(p),particles.X(parents(1)),particles.X(parents(2)),particles.X(parents(3)));
            if(weights.Min()<0) LOG::cout<<"Negative barycentric coordinates on particle "<<p<<" : "<<weights<<", parents : "<<parents<<std::endl;
            solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<T,TV,3>(particles,p,VECTOR<int,3>(parents(1),parents(2),parents(3)),weights));}}
    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(*segmented_curve,solid_body_collection.deformable_body_collection.binding_list,solid_body_collection.deformable_body_collection.soft_bindings);
    tests.Add_Ground();
    
    // add structures and rigid bodies to collisions
    deformable_object.collisions.collision_structures.Append(segmented_curve);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    // correct number nodes
    for(int i=1;i<=deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    solid_body_collection.Add_Force(new GRAVITY<T,TV>(triangulated_area));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(triangulated_area,new DIAGONALIZED_NEO_HOOKEAN_2D<T>((T)5e4,(T).45,(T).01)));

    solid_body_collection.Update_Fragments();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
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
