//#####################################################################
// Copyright 2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLESH_EXAMPLE
//#####################################################################
#ifndef __FLESH_EXAMPLE__
#define __FLESH_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
#include <Solids_And_Fluids/SOLIDS_PARAMETERS_3D.h>
namespace PhysBAM{

template <class T,class RW>
class FLESH_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    typedef VECTOR<T,3> TV;
public:
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::output_directory;using BASE::verbose;using BASE::frame_rate;
    
    std::string keyframe_prefix;
    T keyframe_rate;
    ARRAY<ARRAY<int> > constrained_nodes;
    ARRAY<ARRAY<TV> > constrained_nodes_reference_frame_positions;

    FLESH_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.NONE)
    {
        // File/Directory settings
        output_directory="Flesh/output";

        // Keyframe settings
        keyframe_prefix=data_directory+"/VH_Keyframes/arm_straightening/arm_straightening.11bones.frame";
        keyframe_rate=(T)24;

        // Simulation settings
        solids_parameters.quasistatic=true;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        solids_parameters.implicit_solve_parameters.cg_iterations=500;
        solids_parameters.newton_iterations=20;
        solids_parameters.newton_tolerance=(T)1e-2;
        solids_parameters.use_partially_converged_result=true;
        solids_parameters.perform_collision_body_collisions=true;
        solids_parameters.perform_self_collision=false;
        solid_body_collection.print_residuals=true;
        verbose=true;
    }
    
    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<VECTOR<T,3> > F,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    DEFORMABLE_OBJECT_3D<T>& deformable_object=solid_body_collection.deformable_object;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    solid_body_collection.Add_Force(Create_Quasistatic_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)5e3,(T).4,(T).01,(T).25),true,(T).1));
    solid_body_collection.Update_Fragments();

    // Collision settings
    if(solids_parameters.perform_collision_body_collisions){
        Initialize_Tetrahedron_Collisions();
        COLLISION_PENALTY_FORCES<T> *penalty_force=new COLLISION_PENALTY_FORCES<T>(solid_body_collection.deformable_object.particles);
        penalty_force->Set_Stiffness((T)1e4);
        penalty_force->Set_Separation_Parameter((T)1e-4);
        penalty_force->Set_Collision_Body_List(solids_parameters.collision_body_list);
        penalty_force->Set_Collision_Body_List_ID(1);
        penalty_force->Set_Boundary_Only_Collisions(tetrahedralized_volume.mesh);
        solid_body_collection.deformable_object.collision_penalty_forces.Append(penalty_force);
        solid_body_collection.solid_body_collection.solids_forces.Append(penalty_force);}

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solid_body_collection.deformable_object;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    FILE_UTILITIES::Read_From_File<RW>(data_directory+"/VH_Flesh/union_24_with_nipples.tet",tetrahedralized_volume);
    tetrahedralized_volume.particles.Update_Velocity();tetrahedralized_volume.particles.Store_Mass();
    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(false);
    deformable_object.Add_Structure(tetrahedralized_volume.Append_Particles_And_Create_Copy(deformable_object.particles));
    delete &tetrahedralized_volume;
    deformable_object.structures(1)->Update_Number_Nodes();

    // Rigid Bodies
    Initialize_Rigid_Bodies();

    // Initialize attachments
    FILE_UTILITIES::Read_From_File<RW>(data_directory+"/VH_Flesh/constrained_nodes",constrained_nodes);
    FILE_UTILITIES::Read_From_File<RW>(data_directory+"/VH_Flesh/constrained_nodes_reference_frame_positions",constrained_nodes_reference_frame_positions);
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Rigid_Bodies()
{
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/sternum_ribs",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/clavicle_right",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/scapula_right",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/humerus_right",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/ulna_right",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/radius_right",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/clavicle_left",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/scapula_left",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/humerus_left",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/ulna_left",(T)1);
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/Visible_Human_Bones/radius_left",(T)1);

    Update_Collision_Body_Positions_And_Velocities(1);
}
//#####################################################################
// Function Initialize_Tetrahedron_Collisions
//#####################################################################
void Initialize_Tetrahedron_Collisions()
{
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solid_body_collection.deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();

    TETRAHEDRON_COLLISION_BODY<T>& tetrahedron_collision_body=*(new TETRAHEDRON_COLLISION_BODY<T>(tetrahedralized_volume));
    PARTICLES<TV>& undeformed_particles=*(new PARTICLES<TV>(tetrahedralized_volume.particles));
    TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface=*(new TRIANGULATED_SURFACE<T>(tetrahedron_collision_body.triangulated_surface->mesh,undeformed_particles));
    undeformed_triangulated_surface.Update_Triangle_List();undeformed_triangulated_surface.Initialize_Triangle_Hierarchy();

    // undeformed levelset
    LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    FILE_UTILITIES::Read_From_File<RW>(data_directory+"/VH_Flesh/union_24_with_nipples.phi",undeformed_levelset);
    undeformed_levelset.Update_Box();
    tetrahedron_collision_body.Set_Collision_Thickness((T)1e-3);
    tetrahedron_collision_body.Set_Implicit_Geometry(&undeformed_levelset);
    tetrahedron_collision_body.Set_Undeformed_Triangulated_Surface(&undeformed_triangulated_surface);
    solids_parameters.collision_body_list.Add_Body(&tetrahedron_collision_body);
}
//#####################################################################
// Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int rb=0;rb<constrained_nodes.m;rb++){
        FRAME_3D<T> frame=solids_parameters.rigid_body_parameters.list(rb)->frame;
        for(int i=0;i<constrained_nodes(rb).m;i++) X(constrained_nodes(rb)(i))=frame*constrained_nodes_reference_frame_positions(rb)(i);}
}
//#####################################################################
// Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int rb=0;rb<constrained_nodes.m;rb++) for(int i=0;i<constrained_nodes(rb).m;i++) X(constrained_nodes(rb)(i))=TV();
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{   
    int frame=(int)(time*frame_rate+(T).5);
    std::istream *input=FILE_UTILITIES::Safe_Open_Input(keyframe_prefix+STRING_UTILITIES::string_sprintf(".%d",frame));
    for(int i=0;i<solids_parameters.rigid_body_parameters.list.rigid_bodies.m;i++)
        solids_parameters.rigid_body_parameters.list(i)->frame.template Read<RW>(*input);
    delete input;

    Set_External_Positions(solid_body_collection.deformable_object.particles.X.array,time,1);
}
//#####################################################################
};
}
#endif
