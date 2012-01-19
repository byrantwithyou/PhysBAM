//#####################################################################
// Copyright 2004-2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_EXAMPLE
//#####################################################################
#ifndef __FACE_EXAMPLE__
#define __FACE_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_FACE_3D.h>
#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
namespace PhysBAM{

template <class T,class RW>
class QUASISTATIC_FACE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::solids_parameters;using BASE::fluids_parameters;
    using BASE::output_directory;using BASE::data_directory;using BASE::verbose;using BASE::frame_rate;

    std::string model_directory,input_directory;
    int number_of_muscles,number_of_attachments;
    ARRAY<ARRAY<int> > muscle_tets;
    ARRAY<ARRAY<VECTOR_3D<T> > > muscle_fibers;
    ARRAY<ARRAY<T> > muscle_densities;
    ARRAY<T> muscle_activations;
    ARRAY<ARRAY<int> > attached_nodes;
    ARRAY<VECTOR_3D<T> > X_save;
    ARRAY<T> initial_activations;
    ARRAY<T> final_activations;
    bool collide_with_cranium,collide_with_jaw;

    QUASISTATIC_FACE_EXAMPLE()
        :BASE(0,fluids_parameters.NONE)
    {
        solids_parameters.quasistatic=true;

        //model_directory=data_directory+"/Face_Data/Eftychis_1620k";input_directory=data_directory+"/Face_Data/Eftychis_1620k/Front_720k";
        model_directory=data_directory+"/Face_Data/Eftychis_840k";input_directory=data_directory+"/Face_Data/Eftychis_840k/Front_370k";
        output_directory="Quasistatic_Face/Output";

        last_frame=10*24;
        restart=false;restart_frame=0;
        
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-4;
        solids_parameters.implicit_solve_parameters.cg_iterations=200;
        solids_parameters.newton_tolerance=(T)1e-3;
        solids_parameters.newton_iterations=10;
        solids_parameters.perform_self_collision=false;
        solids_parameters.perform_collision_body_collisions=false;
        solids_parameters.use_partially_converged_result=true;

        solids_parameters.deformable_body_parameters.print_residuals=true;
        verbose=true;

        initial_activations.Resize(32);final_activations.Resize(32);
        final_activations(20)=5e1;
        final_activations(21)=5e1;
        final_activations(22)=5e1;
        final_activations(23)=5e1;

        collide_with_cranium=true;
        collide_with_jaw=true;
    }

    ~QUASISTATIC_FACE_EXAMPLE()
    {}

    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}    
    void Add_External_Forces(ARRAY<VECTOR_3D<T> >& F,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();
    if(solids_parameters.perform_collision_body_collisions) Initialize_Tetrahedron_Collisions();

    solids_parameters.deformable_body_parameters.list.deformable_objects(1)->Add_Force(Create_Quasistatic_Diagonalized_Face(
        *solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume,muscle_tets,muscle_fibers,muscle_densities,muscle_activations));

    if(solids_parameters.perform_collision_body_collisions){
        std::cout<<"Initializing collision forces"<<std::endl;
        COLLISION_PENALTY_FORCES<T> *penalty_force=new COLLISION_PENALTY_FORCES<T>(solids_parameters.deformable_body_parameters.list.deformable_objects(1)->tetrahedralized_volume->particles);
        penalty_force->Set_Stiffness((T)5e2);
        penalty_force->Set_Separation_Parameter((T)3e-4);
        penalty_force->Set_Collision_Body_List(solids_parameters.collision_body_list);
        penalty_force->Set_Collision_Body_List_ID(solids_parameters.rigid_body_parameters.list.rigid_bodies.m+1);
        penalty_force->Set_Boundary_Only_Collisions(solids_parameters.deformable_body_parameters.list.deformable_objects(1)->tetrahedralized_volume->tetrahedron_mesh);
        solids_parameters.deformable_body_parameters.list.deformable_objects(1)->collision_penalty_forces.Append(penalty_force);
        solids_parameters.deformable_body_parameters.list.deformable_objects(1)->solids_forces.Append(penalty_force);}
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    solids_parameters.deformable_body_parameters.list.Add_Deformable_Object();
    solids_parameters.deformable_body_parameters.list(1).Allocate_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T> &tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;
    std::string input_file;
    std::istream *input;

    input_file=input_directory+"/face_simulation.tet";input=FILE_UTILITIES::Safe_Open_Input(input_file);
    tetrahedralized_volume.template Read<RW>(*input);delete input;
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;
    std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Update_Velocity();particles.Store_Mass();
    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);tetrahedralized_volume.Update_Bounding_Box();

    input_file=input_directory+"/muscle_list.txt";input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    (*input)>>number_of_muscles;std::cout<<"muscles = "<<number_of_muscles<<std::endl;
    muscle_tets.Resize(number_of_muscles);muscle_fibers.Resize(number_of_muscles);
    muscle_densities.Resize(number_of_muscles);muscle_activations.Resize(number_of_muscles);
    for(int i=0;i<number_of_muscles;i++){
        std::string muscle_name;(*input)>>muscle_name;input_file=input_directory+"/"+muscle_name+".constitutive_data";
        if(verbose)std::cout<<"Reading muscle data : "<<muscle_name<<".constitutive_data"<<std::endl;
        std::istream *muscle_input=FILE_UTILITIES::Safe_Open_Input(input_file);
        muscle_tets(i).template Read<RW>(*muscle_input);muscle_fibers(i).template Read<RW>(*muscle_input);muscle_densities(i).template Read<RW>(*muscle_input);delete muscle_input;}
    delete input;

    input_file=input_directory+"/attachment_list.txt";input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    (*input)>>number_of_attachments;std::cout<<"attachments = "<<number_of_attachments<<std::endl;
    attached_nodes.Resize(number_of_attachments);
    for(int i=0;i<number_of_attachments;i++){
        std::string attachment_name;(*input)>>attachment_name;input_file=input_directory+"/"+attachment_name+".attached_nodes";
        if(verbose)std::cout<<"Reading attachment data : "<<attachment_name<<".attached_nodes"<<std::endl;
        std::istream *attachment_input=FILE_UTILITIES::Safe_Open_Input(input_file);
        attached_nodes(i).template Read<RW>(*attachment_input);}
    delete input;

    X_save=tetrahedralized_volume.particles.X.array;

    if(solids_parameters.perform_collision_body_collisions){
        if(collide_with_cranium){
            int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>(model_directory+"/eftychis_cranium_collision_surface");
            solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction(0);}
        if(collide_with_jaw){
            int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>(model_directory+"/eftychis_jaw_collision_surface");
            solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction(0);}
        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}
}
//#####################################################################
// Function Initialize_Tetrahedron_Collisions
//#####################################################################
void Initialize_Tetrahedron_Collisions()
{
    std::cout<<"Initializing geometric structures for self collision"<<std::endl;
    TETRAHEDRALIZED_VOLUME<T> &tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;
    TETRAHEDRON_COLLISION_BODY<T> *face_collision_body=new TETRAHEDRON_COLLISION_BODY<T>(tetrahedralized_volume);
    PARTICLES<T,VECTOR_3D<T> >* undeformed_particles=new PARTICLES<T,VECTOR_3D<T> >(tetrahedralized_volume.particles);
    TRIANGULATED_SURFACE<T> *undeformed_triangulated_surface=new TRIANGULATED_SURFACE<T>(tetrahedralized_volume.triangulated_surface->triangle_mesh,*undeformed_particles);
    undeformed_triangulated_surface->Update_Triangle_List();undeformed_triangulated_surface->Initialize_Triangle_Hierarchy();

    std::cout<<"Loading implicit surface for self collision"<<std::endl;
    LEVELSET_IMPLICIT_SURFACE<T>* undeformed_levelset;FILE_UTILITIES::Create_From_File<RW>(model_directory+"/face_full.phi",undeformed_levelset);
    undeformed_levelset->Update_Box();
    face_collision_body->Set_Implicit_Geometry(undeformed_levelset);
    face_collision_body->Set_Collision_Thickness((T)1e-3);
    face_collision_body->Set_Undeformed_Triangulated_Surface(undeformed_triangulated_surface);
    solids_parameters.collision_body_list.Add_Body(face_collision_body);
}
//#####################################################################
// Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY<VECTOR_3D<T> >& X,const T time) PHYSBAM_OVERRIDE {
    switch(id_number){
    case 1: for(int i=0;i<attached_nodes.m;i++)for(int j=1;j<=attached_nodes(i).m;j++) X(attached_nodes(i)(j))=X_save(attached_nodes(i)(j)); break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY<VECTOR_3D<T> >& X,const T time) PHYSBAM_OVERRIDE {
    switch(id_number){
    case 1: for(int i=0;i<attached_nodes.m;i++)for(int j=1;j<=attached_nodes(i).m;j++) X(attached_nodes(i)(j))=VECTOR_3D<T>(); break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {
    switch(id_number){
    case 1: for(int m=0;m<muscle_activations.m;m++)muscle_activations(m)=initial_activations(m)+(time*frame_rate/last_frame)*(final_activations(m)-initial_activations(m)); break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
};
}
#endif
