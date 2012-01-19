//#####################################################################
// Copyright 2006-2007, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_FACE_EXAMPLE
//#####################################################################
#ifndef __EMBEDDED_FACE_EXAMPLE__
#define __EMBEDDED_FACE_EXAMPLE__

#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING_DYNAMIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/FACE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Motion/ACTIVATION_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/ATTACHMENT_FRAME_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T>
class EMBEDDED_FACE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >
{
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;typedef VECTOR<T,3> TV;
public:
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;
    using BASE::frame_rate;using BASE::stream_type;

    SOLIDS_STANDARD_TESTS<TV> tests;

    std::string model_directory,control_directory;
    int number_of_muscles,number_of_attachments;
    ARRAY<ARRAY<int> > muscle_tets;
    ARRAY<ARRAY<VECTOR<T,3> > > muscle_fibers;
    ARRAY<ARRAY<T> > muscle_densities;
    ARRAY<std::string> muscle_names;
    ARRAY<T> peak_isometric_stress;
    ARRAY<ARRAY<int> > attached_nodes;
    ARRAY<int> boundary_nodes;
    int jaw_attachment_index;

    FACE_CONTROL_PARAMETERS<T> control_parameters;
    ACTIVATION_CONTROL_SET<T> *activation_controls;
    ATTACHMENT_FRAME_CONTROL_SET<T> *attachment_frame_controls;
    CONSTITUTIVE_MODEL<T,3> *face_constitutive_model;
    int control_frame_rate;
    int previous_frame;
    T previous_time;
    ARRAY<TV> effective_V;
    SEGMENT_MESH binding_mesh;

    bool print_controls;

    EMBEDDED_FACE_EXAMPLE(const STREAM_TYPE stream_type,const PARSE_ARGS& parse_args)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solids_parameters)
    {
        // simulation parameters
        control_frame_rate=60;
        frame_rate=60;
        solids_parameters.verbose=true;
        solids_parameters.verbose_dt=false;
        print_controls=false;

        // directory parameters

        last_frame=108;
        output_directory="Embedded_Face/output";
        control_directory=data_directory+"/Face_Data/Keyframes/Eftychis_1710k_930k_39muscles/we_can_also_model";

        //last_frame=108;
        //output_directory="Embedded_Face/output_we_can_also_model";
        //control_directory=data_directory+"/Face_Data/Keyframes/Eftychis_1710k_930k_39muscles/we_can_also_model";

        //last_frame=123;
        //output_directory="Embedded_Face/output_interactions_with_objects";
        //control_directory=data_directory+"/Face_Data/Keyframes/Eftychis_1710k_930k_39muscles/interactions_with_objects";

        //last_frame=89;
        //output_directory="Embedded_Face/output_during_speech";
        //control_directory=data_directory+"/Face_Data/Keyframes/Eftychis_1710k_930k_39muscles/during_speech";

        model_directory=data_directory+"/Face_Data/Eftychis_155k_embedded_with_eyelids";

        // evolution parameters
        solids_parameters.cfl=(T)1;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-1;
        solids_parameters.implicit_solve_parameters.cg_iterations=50;
        solids_parameters.perform_self_collision=false;
        solids_parameters.perform_collision_body_collisions=true;
        solid_body_collection.print_residuals=false;
        dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution)->throw_exception_on_backward_euler_failure=false;
    }

    ~EMBEDDED_FACE_EXAMPLE()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    // initialize simulation object
    TETRAHEDRALIZED_VOLUME<T>& input_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    FILE_UTILITIES::Read_From_File(stream_type,model_directory+"/face_embedding_volume.tet",input_tetrahedralized_volume);
    LOG::cout<<"Adding Embedding Volume - Tetrahedra = "<<input_tetrahedralized_volume.mesh.elements.m<<std::endl;
    input_tetrahedralized_volume.Set_Density((T)1000);input_tetrahedralized_volume.Set_Mass_Of_Particles(false);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Copy_And_Add_Structure(input_tetrahedralized_volume);
    tetrahedralized_volume.Update_Number_Nodes();
    int number_of_original_particles=particles.array_collection->Size();

    // determine simulation subset
    ARRAY<bool> particle_is_included(number_of_original_particles);
    for(int p=0;p<particles.array_collection->Size();p++) if((particles.X(p)-TV(0,0,.18)).Magnitude()<(T).18) particle_is_included(p)=true;

    ARRAY<bool> particle_on_simulation_boundary(number_of_original_particles);ARRAY<int> map_to_old_elements;
    for(int t=0;t<tetrahedralized_volume.mesh.elements.m;t++){VECTOR<int,4>& element=tetrahedralized_volume.mesh.elements(t);
        if(!particle_is_included.Subset(element).Contains(false)) map_to_old_elements.Append(t);
        else for(int i=0;i<4;i++) if(particle_is_included(element(i))) particle_on_simulation_boundary(element(i))=true;}

    ARRAY<int> map_to_new_elements(tetrahedralized_volume.mesh.elements.m);map_to_new_elements.Subset(map_to_old_elements)=IDENTITY_ARRAY<>(map_to_old_elements.m);
    for(int p=0;p<number_of_original_particles;p++) if(particle_on_simulation_boundary(p)) boundary_nodes.Append(p);
    ARRAY<VECTOR<int,4> > old_elements=tetrahedralized_volume.mesh.elements;
    tetrahedralized_volume.mesh.elements=old_elements.Subset(map_to_old_elements);

    // initialize embedding volume bindings (T-junctions)
    ARRAY<int> bound_particles;ARRAY<ARRAY<int> > bound_particle_parents;ARRAY<ARRAY<T> > bound_particle_weights;
    FILE_UTILITIES::Read_From_File(stream_type,model_directory+"/face_embedding_volume_bindings",bound_particles,bound_particle_parents,bound_particle_weights);
    for(int i=0;i<bound_particles.m;i++)
        if(particle_is_included.Subset(bound_particle_parents(i)).Contains(false)){
            if(particle_is_included.Subset(bound_particle_parents(i)).Contains(true)){            
                particle_on_simulation_boundary(bound_particles(i))=true;boundary_nodes.Append(bound_particles(i));}}
        else switch(bound_particle_parents(i).m){
            case 2:
                binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,bound_particles(i),
                        VECTOR<int,2>(bound_particle_parents(i)(1),bound_particle_parents(i)(2)),
                        VECTOR<T,2>(bound_particle_weights(i)(1),bound_particle_weights(i)(2))));
                break;
            case 3:
                binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,bound_particles(i),
                        VECTOR<int,3>(bound_particle_parents(i)(1),bound_particle_parents(i)(2),bound_particle_parents(i)(3)),
                        VECTOR<T,3>(bound_particle_weights(i)(1),bound_particle_weights(i)(2),bound_particle_weights(i)(3))));
            break;
            case 4:
                binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,bound_particles(i),
                        VECTOR<int,4>(bound_particle_parents(i)(1),bound_particle_parents(i)(2),bound_particle_parents(i)(3),bound_particle_parents(i)(4)),
                        VECTOR<T,4>(bound_particle_weights(i)(1),bound_particle_weights(i)(2),bound_particle_weights(i)(3),bound_particle_weights(i)(4))));
                break;
            default: PHYSBAM_FATAL_ERROR();}

    // distribute mass just for T-junction bindings
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);

    // initialize embedded surface
    EMBEDDING<TV>& embedding=*new EMBEDDING<TV>(particles);
    ARRAY<ARRAY<int> > embedding_parents;ARRAY<ARRAY<T> > embedding_weights;
    FILE_UTILITIES::Read_From_File(stream_type,model_directory+"/face_embedded_surface_bindings",embedding_parents,embedding_weights);
    int particle_offset=particles.array_collection->Size();particles.array_collection->Add_Elements(embedding_parents.m);
    tetrahedralized_volume.Update_Number_Nodes();
    for(int i=0;i<embedding_parents.m;i++){
        LINEAR_BINDING_DYNAMIC<TV>* binding=new LINEAR_BINDING_DYNAMIC<TV>(particles,particle_offset+i,embedding_parents(i).m);
        binding->parents=embedding_parents(i);binding->weights=embedding_weights(i);
        binding->Clamp_To_Embedded_Position();binding->Clamp_To_Embedded_Velocity();
        if(!particle_is_included.Subset(embedding_parents(i)).Contains(false)) binding_list.Add_Binding(binding);}
    FILE_UTILITIES::Read_From_File(stream_type,model_directory+"/face_embedded_surface_mesh",embedding.material_surface_mesh);
    for(int t=0;t<embedding.material_surface_mesh.elements.m;t++) embedding.material_surface_mesh.elements(t)+=particle_offset*VECTOR<int,3>::All_Ones_Vector();
    embedding.Update_Number_Nodes();
    deformable_object.Add_Structure(&embedding);

    // make soft bindings
    particles.Increase_Array_Size(binding_list.bindings.m);
    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,solid_body_collection.deformable_body_collection.soft_bindings);
    //soft_bindings.use_impulses_for_collisions.Fill(false);

    // initialize muscles
    std::string input_file=model_directory+"/muscle_list.txt";std::istream* input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    (*input)>>number_of_muscles;LOG::cout<<"muscles = "<<number_of_muscles<<std::endl;
    muscle_tets.Resize(number_of_muscles);muscle_fibers.Resize(number_of_muscles);muscle_densities.Resize(number_of_muscles);
    muscle_names.Resize(number_of_muscles),peak_isometric_stress.Resize(number_of_muscles);
    for(int i=0;i<number_of_muscles;i++){
        (*input)>>muscle_names(i);input_file=model_directory+"/"+muscle_names(i)+".constitutive_data";
        if(solids_parameters.verbose) LOG::cout<<"Reading muscle data : "<<muscle_names(i)<<".constitutive_data"<<std::endl;
        std::istream *muscle_input=FILE_UTILITIES::Safe_Open_Input(input_file);TYPED_ISTREAM typed_input(*muscle_input,stream_type);
        Read_Binary(typed_input,muscle_tets(i));Read_Binary(typed_input,muscle_fibers(i));Read_Binary(typed_input,muscle_densities(i));
        //muscle_tets(i).Read(stream_type,*muscle_input);muscle_fibers(i).Read(stream_type,*muscle_input);muscle_densities(i).Read(stream_type,*muscle_input);
        delete muscle_input;
        for(int t=muscle_tets(i).m;t>=1;t--)
            if(map_to_new_elements(muscle_tets(i)(t))) muscle_tets(i)(t)=map_to_new_elements(muscle_tets(i)(t));
            else{muscle_tets(i).Remove_Index_Lazy(t);muscle_fibers(i).Remove_Index_Lazy(t);muscle_densities(i).Remove_Index_Lazy(t);}}
    delete input;

    // initialize peak isometric stress
    peak_isometric_stress.Fill((T)3e5);
    input_file=model_directory+"/peak_isometric_stress_list.txt";input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    int peak_isometric_stress_entries;(*input)>>peak_isometric_stress_entries;
    for(int i=0;i<peak_isometric_stress_entries;i++){
        std::string muscle_name;T peak_isometric_stress_value;(*input)>>muscle_name>>peak_isometric_stress_value;
        int index;if(muscle_names.Find(muscle_name,index)) peak_isometric_stress(index)=peak_isometric_stress_value;}
    delete input;
    if(solids_parameters.verbose) for(int i=0;i<number_of_muscles;i++) LOG::cout<<muscle_names(i)<<" : Peak isometric stress = "<<peak_isometric_stress(i)<<std::endl;

    // initialize attachments
    ARRAY<ARRAY<int> > input_attached_nodes;
    input_file=model_directory+"/attachment_list.txt";input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    (*input)>>number_of_attachments;LOG::cout<<"attachments = "<<number_of_attachments<<std::endl;
    input_attached_nodes.Resize(number_of_attachments);
    jaw_attachment_index=0;
    for(int i=0;i<number_of_attachments;i++){
        std::string attachment_name;(*input)>>attachment_name;input_file=model_directory+"/"+attachment_name+".attached_nodes";
        if(solids_parameters.verbose) LOG::cout<<"Reading attachment data : "<<attachment_name<<".attached_nodes"<<std::endl;
        std::istream *attachment_input=FILE_UTILITIES::Safe_Open_Input(input_file);TYPED_ISTREAM typed_input(*attachment_input,stream_type);
        Read_Binary(typed_input,input_attached_nodes(i));if(attachment_name=="jaw") jaw_attachment_index=i;}
    delete input;
    if(!jaw_attachment_index) {LOG::cerr<<"ERROR: No jaw attachment specified"<<std::endl;exit(1);}

    // remap attachments to soft bindings, and create collision particles
    ARRAY<int> binding_attachments(embedding_parents.m);
    FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create(particles);
    for(int i=0;i<input_attached_nodes.m;i++) binding_attachments.Subset(input_attached_nodes(i)).Fill(i);
    attached_nodes.Resize(input_attached_nodes.m);
    for(int b=0;b<soft_bindings.bindings.m;b++){
        int particle_index,parent_index;soft_bindings.bindings(b).Get(particle_index,parent_index);
        if(int attachment=binding_attachments(parent_index-number_of_original_particles)){
            attached_nodes(attachment).Append(particle_index);
            soft_bindings.use_impulses_for_collisions(b)=false;}
        else free_particles.nodes.Append(particle_index);}
    deformable_object.Add_Structure(&free_particles);

    // initialize binding mesh for boundary conditions
    soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    soft_bindings.Initialize_Binding_Mesh(true);binding_mesh.Initialize_Mesh(*soft_bindings.binding_mesh);

    // correct number nodes (one last time)
    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // add structures and rigid bodies to collisions
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body(stream_type,model_directory+"/eftychis_cranium_collision_surface");
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body(stream_type,model_directory+"/eftychis_jaw_collision_surface");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Set_Coefficient_Of_Friction(0);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->Set_Coefficient_Of_Friction(0);
    //solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    //deformable_object.collisions.collision_structures.Append(&embedding);
    //deformable_object.collisions.collision_structures.Append(&tetrahedralized_volume);
    deformable_object.collisions.collision_structures.Append(&free_particles);

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    Get_Initial_Data();

    // initialize controls

    activation_controls=new ACTIVATION_CONTROL_SET<T>;for(int i=0;i<number_of_muscles;i++) activation_controls->Add_Activation(muscle_names(i));
    control_parameters.list.Append(activation_controls);
    FILE_UTILITIES::Write_To_File<T>(output_directory+"/face_control_set_1.muscle_names",muscle_names);
    attachment_frame_controls=new ATTACHMENT_FRAME_CONTROL_SET<T>(deformable_object.particles.X.array,attached_nodes,jaw_attachment_index);
    attachment_frame_controls->Read_Jaw_Joint_From_File(stream_type,model_directory+"/jaw_joint_parameters");
    control_parameters.list.Append(attachment_frame_controls);
    attachment_frame_controls->Write_Jaw_Joint_To_File(stream_type,output_directory+"/face_control_set_2.jaw_joint_parameters");
    std::ostream *output=FILE_UTILITIES::Safe_Open_Output(output_directory+"/face_control_set_types");
    control_parameters.Write_Configuration_To_File(stream_type,output_directory+"/face_control_parameters_configuration");
    control_parameters.Write_Control_Set_Types(stream_type,*output);delete output;
    attachment_frame_controls->Set_Original_Attachment_Configuration(solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Frame(),solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->Frame());

    // initialize forces
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    FINITE_VOLUME<TV,3>& fvm=*Create_Face(tetrahedralized_volume,muscle_tets,muscle_fibers,muscle_densities,activation_controls->activations,&peak_isometric_stress,(T)6e4,(T)2e4,(T)1e5,(T).15,(T).1,true,(T).1,true,false,true);
    //FINITE_VOLUME<TV,3>& fvm=*Create_Face(tetrahedralized_volume,muscle_tets,muscle_fibers,muscle_densities,activation_controls->activations,&peak_isometric_stress,(T)6e3,(T)2e3,(T)1e4,(T).05,(T).1,true,(T).1,true,true,true);
    solid_body_collection.Add_Force(&fvm);

    //solid_body_collection.Add_Force(new GRAVITY<T,TV>(tetrahedralized_volume));
    solid_body_collection.Add_Force(Create_Edge_Binding_Springs(binding_mesh,particles,(T)1e1,(T)1));

    activation_controls->muscle_force=&fvm;attachment_frame_controls->muscle_force=&fvm;

    // TODO: initialize MPI

    // update fragments and snap embedded surface to its target location
    solid_body_collection.Update_Fragments();
    if(deformable_object.particles_of_fragment.Size()!=FRAGMENT_ID(1)) PHYSBAM_FATAL_ERROR();
    soft_bindings.Clamp_Particles_To_Embedded_Positions(FRAGMENT_ID(1),false);soft_bindings.Clamp_Particles_To_Embedded_Velocities(FRAGMENT_ID(1),false);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Write_Output_Files(frame);
    Write_Controls(output_directory+"/",frame);
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Read_Output_Files_Solids(frame);
    Read_Controls(output_directory+"/",frame);
}
//#####################################################################
// Function Write_Controls
//#####################################################################
void Write_Controls(const std::string& prefix,const int frame) const
{
    std::string f=STRING_UTILITIES::Value_To_String(frame);FILE_UTILITIES::Write_To_File(stream_type,prefix+"controls."+f,control_parameters);
}
//#####################################################################
// Function Read_Controls
//#####################################################################
void Read_Controls(const std::string& prefix,const int frame)
{
    std::string f=STRING_UTILITIES::Value_To_String(frame);FILE_UTILITIES::Read_From_File(stream_type,prefix+"controls."+f,control_parameters);
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{
    T control_frame=time*control_frame_rate;
    int previous_control_frame=1+(int)control_frame;
    T interpolation_fraction=control_frame-(int)control_frame;
    LOG::cout<<"Setting controls for time = "<<time<<std::endl;
    int number_of_particles=solid_body_collection.deformable_object.particles.array_collection->Size();
    ARRAY<TV> X_initial(number_of_particles),X_final(number_of_particles);effective_V.Resize(number_of_particles,false,false);
    Read_Controls(control_directory+"/",previous_control_frame);attachment_frame_controls->Set_Attachment_Positions(X_initial);control_parameters.Save_Controls();
    Read_Controls(control_directory+"/",previous_control_frame+1);attachment_frame_controls->Set_Attachment_Positions(X_final);control_parameters.Interpolate(interpolation_fraction);
    for(int i=0;i<attached_nodes.m;i++) for(int j=1;j<=attached_nodes(i).m;j++)
        effective_V(attached_nodes(i)(j))=(T)control_frame_rate*(X_final(attached_nodes(i)(j))-X_initial(attached_nodes(i)(j)));
    if(print_controls) control_parameters.Print_Diagnostics();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<boundary_nodes.m;i++) V(boundary_nodes(i))=TV();
    for(int i=0;i<attached_nodes.m;i++) for(int j=1;j<=attached_nodes(i).m;j++) V(attached_nodes(i)(j))=effective_V(attached_nodes(i)(j));
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    attachment_frame_controls->Set_Attachment_Positions(X);
}
//##########################################y###########################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_velocity_time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<boundary_nodes.m;i++) V(boundary_nodes(i))=TV();
    for(int i=0;i<attached_nodes.m;i++) for(int j=1;j<=attached_nodes(i).m;j++) V(attached_nodes(i)(j))=TV();
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{   
    solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Set_Frame(attachment_frame_controls->Cranium_Frame());
    solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->Set_Frame(attachment_frame_controls->Jaw_Frame());
}
//#####################################################################
};
}
#endif
