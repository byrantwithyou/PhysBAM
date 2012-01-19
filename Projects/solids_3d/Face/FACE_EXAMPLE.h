//#####################################################################
// Copyright 2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_EXAMPLE
//#####################################################################
#ifndef __FACE_EXAMPLE__
#define __FACE_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Motion/ACTIVATION_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/ATTACHMENT_FRAME_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_FACE_3D.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT_EVOLUTION_BACKWARD_EULER.h>
#include <Deformable_Objects/DEFORMABLE_OBJECT_EVOLUTION_QUASISTATIC.h>
namespace PhysBAM{

template<class T,class RW>
class FACE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;typedef VECTOR<T,3> TV;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::verbose;using BASE::verbose_dt;using BASE::frame_rate;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;

    bool use_quasistatics,use_backward_euler,use_tetrahedron_collisions;

    std::string input_directory,model_directory,control_directory;
    int number_of_muscles,number_of_attachments;
    ARRAY<ARRAY<int> > muscle_tets;
    ARRAY<ARRAY<VECTOR<T,3> > > muscle_fibers;
    ARRAY<ARRAY<T> > muscle_densities;
    ARRAY<std::string> muscle_names;
    ARRAY<T> peak_isometric_stress;
    ARRAY<ARRAY<int> > attached_nodes;
    int jaw_attachment_index;

    FACE_CONTROL_PARAMETERS<T> control_parameters;
    ACTIVATION_CONTROL_SET<T> *activation_controls;
    ATTACHMENT_FRAME_CONTROL_SET<T> *attachment_frame_controls;
    DIAGONALIZED_CONSTITUTIVE_MODEL_3D<T> *face_constitutive_model;
    int control_frame_rate;
    int previous_frame;
    T previous_time;
    ARRAY<TV> effective_V;

    FACE_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),tests(*this,solids_parameters)
    {
        // simulation parameters
        use_quasistatics=true;
        use_backward_euler=false;
        use_tetrahedron_collisions=false;
        control_frame_rate=24;
        frame_rate=24;
        last_frame=24;
        verbose=true;
        verbose_dt=true;

        // directory parameters
        output_directory="Face/output";
        model_directory=data_directory+"/Face_Data/Eftychis_1710k";
        input_directory=model_directory+"/Front_930k";
        control_directory=data_directory+"/Face_Data/Keyframes/Eftychis_1710k_930k_39muscles/phonemesA1";

        // evolution parameters
        if(use_quasistatics){
            solids_parameters.Set_Evolution(new DEFORMABLE_OBJECT_EVOLUTION_QUASISTATIC<TV>(solids_parameters));
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)5e-2;
            solids_parameters.implicit_solve_parameters.cg_iterations=200;
            solids_parameters.newton_tolerance=(T)1e-1;
            solids_parameters.newton_iterations=5;
            solids_parameters.use_partially_converged_result=true;
            solids_parameters.perform_self_collision=false;
            solids_parameters.perform_collision_body_collisions=false;
            solid_body_collection.print_residuals=true;}
        else if(use_backward_euler){
            solids_parameters.Set_Evolution(new DEFORMABLE_OBJECT_EVOLUTION_BACKWARD_EULER<TV>(solids_parameters));
            solids_parameters.cfl=200;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)2e3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100;
            solids_parameters.newton_tolerance=(T)5e3;
            solids_parameters.newton_iterations=5;
            solids_parameters.use_partially_converged_result=true;
            solids_parameters.perform_self_collision=false;
            solids_parameters.perform_collision_body_collisions=false;
            solid_body_collection.print_residuals=true;}
        else{
            solids_parameters.cfl=10;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-1;
            solids_parameters.implicit_solve_parameters.cg_iterations=100;
            solids_parameters.perform_self_collision=false;
            solids_parameters.perform_collision_body_collisions=false;
            solid_body_collection.print_residuals=true;}
    }

    ~FACE_EXAMPLE()
    {}

    // Unused callbacks
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    std::string input_file;
    std::istream *input;
    int peak_isometric_stress_entries;

    // initialize simulation object
    TETRAHEDRALIZED_VOLUME<T>& input_tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    FILE_UTILITIES::Read_From_File<RW>(input_directory+"/face_simulation.tet",input_tetrahedralized_volume);
    LOG::cout<<"Adding Tetrahedralized Volume - Tetrahedra = "<<input_tetrahedralized_volume.mesh.elements.m<<std::endl;
    input_tetrahedralized_volume.Set_Density((T)1000);input_tetrahedralized_volume.Set_Mass_Of_Particles(false);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Copy_And_Add_Structure(input_tetrahedralized_volume);
    tetrahedralized_volume.Update_Number_Nodes();

    // initialize muscles
    input_file=input_directory+"/muscle_list.txt";input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    (*input)>>number_of_muscles;LOG::cout<<"muscles = "<<number_of_muscles<<std::endl;
    muscle_tets.Resize(number_of_muscles);muscle_fibers.Resize(number_of_muscles);muscle_densities.Resize(number_of_muscles);
    muscle_names.Resize(number_of_muscles),peak_isometric_stress.Resize(number_of_muscles);
    for(int i=0;i<number_of_muscles;i++){
        (*input)>>muscle_names(i);input_file=input_directory+"/"+muscle_names(i)+".constitutive_data";
        if(verbose) LOG::cout<<"Reading muscle data : "<<muscle_names(i)<<".constitutive_data"<<std::endl;
        std::istream *muscle_input=FILE_UTILITIES::Safe_Open_Input(input_file);
        muscle_tets(i).template Read<RW>(*muscle_input);muscle_fibers(i).template Read<RW>(*muscle_input);muscle_densities(i).template Read<RW>(*muscle_input);delete muscle_input;}
    delete input;

    // initialize peak isometric stress
    ARRAY<T>::copy((T)3e5,peak_isometric_stress);
    input_file=input_directory+"/peak_isometric_stress_list.txt";input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    (*input)>>peak_isometric_stress_entries;
    for(int i=0;i<peak_isometric_stress_entries;i++){
        std::string muscle_name;T peak_isometric_stress_value;(*input)>>muscle_name>>peak_isometric_stress_value;
        int index;if(muscle_names.Find(muscle_name,index)) peak_isometric_stress(index)=peak_isometric_stress_value;}
    delete input;
    if(verbose) for(int i=0;i<number_of_muscles;i++) LOG::cout<<muscle_names(i)<<" : Peak isometric stress = "<<peak_isometric_stress(i)<<std::endl;

    // initialize attachments
    input_file=input_directory+"/attachment_list.txt";input=FILE_UTILITIES::Safe_Open_Input(input_file,false);
    (*input)>>number_of_attachments;LOG::cout<<"attachments = "<<number_of_attachments<<std::endl;
    attached_nodes.Resize(number_of_attachments);
    jaw_attachment_index=0;
    for(int i=0;i<number_of_attachments;i++){
        std::string attachment_name;(*input)>>attachment_name;input_file=input_directory+"/"+attachment_name+".attached_nodes";
        if(verbose) LOG::cout<<"Reading attachment data : "<<attachment_name<<".attached_nodes"<<std::endl;
        std::istream *attachment_input=FILE_UTILITIES::Safe_Open_Input(input_file);
        attached_nodes(i).template Read<RW>(*attachment_input);if(attachment_name=="jaw") jaw_attachment_index=i;}
    delete input;
    if(!jaw_attachment_index) {LOG::cerr<<"ERROR: No jaw attachment specified"<<std::endl;exit(1);}

    // collisions
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(model_directory+"/eftychis_cranium_collision_surface");
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(model_directory+"/eftychis_jaw_collision_surface");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Set_Coefficient_Of_Friction(0);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->Set_Coefficient_Of_Friction(0);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    if(use_tetrahedron_collisions)
        tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume);
    else
        solids_parameters.perform_self_collision=false;

    // complete computation of mass and auxiliary structures
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
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;

    Get_Initial_Data();

    // initialize controls
    activation_controls=new ACTIVATION_CONTROL_SET<T>;for(int i=0;i<number_of_muscles;i++)activation_controls->Add_Activation(muscle_names(i));
    control_parameters.list.Append(activation_controls);
    FILE_UTILITIES::Write_To_File<T>(output_directory+"/face_control_set_1.muscle_names",muscle_names);
    attachment_frame_controls=new ATTACHMENT_FRAME_CONTROL_SET<T>(deformable_object.particles.X.array,attached_nodes,jaw_attachment_index);
    attachment_frame_controls->template Read_Jaw_Joint_From_File<RW>(input_directory+"/jaw_joint_parameters");
    control_parameters.list.Append(attachment_frame_controls);
    attachment_frame_controls->template Write_Jaw_Joint_To_File<RW>(output_directory+"/face_control_set_2.jaw_joint_parameters");
    std::ostream *output=FILE_UTILITIES::Safe_Open_Output(output_directory+"/face_control_set_types");
    control_parameters.template Write_Configuration_To_File<RW>(output_directory+"/face_control_parameters_configuration");
    control_parameters.template Write_Control_Set_Types<RW>(*output);delete output;
    attachment_frame_controls->Set_Original_Attachment_Configuration(solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->frame,solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->frame);

    // initialize forces
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    if(use_quasistatics)
        solid_body_collection.Add_Force(Create_Quasistatic_Diagonalized_Face(tetrahedralized_volume,muscle_tets,muscle_fibers,muscle_densities,
                activation_controls->activations,&activation_controls->single_activation_used_for_force_derivative,&peak_isometric_stress,&face_constitutive_model,(T)6e4,(T)2e4,(T)1e7,(T).1));
    else{
        DIAGONALIZED_FINITE_VOLUME_3D<T>& fvm=*Create_Diagonalized_Face(tetrahedralized_volume,muscle_tets,muscle_fibers,muscle_densities,
            activation_controls->activations,&peak_isometric_stress,(T)6e4,(T)2e4,(T)1e7,(T).05,(T).1,true,(T).1,true,true,true);
        solid_body_collection.Add_Force(&fvm);
        if(use_backward_euler) fvm.Use_Quasistatics();}

    DIAGONALIZED_FINITE_VOLUME_3D<T>* diagonalized_finite_volume_3d=solid_body_collection.template Find_Force<DIAGONALIZED_FINITE_VOLUME_3D<T>*>();
    activation_controls->muscle_force=diagonalized_finite_volume_3d;
    attachment_frame_controls->muscle_force=diagonalized_finite_volume_3d;

    // initialize collision penalty forces
    if(solids_parameters.perform_collision_body_collisions){
        COLLISION_PENALTY_FORCES<TV> *penalty_force=new COLLISION_PENALTY_FORCES<TV>(particles);
        penalty_force->Set_Stiffness((T)1e4);
        penalty_force->Set_Separation_Parameter((T)1e-4);
        penalty_force->Set_Collision_Body_List(solids_parameters.collision_body_list);
        if(use_tetrahedron_collisions) penalty_force->Set_Collision_Body_List_ID(solids_parameters.rigid_body_parameters.list.rigid_bodies.m+1);
        penalty_force->Set_Boundary_Only_Collisions(tetrahedralized_volume.mesh);
        solid_body_collection.Add_Force(penalty_force);}

    // this is only correct in the absence of bindings (comment still necessary?)
    if(solid_body_collection.deformable_body_collection.mpi_solids){
        MPI_SOLIDS& mpi_solids=*solid_body_collection.deformable_body_collection.mpi_solids;
        FILE_UTILITIES::Read_From_File<T>(input_directory+STRING_UTILITIES::string_sprintf("/Partition_Data/partition_nodes_%d",mpi_solids.number_of_processes),
            mpi_solids.dynamic_particles_of_partition);}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Write_Output_Files(frame);Write_Controls(output_directory+"/",frame);
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Read_Output_Files_Solids(frame);Read_Controls(output_directory+"/",frame);
}
//#####################################################################
// Function Write_Controls
//#####################################################################
void Write_Controls(const std::string& prefix,const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);FILE_UTILITIES::Write_To_File<RW>(prefix+"controls"+f,control_parameters);
}
//#####################################################################
// Function Read_Controls
//#####################################################################
void Read_Controls(const std::string& prefix,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf(".%d",frame);FILE_UTILITIES::Read_From_File<RW>(prefix+"controls"+f,control_parameters);
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    previous_frame=frame-1;previous_time=previous_frame/frame_rate;
    Read_Controls(control_directory+"/",frame);control_parameters.Save_Controls();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
    T control_frame=time*control_frame_rate;
    int previous_control_frame=1+(int)control_frame;
    T interpolation_fraction=control_frame-(int)control_frame;
    LOG::cout<<"Setting controls for time = "<<time<<std::endl;
    int number_of_particles=solid_body_collection.deformable_object.particles.array_collection->Size();
    ARRAY<TV> X_initial(number_of_particles),X_final(number_of_particles);effective_V.Resize(number_of_particles,false,false);
    Read_Controls(control_directory+"/",previous_control_frame);attachment_frame_controls->Set_Attachment_Positions(X_initial);control_parameters.Save_Controls();
    Read_Controls(control_directory+"/",previous_control_frame+1);attachment_frame_controls->Set_Attachment_Positions(X_final);control_parameters.Interpolate(interpolation_fraction);
    for(int i=0;i<attached_nodes.m;i++) for(int j=0;j<attached_nodes(i).m;j++)
        effective_V(attached_nodes(i)(j))=(T)control_frame_rate*(X_final(attached_nodes(i)(j))-X_initial(attached_nodes(i)(j)));
    control_parameters.Print_Diagnostics();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<attached_nodes.m;i++) for(int j=0;j<attached_nodes(i).m;j++) V(attached_nodes(i)(j))=effective_V(attached_nodes(i)(j));
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
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<attached_nodes.m;i++) for(int j=0;j<attached_nodes(i).m;j++) V(attached_nodes(i)(j))=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<attached_nodes.m;i++) for(int j=0;j<attached_nodes(i).m;j++) X(attached_nodes(i)(j))=TV();
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{   
    solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->frame=attachment_frame_controls->Cranium_Frame();
    solids_parameters.rigid_body_parameters.list.rigid_bodies(2)->frame=attachment_frame_controls->Jaw_Frame();
}
//#####################################################################
};
}
#endif
