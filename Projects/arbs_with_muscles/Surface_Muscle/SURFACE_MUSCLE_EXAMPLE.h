//#####################################################################
// Copyright 2006-2007, Kevin Der, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_MUSCLE_EXAMPLE
//##################################################################### 
#ifndef __SURFACE_MUSCLE_EXAMPLE__
#define __SURFACE_MUSCLE_EXAMPLE__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/FACE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ANALYTIC_SURFACE_MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
#include "../ARB_PARAMETERS.h"
namespace PhysBAM{

template<class T>
class SURFACE_MUSCLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,3> > >
{
    typedef VECTOR<T,3> TV;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE; 
public:
    typedef typename MUSCLE_SEGMENT<TV>::MUSCLE_SEGMENT_TYPE T_MUSCLE_SEGMENT_TYPE;
    typedef typename ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_TYPE T_MUSCLE_SEGMENT_CURVE_TYPE;
    typedef TRIPLE<T_MUSCLE_SEGMENT_TYPE,T_MUSCLE_SEGMENT_CURVE_TYPE,ARRAY<T> > T_MUSCLE_SEGMENT_DATA;

    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::output_directory;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::stream_type;using BASE::solid_body_collection;

    SOLIDS_STANDARD_TESTS<TV> tests;

    ARTICULATED_RIGID_BODY<TV>* arb;
    MUSCLE<TV>* muscle;
    bool add_ground;
    PARAMETER_LIST parameter_list;
    T peak_force;
    GRID<TV> mattress_grid;

    // Data pertaining to enslaving flesh to bone
    int num_planks;
    ARRAY<int> plank_ids;
    ARRAY<ARRAY<int> > enslaved_nodes; // indicies of enslaved nodes. 
    ARRAY<ARRAY<TV> > positions_relative_to_plank_frames; // corresponding positions

    ARRAY<int> muscle_segment_particle_ids;
    ARRAY<int> tet_particle_ids;
    ARRAY<TV> muscle_segment_particle_original_positions;
    int number_of_muscle_segment_particles;

    INCOMPRESSIBLE_FINITE_VOLUME<TV,3>* finite_volume;
    ARRAY<ARRAY<int> > muscle_tets;
    ARRAY<ARRAY<TV> > muscle_fibers;
    ARRAY<ARRAY<T> > muscle_densities;
    ARRAY<T> tet_activations;
    ARRAY<MATRIX<T,3> > tet_F_o;
    ARRAY<T> peak_isometric_stress;
    ARRAY<T> activation_delta;
    ARRAY<T> previous_elongation;
    T previous_dt;

    ARRAY<bool> enlarge_nodes;

    SURFACE_MUSCLE_EXAMPLE(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,0,fluids_parameters.NONE),tests(*this,solids_parameters),muscle(0),
        mattress_grid(15,35,15,(T)-.2,(T).4,(T)-2.2,(T)1.10,(T)-.3,(T).3),number_of_muscle_segment_particles(0){
        last_frame=3000;
        frame_rate=96;
        output_directory="Surface_Muscle/output";

        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-4;
        solids_parameters.perform_self_collision=false;
        solids_parameters.perform_collision_body_collisions=false;
        //solids_parameters.cfl=(T).1;

        arb=new ARTICULATED_RIGID_BODY<TV>(solid_body_collection.deformable_body_collection.particles,solids_parameters.rigid_body_parameters.list);
        solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);
        arb->Set_Do_Final_Pass(false);
        arb->Use_Epsilon_Scale(false);
        arb->Set_Use_Shock_Propagation(false);
        arb->Use_Muscle_Actuators();

        ARB_PARAMETERS::Read_Common_Parameters("Surface_Muscle/example.param",*this,parameter_list);
        peak_force=parameter_list.Get_Parameter("peak_force",(T)10);
    }

    ~SURFACE_MUSCLE_EXAMPLE()
    {}

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    // Initialize rigid bodies and muscles
    arb->muscle_list->muscle_force_curve.Initialize(data_directory); // initialize here rather than constructor since data directory might be set after constructor

    Simple_Muscle_Across_Joint();
    std::cout <<"\nnumber of muscles is: "<<arb->muscle_list->muscles.m<<std::endl;

    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solids_parameters.rigid_body_parameters.list;
    tests.Add_Gravity();

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();

    // For now create all surface muscle segments
    for(int i=1;i<=arb->muscle_list->muscles.m;i++){
        ARRAY<T_MUSCLE_SEGMENT_DATA> segment_data;
        for (int j=1;j<=arb->muscle_list->muscles(i)->via_points.m+1;j++) {
            ARRAY<T> parameters;parameters.Append(0.04);parameters.Append(0.05);
            // tendon parameters
//                parameters.Append(1.0);parameters.Append(0.2);parameters.Append(0.3);
            parameters.Append(1.0);parameters.Append(0.01);parameters.Append(0.01);
            if(j==2) segment_data.Append(T_MUSCLE_SEGMENT_DATA(MUSCLE_SEGMENT<TV>::LINEAR_SEGMENT,ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_NONE,parameters));
            else segment_data.Append(T_MUSCLE_SEGMENT_DATA(MUSCLE_SEGMENT<TV>::ANALYTIC_SURFACE_SEGMENT,ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_COSINE,parameters));}
        T total_length=arb->muscle_list->muscles(i)->Total_Length();
        arb->muscle_list->muscles(i)->Set_Optimal_Length((T).8*total_length);
        arb->muscle_list->muscles(i)->Set_Tendon_Slack_Length((T).2*total_length);
        arb->muscle_list->muscles(i)->Set_Peak_Force(peak_force);
        arb->muscle_list->muscles(i)->Set_Max_Shortening_Velocity(1);
        arb->muscle_list->muscles(i)->Initialize(segment_data);}

    std::cout<<"number of muslce segments: "<<muscle->muscle_segments.m<<std::endl;
    std::cout<<"muscle segment type:       "<< muscle->muscle_segments(1)->segment_type<<std::endl;

    // Add in skin mesh
    Add_Skin_Mesh();

    Get_Constrained_Particle_Data();
}
//#####################################################################
// Function Add_Basic_Muscle
//#####################################################################
MUSCLE<TV>* Add_Basic_Muscle(const std::string& name,RIGID_BODY<TV>& origin_body,const TV& origin,RIGID_BODY<TV>& insertion_body,const TV& insertion,const T force_factor)
{
    MUSCLE<TV>* new_muscle=new MUSCLE<TV>(arb->muscle_list->particles,arb->muscle_list->muscle_force_curve);
    new_muscle->Set_Name(name);
// MAKE THESE RIGID BODY BINDINGS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    new_muscle->Set_Attachment_Point_1(new RIGID_BODY_BINDING<TV>(arb->muscle_list->particles,0,origin_body.rigid_body_particles,origin_body.particle_index,origin));
    new_muscle->Set_Attachment_Point_2(new RIGID_BODY_BINDING<TV>(arb->muscle_list->particles,0,insertion_body.rigid_body_particles,insertion_body.particle_index,insertion));

    //new_muscle->Set_Attachment_Point_1(new T_CONSTRAINED_POINT_IN_RIGID_BODY(origin_body,origin));
    //new_muscle->Set_Attachment_Point_2(new T_CONSTRAINED_POINT_IN_RIGID_BODY(insertion_body,insertion));
    arb->muscle_list->Add_Muscle(new_muscle);
    return new_muscle;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
//    enlarge_nodes.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
//    enlarge_nodes.Fill(false);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    LOG::cout<<"PREPROCESS_SOLIDS_SUBSTEP ---------------------------"<<std::endl;
    PARTICLES<TV> particles=solid_body_collection.deformable_body_collection.particles;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    ARRAY<bool> marked(particles.array_collection->Size());marked.Fill(false);

    FACE_3D<T>& face_constitutive_model=dynamic_cast<FACE_3D<T>&>(finite_volume->constitutive_model);

    // these parameters may be on a per muscle basis?

    T max_change_delta=5,min_change_delta=-5;

    for(int t=0;t<muscle_tets.m;t++){

        // calculate elongation
        MATRIX<T,3> F=finite_volume->strain_measure.F(muscle_tets(t)(1));
        F*=pow(F.Determinant(),-(T)one_third);

        int i,j,k,l;tetrahedralized_volume.mesh.elements.Get(muscle_tets(t)(1),i,j,k,l);
        ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>* analytic_segment=(ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>*)muscle->muscle_segments(1); // there should be only one segment for this muscle

        T elongation_i=analytic_segment->Get_Elongation_At_Local_Position(analytic_segment->inside_particle_rest_positions(i),analytic_segment->inside_particle_segments(i));
        T elongation_j=analytic_segment->Get_Elongation_At_Local_Position(analytic_segment->inside_particle_rest_positions(j),analytic_segment->inside_particle_segments(j));
        T elongation_k=analytic_segment->Get_Elongation_At_Local_Position(analytic_segment->inside_particle_rest_positions(k),analytic_segment->inside_particle_segments(k));
        T elongation_l=analytic_segment->Get_Elongation_At_Local_Position(analytic_segment->inside_particle_rest_positions(l),analytic_segment->inside_particle_segments(l));
        T desired_elongation=(T).25*(elongation_i+elongation_j+elongation_k+elongation_l);


        TV material_space_fiber_direction=face_constitutive_model.tet_fibers(muscle_tets(t)(1))(1).Normalized();
        T current_elongation=(F*material_space_fiber_direction).Magnitude();
        
        T v_target_times_dt=desired_elongation-current_elongation;
        T v_times_dt=current_elongation-previous_elongation(t);
        
        // TODO: damping coefficient--there would presumably be a damping coefficient in here as well which would be set on a per tet/muscle basis
        activation_delta(t)=5*(v_target_times_dt-v_times_dt);
        
        activation_delta(t)=clamp(activation_delta(t),min_change_delta,max_change_delta);

        LOG::cout<<"desired_elongation="<<desired_elongation<<std::endl;
        LOG::cout<<"current_elongation="<<current_elongation<<std::endl;

        if(current_elongation>desired_elongation) tet_activations(t)-=activation_delta(t);
        if(current_elongation<desired_elongation) tet_activations(t)+=activation_delta(t);
        // if current elongation is desired elongation, do nothing

        tet_activations(t)=clamp(tet_activations(t),(T)0,(T)1);
        previous_elongation(t)=current_elongation;
        LOG::cout<<"Activation "<<t<<": "<<tet_activations(t)<<std::endl;
    }

}
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    //TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
/*    if(finite_volume->rest_volumes.m>0 && frame < 2)
        for(int e=0;e<enlarge_nodes.m;e++){
            finite_volume->rest_volumes(e)*=1.1;
        }
*/
}

//#####################################################################
// Function Add_Skin_Mesh
//#####################################################################
void Add_Skin_Mesh()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Mattress(mattress_grid);
    tetrahedralized_volume.Update_Number_Nodes();

    // collisions
    // tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume);
    deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures.Last());

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // compute muscle tets
    ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>& muscle_segment=*static_cast<ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>*>(muscle->muscle_segments(1));
    for(int t=0;t<tetrahedralized_volume.mesh.elements.m;t++) for(int v=0;v<4;v++){
        if(muscle_segment.Analytic_Inside_Test(muscle_segment.frame.Inverse()*particles.X(tetrahedralized_volume.mesh.elements(v,t)))){
            muscle_tets.Resize(muscle_tets.m+1);muscle_tets.Last().Append(t);break;}}
    muscle_fibers.Resize(muscle_tets.m);muscle_densities.Resize(muscle_tets.m);tet_activations.Resize(muscle_tets.m);peak_isometric_stress.Resize(muscle_tets.m);
    activation_delta.Resize(muscle_tets.m);activation_delta.Fill(.1);previous_elongation.Resize(muscle_tets.m);
    // TODO: each muscle should contain each of these arrays; then muscles can set activations based on their own calculations
    for(int t=0;t<muscle_tets.m;t++){
        // SET FIBER DIRECTION
        muscle_fibers(t).Append(TV(0,1,0)); // just for testing--Should be specified in the muscle, perhaps based on pennation angle?
        muscle_densities(t).Append(1);tet_activations(t)=0;peak_isometric_stress(t)=(T)10e6;}

    finite_volume=Create_Incompressible_Face(tetrahedralized_volume,muscle_tets,muscle_fibers,muscle_densities,tet_activations,&peak_isometric_stress,(T)6e4,(T)2e4);
    solid_body_collection.deformable_body_collection.Add_Force(finite_volume);
    //solid_body_collection.deformable_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));

    muscle_segment.Initialize_Inside_Particles(tetrahedralized_volume);
    solid_body_collection.deformable_body_collection.Update_Fragments();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int p=0;p<num_planks;p++){
        //if(p==1) continue;
        RIGID_BODY<TV>& plank_rigid_body=*arb->rigid_body_list.rigid_bodies(plank_ids(p));
        for(int i=0;i<enslaved_nodes(p).m;i++){
//            X(enslaved_nodes(p)(i))=inverted_frame.r.Inverse_Rotate(positions_relative_to_plank_frames(p)(i)-inverted_frame.t);}}
            X(enslaved_nodes(p)(i))=plank_rigid_body.Frame()*positions_relative_to_plank_frames(p)(i);}}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int p=0;p<num_planks;p++){
        //if(p==1) continue;
        RIGID_BODY<TV>& plank_rigid_body=*arb->rigid_body_list.rigid_bodies(plank_ids(p));
        for(int i=0;i<enslaved_nodes(p).m;i++){
//            V(enslaved_nodes(p)(i))=TV(1,0,0);}}
            V(enslaved_nodes(p)(i))=plank_rigid_body.Pointwise_Object_Velocity(plank_rigid_body.World_Space_Point(positions_relative_to_plank_frames(p)(i)));}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    for(int p=0;p<num_planks;p++){
        //if(p==1) continue;
        for(int i=0;i<enslaved_nodes(p).m;i++){
            V(enslaved_nodes(p)(i))=TV();}}
}
//#####################################################################
// Function Get_Constrained_Particle_Data
//#####################################################################
void Get_Constrained_Particle_Data()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
 
    enslaved_nodes.Resize(num_planks);
    positions_relative_to_plank_frames.Resize(num_planks);

    for(int i=1;i<=mattress_grid.m*mattress_grid.n*mattress_grid.mn;i++){
        for(int p=0;p<num_planks;p++){
            RIGID_BODY<TV>& plank_rigid_body=*arb->rigid_body_list.rigid_bodies(plank_ids(p));
            if (!plank_rigid_body.Implicit_Geometry_Lazy_Outside(particles.X(i))) {
                enslaved_nodes(p).Append(i);
//                positions_relative_to_plank_frames(p).Append(plank_rigid_body.Rotation().Inverse_Rotate(particles.X(i)-plank_rigid_body.X()));}}}
                positions_relative_to_plank_frames(p).Append(plank_rigid_body.Frame().Inverse_Times(particles.X(i)));}}}
}
//#####################################################################
// Function Simple_Muscle_Across_Joint
//#####################################################################
void Simple_Muscle_Across_Joint()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solids_parameters.rigid_body_parameters.list;
    assert(rigid_body_particles.array_collection->Size()==0);

    add_ground=false;
    solids_parameters.gravity=0;

    int num_bodies=num_planks=parameter_list.Get_Parameter("num_bodies",(int)2);
    plank_ids.Resize(num_planks);
    bool use_bend_joint=parameter_list.Get_Parameter("use_bend_joint",false);
    TV plank_rescale=parameter_list.Get_Parameter("plank_rescale",TV(1,1,1));
    T target_angle=parameter_list.Get_Parameter("target_angle",-(T)pi/2);
    T initial_angle=parameter_list.Get_Parameter("initial_angle",0);
    T k_p=parameter_list.Get_Parameter("k_p",(T)100000);

    for(int i=0;i<num_bodies;i++){
        int id=plank_ids(i)=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/plank",(T).2);
        RIGID_BODY<TV>* rigid_body=&rigid_body_particles.Rigid_Body(id);
        rigid_body->simplicial_object->Rescale(plank_rescale.x,plank_rescale.y,plank_rescale.z);
        rigid_body->Rotation()=QUATERNION<T>(pi/2,TV(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(0.5);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name(STRING_UTILITIES::string_sprintf("body_%d",i));
        rigid_body->Set_Mass(1);

        if(i>1){
            JOINT<TV>* joint=0;
            if(use_bend_joint) joint=new ANGLE_JOINT<TV>();
            else joint=new POINT_JOINT<TV>();
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,0,1.05)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,.0,-1.05)));
            arb->Add_Articulation(i-1,i,joint);
            JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,&rigid_body_particles.Rigid_Body(i-1),&rigid_body_particles.Rigid_Body(i));
            joint->Set_Joint_Function(jfunc);
            jfunc->Set_Target_Angle(QUATERNION<T>(target_angle,TV(1,0,0)));
            jfunc->muscle_control=true;
            joint->joint_function->Set_k_p(k_p);
            joint->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(initial_angle,TV(1,0,0))));
        }
    }

    rigid_body_particles.frame(1).r=QUATERNION<T>(-pi/2,TV(0,0,1))*rigid_body_particles.frame(1).r;

    for(int i=2;i<=num_bodies;i++){
        std::string suffix=STRING_UTILITIES::string_sprintf("_j%d",(i-1));
        RIGID_BODY<TV> *parent=&rigid_body_particles.Rigid_Body(i-1),*child=&rigid_body_particles.Rigid_Body(i);

        muscle=Add_Basic_Muscle("flexor"+suffix,*parent,TV(0,0.11,-.8),*child,TV(0,0.11,-0.8),.01);

        //if(parameter_list.Get_Parameter("add_extensor",true)){
            // Add_Basic_Muscle("extensor"+suffix,*parent,TV(0,-0.05,0.5),*child,TV(0,-0.05,-0.5),0.1);
             // MAKE THESE RIGID BODY BINDINGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//             extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*parent,TV(0,-0.2,1)));
//             extensor->Add_Via_Point(new T_CONSTRAINED_POINT_IN_RIGID_BODY(*child,TV(0,-0.2,-1)));}
        //}

    }
    arb->Update_With_Breadth_First_Directed_Graph(1);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Write_Output_Files(frame);
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(STRING_UTILITIES::string_sprintf("%s/%d/muscle_info",output_directory.c_str(),frame));
    TYPED_OSTREAM typed_output(*output,stream_type);
    Write_Binary(typed_output,muscle_segment_particle_ids.m+tet_particle_ids.m);
    for(int i=0;i<muscle_segment_particle_ids.m;i++) Write_Binary(typed_output,solid_body_collection.deformable_body_collection.particles.X(muscle_segment_particle_ids(i)));
    for(int i=0;i<tet_particle_ids.m;i++) Write_Binary(typed_output,solid_body_collection.deformable_body_collection.particles.X(tet_particle_ids(i)));
    delete output;
}
};
}
#endif
