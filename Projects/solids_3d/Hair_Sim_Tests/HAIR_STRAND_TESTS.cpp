//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// HAIR_STRAND_TESTS
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_COMBINED.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fragments/PARTICLE_CONNECTIVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include "HAIR_STRAND_TESTS.h"

/* STRAND TESTS
 * 1) Twist example
 * 2) Push example
 * 3) Buckling example
 * 4) Stretch example
 * 5) Adhesion example
 * 6) Wind
 * 7) Wind stop
 * 8) Simple Test
 */

using namespace PhysBAM;
//#####################################################################
// Function HAIR_STRAND_TESTS
//#####################################################################
template<class T_input> HAIR_STRAND_TESTS<T_input>::
HAIR_STRAND_TESTS(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),use_adhesion(false),reset(false)
{
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Register_Options()
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-hairsim","","the hair sime to run");
    parse_args->Add_String_Argument("-modelname","","the rigid model to bind to");
    parse_args->Add_String_Argument("-guide","","the guide hair sim to read from");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Parse_Options()
{
    if(parse_args->Is_Value_Set("-d")) data_directory=parse_args->Get_String_Value("-d");
    sim_folder=parse_args->Get_String_Value("-hairsim");

    std::string parameter_file=(data_directory+"/"+sim_folder+"/"+parse_args->Get_String_Value("-params"));
    LOG::cout<<"PARAM FILE is "<<parameter_file<<std::endl;
    parameter_list.Begin_Parse(parameter_file);
    std::string test_name=parameter_list.Get_Parameter("test_name",(std::string)"Test");
    output_directory=STRING_UTILITIES::string_sprintf("Hair_Sim_Tests/%s_%d",test_name.c_str(),test_number);
    solids_parameters.cfl=parameter_list.Get_Parameter("cfl",(T)10);
    cfl_strain_rate=parameter_list.Get_Parameter("cfl_strain_rate",(T)0.1);
    solids_parameters.implicit_solve_parameters.cg_iterations=parameter_list.Get_Parameter("cg_iterations",(int)200);
    overdamping_fraction=parameter_list.Get_Parameter("overdamping_fraction",(T)40);
    restlength_clamp=parameter_list.Get_Parameter("restlength_clamp",(T)1e-4);
    edge_stiffness=parameter_list.Get_Parameter("edge_stiffness",(T)1e3);
    bending_stiffness=parameter_list.Get_Parameter("bending_stiffness",(T)1e3);
    torsion_stiffness=parameter_list.Get_Parameter("torsion_stiffness",(T)1e3);
    altitude_stiffness=parameter_list.Get_Parameter("altitude_stiffness",(T)1e3);
    ether_drag_wind=parameter_list.Get_Parameter("ether_drag",(T)20);
    momentum_conserving_projection_iterations=parameter_list.Get_Parameter("momentum_conserving_projection_iterations",(int)0);
    use_momentum_conserving_before=parameter_list.Get_Parameter("use_momentum_conserving_before",(bool)true);
    use_non_momentum_conserving_before=parameter_list.Get_Parameter("use_non_momentum_conserving_before",(bool)false);
    use_momentum_conserving_after=parameter_list.Get_Parameter("use_momentum_conserving_after",(bool)false);
    use_non_momentum_conserving_after=parameter_list.Get_Parameter("use_non_momentum_conserving_after",(bool)true);
    solids_parameters.triangle_collision_parameters.perform_self_collision=parameter_list.Get_Parameter("perform_self_collision",(bool)false);
    solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=parameter_list.Get_Parameter("self_collision_friction_coefficient",(T).4);
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=parameter_list.Get_Parameter("repulsion_thickness",(T)0.0005);
    use_implicit=parameter_list.Get_Parameter("use_implicit",(bool)false);
    parameter_list.End_Parse();
    LOG::cout<<"Parameters were"<<std::endl;
    parameter_list.Write(LOG::cout);
    parameter_list.Write(output_directory+"/parameters");
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Initialize_Bodies()
{
    //helper references
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    ARRAY<T> masses;

    momentum_conserving_projection_iterations=0;

    if (test_number==5){
        use_adhesion=true;
        adhesion_stiffness=(T)20.;
        adhesion_start_radius=(T).003;
        adhesion_stop_radius=(T).02;}
    
    last_frame=1000;
    frame_rate=30;

    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.enforce_repulsions_in_cg=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.verbose_dt=true;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-5;
    solids_parameters.triangle_collision_parameters.collisions_disable_repulsions_based_on_proximity_factor=(T)1.5;

    wind_start_time=(T)1;
    wind_stop_time=(T)5;

    //################################################################
    // Geometry Phase
    //################################################################
    SEGMENTED_CURVE<TV>& edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& fixed_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& extra_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& bending_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& torsion_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);

    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/masses",masses);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/particles",particles);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/edges.curve",edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/fixed_edges.curve",fixed_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/extra_edges.curve",extra_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/bending_edges.curve",bending_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/torsion_edges.curve",torsion_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/tets.gz",volume->mesh);
    //FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/project_mesh.gz",project_mesh);

    particles.Store_Velocity(true);
    
    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(edges,density,true);
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(extra_edges,density,true);

    fixed_edges.Update_Number_Nodes();
    edges.Update_Number_Nodes();
    extra_edges.Update_Number_Nodes();
    bending_edges.Update_Number_Nodes();
    torsion_edges.Update_Number_Nodes();
    assert(masses.m==particles.array_collection->Size());
    for(int i=0;i<masses.m;i++) {assert(masses(i));particles.mass(i)=masses(i);}
    //for(int i=1;i<=particles.array_collection->Size();i++) {particles.mass(i)=1.;}
    //for(int i=1;i<fixed_nodes.m+1;i++) particles.mass(fixed_nodes(i))=FLT_MAX;
    for(int i=1;i<fixed_nodes_start.m+1;i++) particles.mass(fixed_nodes_start(i))=FLT_MAX;
    for(int i=1;i<fixed_nodes_start.m+1;i++) particles.mass(fixed_nodes_end(i))=FLT_MAX;

    // Fix tetrahedra orientation
    int i=1;
    for(int t=0;t<volume->mesh.elements.m;t++){
        VECTOR<int,4>& nodes=volume->mesh.elements(t);
        if (TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) exchange(nodes[3],nodes[4]);
        i++;}
    for(int t=0;t<volume->mesh.elements.m;t++){
        VECTOR<int,4>& nodes=volume->mesh.elements(t);
        if (TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) PHYSBAM_FATAL_ERROR();}
    volume->Update_Number_Nodes();
    
    //################################################################
    // Mass Update Phase
    //################################################################
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    //################################################################
    // Forces
    //################################################################
    bool strain_limit=false,use_implicit=false;
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,true,true));

    LINEAR_SPRINGS<TV>* edge_springs=Create_Edge_Springs(edges,edge_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    edge_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(edge_springs);
    LINEAR_SPRINGS<TV>* extra_edge_springs=Create_Edge_Springs(extra_edges,edge_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    extra_edge_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(extra_edge_springs);
    LINEAR_SPRINGS<TV>* bending_springs=Create_Edge_Springs(bending_edges,bending_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    bending_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(bending_springs);
    LINEAR_SPRINGS<TV>* torsion_springs=Create_Edge_Springs(torsion_edges,torsion_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    torsion_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(torsion_springs);
    LINEAR_TET_SPRINGS<T> *tet_springs=Create_Tet_Springs(*volume,altitude_stiffness,overdamping_fraction,false,(T).1,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    tet_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(tet_springs);
    // drag forces
#if 0
    if(test_number>1 && test_number<8){
        WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,volume->Get_Boundary_Object());solid_body_collection.Add_Force(drag);
        drag->Use_Linear_Normal_Viscosity(1);drag->Use_Constant_Wind(0,TV());}
#endif
    //if(test_number==8){
        ETHER_DRAG<GRID<TV> >* ether_drag=new ETHER_DRAG<GRID<TV> >(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true);
        solid_body_collection.Add_Force(ether_drag);
        ether_drag->Use_Constant_Wind(ether_drag_wind);//}

    //hairs
    deformable_body_collection.deformable_geometry.Add_Structure(&edges);
    deformable_body_collection.deformable_geometry.Add_Structure(&fixed_edges);
    deformable_body_collection.deformable_geometry.Add_Structure(&extra_edges);
    deformable_body_collection.deformable_geometry.Add_Structure(&bending_edges);
    deformable_body_collection.deformable_geometry.Add_Structure(&torsion_edges);
    FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/%s/fixed_nodes_start",data_directory.c_str(),sim_folder.c_str()),fixed_nodes_start);
    FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/%s/fixed_nodes_end",data_directory.c_str(),sim_folder.c_str()),fixed_nodes_end);

    PHYSBAM_FATAL_ERROR();
#if 0
    SPARSE_UNION_FIND<> particle_connectivity(particles.array_collection->Size()+rigid_body_collection.rigid_body_particle.array_collection->Size());
    particle_to_spring_id.Resize(particles.array_collection->Size());
    edge_springs->Add_Fragment_Connectivity(particle_connectivity);extra_edge_springs->Add_Fragment_Connectivity(particle_connectivity);
    if(torsion_springs) torsion_springs->Add_Fragment_Connectivity(particle_connectivity);bending_springs->Add_Fragment_Connectivity(particle_connectivity);//guide_springs->Add_Fragment_Connectivity(union_find);}
    HAIR_ID next_segment_id(0);
    for(int p=1;p<=particles.array_collection->Size();p++){
        int root=particle_connectivity.Find(p);
        if(particle_to_spring_id(root)==HAIR_ID(0)){
            next_segment_id++;
            particle_to_spring_id(root)=next_segment_id;}
        particle_to_spring_id(p)=particle_to_spring_id(root);}
#endif
    
    if(use_adhesion){
        int max_connections=1;
        segment_adhesion=new SEGMENT_ADHESION<TV>(particles,edges.mesh,particle_to_spring_id,solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.intersecting_edge_edge_pairs);
        segment_adhesion->Set_Parameters(adhesion_stiffness,(T)10,adhesion_start_radius,adhesion_stop_radius,max_connections);
        solid_body_collection.Add_Force(segment_adhesion);
        segment_adhesion->Write_State(stream_type,output_directory+"/adhesion.0");
        segment_adhesion->Update_Partitions(false,solid_body_collection.deformable_body_collection.mpi_solids,output_directory);}

    for(int i=0;i<fixed_nodes_start.m;i++) init_positions_start.Append(particles.X(fixed_nodes_start(i)));
    for(int i=0;i<fixed_nodes_end.m;i++) init_positions_end.Append(particles.X(fixed_nodes_end(i)));

    solid_body_collection.deformable_body_collection.collisions.collision_structures.Append(&edges);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(&edges);

    // initial projection rest lengths (needs to happen after MPI restriction)
    project_restlengths.Resize(project_mesh.elements.m);
    for(int i=0;i<project_mesh.elements.m;i++){
        const VECTOR<int,2>& nodes=project_mesh.elements(i);
        project_restlengths(i)=(deformable_body_collection.particles.X(nodes[1])-deformable_body_collection.particles.X(nodes[2])).Magnitude();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    for(int i=0;i<fixed_nodes_start.m;i++) V(fixed_nodes_start(i))=TV();
    if(test_number!=4||velocity_time>(T)2.) for(int i=0;i<fixed_nodes_end.m;i++) V(fixed_nodes_end(i))=TV();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    INTERPOLATION_CURVE<T,FRAME<TV> > interp_start,interp_start2,interp_end;
    FRAME<TV> frame;
    T duration=(T)last_frame/(T)frame_rate;
    T length;
    T factor=2.;
    switch(test_number){
        case 1:
            length=init_positions_end(1)[1]-init_positions_start(init_positions_start.m)[1];
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point(duration,FRAME<TV>());
            if(velocity_time<1.) for(int i=0;i<fixed_nodes_end.m;i++) V(fixed_nodes_end(i))=TV((T)-.3*length,0,0);
            else if(velocity_time<2.) for(int i=0;i<fixed_nodes_end.m;i++) V(fixed_nodes_end(i))=TV();
            else{
                for(int i=0;i<fixed_nodes_end.m;i++){
                    TV rotation_axis=TV((T)1.,0,0);
                    V(fixed_nodes_end(i))=(T)pi*TV::Cross_Product(rotation_axis,particles.X(fixed_nodes_end(i))-TV(particles.X(fixed_nodes_end(i))(1),0,0));}}
            break;
        case 2:
            length=init_positions_end(1)[1]-init_positions_start(init_positions_start.m)[1];
            frame.t=TV((T)-.5*length,0,0);
            interp_end.Add_Control_Point(0,FRAME<TV>());
            interp_end.Add_Control_Point((T)1.,FRAME<TV>());
            interp_end.Add_Control_Point((T)2.,frame);
            interp_end.Add_Control_Point(duration,frame);
            frame.t=TV((T).5*length,0,0);
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point((T)1.,FRAME<TV>());
            interp_start.Add_Control_Point((T)2.,frame);
            interp_start.Add_Control_Point(duration,frame);
            break;
        case 3:
            interp_start.Add_Control_Point(0,FRAME<TV>());
            for(int i=1;i<=duration*factor;i++){
                frame.t=TV(0,(i%2==1)*(T)-1.*init_positions_start(1)[2],0);
                interp_start.Add_Control_Point((T)i/factor,frame);
            }
            break;
        case 4:
            length=init_positions_end(1)[2]-init_positions_start(init_positions_start.m)[2];
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point(duration,FRAME<TV>());
            if(velocity_time<(T)3.){
                frame.t=TV(0,(T).5*length,0);
                interp_end.Add_Control_Point(0,FRAME<TV>());
                interp_end.Add_Control_Point((T)2.,FRAME<TV>());
                interp_end.Add_Control_Point((T)4.,frame);
                interp_end.Add_Control_Point(duration,frame);}
            else if(!reset) {fixed_nodes_end.Remove_All();reset=true;}
            break;
        case 5:
            length=init_positions_end(1)[1]-init_positions_start(init_positions_start.m)[1];
            frame.t=TV(2*length,0,0);
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point((T)3.,FRAME<TV>());
            interp_start.Add_Control_Point((T)6.,frame);
            interp_start.Add_Control_Point(duration,frame);
            interp_end.Add_Control_Point(0,FRAME<TV>());
            interp_end.Add_Control_Point(duration,FRAME<TV>());
            if(velocity_time>2.&&!reset){reset=true;fixed_nodes_end.Remove_End();}
            interp_start2.Add_Control_Point((T)0,FRAME<TV>());
            frame.t=TV(length*(T).01,0,0);
            interp_start2.Add_Control_Point((T)1.,frame);
            interp_start2.Add_Control_Point(duration,frame);
            break;
        case 8:
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point((T)1,FRAME<TV>());
            interp_start.Add_Control_Point((T)2,FRAME<TV>(TV(1,0,0)));
            interp_start.Add_Control_Point((T)3,FRAME<TV>(TV(0,0,0)));
            interp_start.Add_Control_Point((T)4,FRAME<TV>(TV(0,0,0)));
            interp_start.Add_Control_Point((T)5,FRAME<TV>(TV(0,1,0)));
            interp_start.Add_Control_Point((T)6,FRAME<TV>(TV(1,1,0)));
            interp_start.Add_Control_Point((T)7,FRAME<TV>(TV(0,0,0)));
            break;
        default:
            break;}
    if(test_number==5) for(int i=0;i<fixed_nodes_start.m;i++){if(i==fixed_nodes_start.m) V(fixed_nodes_start(i))=interp_start.Derivative(velocity_time).linear; else V(fixed_nodes_start(i))=interp_start2.Derivative(velocity_time).linear;}
    else for(int i=0;i<fixed_nodes_start.m;i++) V(fixed_nodes_start(i))=interp_start.Derivative(velocity_time).linear;
    if(test_number!=1&&(test_number!=4||velocity_time>(T)2.)) for(int i=0;i<fixed_nodes_end.m;i++) V(fixed_nodes_end(i))=interp_end.Derivative(velocity_time).linear;
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    INTERPOLATION_CURVE<T,FRAME<TV> > interp_start,interp_end;
    INTERPOLATION_CURVE<T,T> interp_angle;
    FRAME<TV> frame;
    T duration=(T)last_frame/(T)frame_rate;
    T length;
    return;
    switch(test_number){
        case 1:
            length=init_positions_end(1)[1]-init_positions_start(init_positions_start.m)[1];
            interp_angle.Add_Control_Point(0,0);
            interp_angle.Add_Control_Point((T)1.,0);
            interp_angle.Add_Control_Point(duration,duration*(T)pi);
            if (time<(T)2.){
                frame.t=TV((T)-.3*length,0,0);
                interp_end.Add_Control_Point(0,FRAME<TV>());
                interp_end.Add_Control_Point((T)1.,frame);
                interp_end.Add_Control_Point((T)2.,frame);}
            else{
                if(!reset){
                    for(int i=0;i<fixed_nodes_end.m;i++) init_positions_end(i)+=TV((T)-.3*length,0,0);
                    reset=true;}
                frame.r=ROTATION<TV>(interp_angle.Value(time),TV((T)1.,0,0));
                interp_end.Add_Control_Point(0,frame);
                interp_end.Add_Control_Point(duration,frame);}
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point(duration,FRAME<TV>());
            break;
        case 2:
            length=init_positions_end(1)[1]-init_positions_start(init_positions_start.m)[1];
            frame.t=TV((T)-.5*length,0,0);
            interp_end.Add_Control_Point(0,FRAME<TV>());
            interp_end.Add_Control_Point((T)1.,FRAME<TV>());
            interp_end.Add_Control_Point((T)2.,frame);
            interp_end.Add_Control_Point(duration,frame);
            frame.t=TV((T).5*length,0,0);
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point((T)1.,FRAME<TV>());
            interp_start.Add_Control_Point((T)2.,frame);
            interp_start.Add_Control_Point(duration,frame);
            break;
        case 3:
            interp_start.Add_Control_Point(0,FRAME<TV>());
            for(int i=0;i<=duration;i++){
                frame.t=TV(0,i%2==0*init_positions_start(0)[0],0);
                interp_start.Add_Control_Point((T)i,frame);
            }
            break;
        case 4:
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point(duration,FRAME<TV>());
            if(time<(T)2.){
                frame.t=TV(0,(T).5,0);
                interp_end.Add_Control_Point(0,FRAME<TV>());
                interp_end.Add_Control_Point(duration,frame);}
            else if(!reset) {fixed_nodes_end.Remove_All(); reset=true;}
            break;
        case 5:
            break;
        case 8:
            interp_start.Add_Control_Point(0,FRAME<TV>());
            interp_start.Add_Control_Point((T)1,FRAME<TV>());
            interp_start.Add_Control_Point((T)2,FRAME<TV>(TV(1,0,0)));
            interp_start.Add_Control_Point((T)3,FRAME<TV>(TV(0,0,0)));
            interp_start.Add_Control_Point((T)4,FRAME<TV>(TV(0,0,0)));
            interp_start.Add_Control_Point((T)5,FRAME<TV>(TV(0,1,0)));
            interp_start.Add_Control_Point((T)6,FRAME<TV>(TV(1,1,0)));
            interp_start.Add_Control_Point((T)7,FRAME<TV>(TV(0,0,0)));
            break;
        default:
            break;}
    if(test_number==5) for(int i=0;i<fixed_nodes_start.m;i++){if(i>fixed_nodes_start.m/2) X(fixed_nodes_start(i))=interp_start.Value(time)*init_positions_start(i); else X(fixed_nodes_start(i))=init_positions_start(i);}
    else for(int i=0;i<fixed_nodes_start.m;i++) X(fixed_nodes_start(i))=interp_start.Value(time)*init_positions_start(i);
    for(int i=0;i<fixed_nodes_end.m;i++) X(fixed_nodes_end(i))=interp_end.Value(time)*init_positions_end(i);
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Update_Time_Varying_Material_Properties(const T time)
{
    return;
    //start wind
    if(test_number>=6 && test_number<8){
        if(time>wind_start_time){
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(17);drag.Use_Constant_Wind(0,TV(0,0,(T)-1.));}}
    //stop wind
    if(test_number==7){
        if(time>wind_stop_time){
            INTERPOLATION_CURVE<T,TV> wind_stop;
            wind_stop.Add_Control_Point(wind_stop_time-(T).3,TV(0,0,(T)-.25));
            wind_stop.Add_Control_Point(wind_stop_time,TV());
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(0);drag.Use_Constant_Wind(0,TV());}}
}
//#####################################################################
// Function Add_External_Impulses_Before
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt)
{
    Add_External_Impulses_Helper(V,time,dt,use_momentum_conserving_before,use_non_momentum_conserving_before);
}
//#####################################################################
// Function Add_External_Impulses_Before
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt)
{
    Add_External_Impulses_Helper(V,time,dt,use_momentum_conserving_after,use_non_momentum_conserving_after);
}
//#####################################################################
// Function Add_External_Impulses_Helper
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Add_External_Impulses_Helper(ARRAY_VIEW<TV> V,const T time,const T dt,bool use_momentum_conserving,bool use_non_momentum_conserving)
{
    //if(time<start_time) return;
    return;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    if(write_substeps_level>-1){
        ARRAY<TV> positions_save(particles.X);
        for(int i=1;i<=particles.array_collection->Size();i++) particles.X(i)=particles.X(i)+dt*V(i);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before impulse",2,2);
        particles.X=positions_save;}
    
    if(use_momentum_conserving){
        for(int iteration=0;iteration<momentum_conserving_projection_iterations;iteration++){
            for(int i=0;i<project_mesh.elements.m;i++){
                int child,parent;project_mesh.elements(i).Get(child,parent);
                T restlength=project_restlengths(i);
                TV X_child_new=particles.X(child)+dt*V(child);
                TV X_parent_new=particles.X(parent)+dt*V(parent);
                TV direction=(X_child_new-X_parent_new);T current_length=direction.Normalize();
                TV impulse=Pseudo_Divide((current_length-restlength)*direction/dt,particles.one_over_effective_mass(child)+particles.one_over_effective_mass(parent));
                V(parent)+=impulse*particles.one_over_mass(parent);
                V(child)-=impulse*particles.one_over_mass(child);}
            //T max_error=0;
            //for(int i=0;i<project_mesh.elements.m;i++){
            //    int child,parent;project_mesh.elements(i).Get(child,parent);
            //    T restlength=project_restlengths(i);
            //    TV X_child_new=particles.X(child)+dt*V(child);
            //    TV X_parent_new=particles.X(parent)+dt*V(parent);
            //    TV direction=(X_child_new-X_parent_new);T current_length=direction.Normalize();
            //    max_error=max(abs(current_length-restlength),max_error);}
            //LOG::cout<<"iteration "<<iteration<<" max error "<<max_error<<std::endl;
        }
    }

    if(use_non_momentum_conserving){
        for(int i=0;i<project_mesh.elements.m;i++){
            int child,parent;project_mesh.elements(i).Get(child,parent);
            T restlength=project_restlengths(i);
            TV X_child_new=particles.X(child)+dt*V(child);
            TV X_parent_new=particles.X(parent)+dt*V(parent);
            TV direction=(X_child_new-X_parent_new);T current_length=direction.Normalize();
            TV X_child_new_prime=X_child_new;
            if(current_length>restlength) X_child_new_prime=X_parent_new+restlength*direction;
            V(child)=(X_child_new_prime-particles.X(child))/dt;}}

    if(write_substeps_level>-1){
        ARRAY<TV> positions_save(particles.X);
        for(int i=1;i<=particles.array_collection->Size();i++) particles.X(i)=particles.X(i)+dt*V(i);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after impulse",2,2);
        particles.X=positions_save;}
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Preprocess_Solids_Substep(const T time,const int substep) {
    if(segment_adhesion) segment_adhesion->Update_Springs(true);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Postprocess_Frame(const int frame)
{
    if(segment_adhesion) segment_adhesion->Write_State(stream_type,output_directory+STRING_UTILITIES::string_sprintf("/adhesion.%d",frame));
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_input> void HAIR_STRAND_TESTS<T_input>::
Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    if(segment_adhesion) segment_adhesion->Write_State(stream_type,output_directory+STRING_UTILITIES::string_sprintf("/adhesion.%d",frame));
}
//#####################################################################
template class HAIR_STRAND_TESTS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HAIR_STRAND_TESTS<double>;
#endif
