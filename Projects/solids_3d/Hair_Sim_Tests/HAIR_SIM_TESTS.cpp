//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// HAIR_SIM_TESTS
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_COMBINED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_COMBINED_EULERIAN.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fragments/PARTICLE_CONNECTIVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include "HAIR_SIM_TESTS.h"

template<class T> int round_number(const T x)
{
    return (int)(((x-(int)x)>=0.5)?(int)x+1:(int)x);
}

using namespace PhysBAM;
//#####################################################################
// Function HAIR_SIM_TESTS
//#####################################################################
template<class T_input> HAIR_SIM_TESTS<T_input>::
HAIR_SIM_TESTS(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),start_time(2.),tests(*this,solid_body_collection),segment_adhesion(0),guide_adhesion(0),guide_object1(0),guide_object2(0),current_levelset(0)
{
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
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
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Parse_Options()
{
    BASE::Parse_Options();
    if(parse_args->Is_Value_Set("-d")) data_directory=parse_args->Get_String_Value("-d");
    sim_folder=parse_args->Get_String_Value("-hairsim");
    rigid_model=parse_args->Get_String_Value("-modelname");
    guide_sim_folder=parse_args->Get_String_Value("-guide");

    std::string parameter_file=(data_directory+"/"+sim_folder+"/"+parse_args->Get_String_Value("-params"));
    if(!FILE_UTILITIES::File_Exists(parameter_file)){
        LOG::cerr<<"Parameter file '"<<parameter_file<<"' does not exist"<<std::endl;
        exit(1);}
    LOG::cout<<"PARAM FILE is "<<parameter_file<<std::endl;
    parameter_list.Begin_Parse(parameter_file);
    std::string test_name=parameter_list.Get_Parameter("test_name",(std::string)"Test");
    output_directory=STRING_UTILITIES::string_sprintf("Hair_Sim_Tests/%s_%d",test_name.c_str(),test_number);
    use_wind=parameter_list.Get_Parameter("use_wind",(bool)false);
    use_drag=parameter_list.Get_Parameter("use_drag",(bool)false);
    drag_viscosity=parameter_list.Get_Parameter("drag_viscosity",(T)5);
    solids_parameters.cfl=parameter_list.Get_Parameter("cfl",(T)4);
    cfl_strain_rate=parameter_list.Get_Parameter("cfl_strain_rate",(T)0.1);
    solids_parameters.implicit_solve_parameters.cg_iterations=parameter_list.Get_Parameter("cg_iterations",(int)200);
    solids_parameters.implicit_solve_parameters.cg_tolerance=parameter_list.Get_Parameter("cg_tolerance",(T)5e-2);
    overdamping_fraction=parameter_list.Get_Parameter("overdamping_fraction",(T)40);
    restlength_clamp=parameter_list.Get_Parameter("restlength_clamp",(T)1e-4);
    edge_stiffness=parameter_list.Get_Parameter("edge_stiffness",(T)1e4);
    bending_stiffness=parameter_list.Get_Parameter("bending_stiffness",(T)1e3);
    torsion_stiffness=parameter_list.Get_Parameter("torsion_stiffness",(T)1e3);
    altitude_stiffness=parameter_list.Get_Parameter("altitude_stiffness",(T)1e3);
    use_adhesion=parameter_list.Get_Parameter("use_adhesion",(bool)false);
    adhesion_stiffness=parameter_list.Get_Parameter("adhesion_stiffness",(T).001);
    adhesion_start_radius=parameter_list.Get_Parameter("adhesion_start_radius",(T).002);
    adhesion_stop_radius=parameter_list.Get_Parameter("adhesion_stop_radius",(T).005);
    max_connections=parameter_list.Get_Parameter("adhesion_connections",(int)5);
    solids_parameters.triangle_collision_parameters.perform_self_collision=parameter_list.Get_Parameter("perform_self_collision",(bool)false);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.allow_intersections=parameter_list.Get_Parameter("allow_intersections",(bool)false);
    solids_parameters.triangle_collision_parameters.turn_off_all_collisions=parameter_list.Get_Parameter("turn_off_all_collisions",(bool)false);
    solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=parameter_list.Get_Parameter("self_collision_friction_coefficient",(T).4);
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=parameter_list.Get_Parameter("repulsion_thickness",(T)1e-3);
    solids_parameters.enforce_repulsions_in_cg=parameter_list.Get_Parameter("enforce_repulsions_in_cg",(bool)false);
    use_spring_guide=parameter_list.Get_Parameter("use_spring_guide",(bool)false);
    use_guide=parameter_list.Get_Parameter("use_guide",(bool)false);
    guide_edge_stiffness=parameter_list.Get_Parameter("guide_edge_stiffness",(T)1e5);
    guide_altitude_stiffness=parameter_list.Get_Parameter("guide_altitude_stiffness",(T)1e4);
    guide_stiffness=parameter_list.Get_Parameter("guide_stiffness",(T)1e5);
    max_connections=parameter_list.Get_Parameter("guide_max_connections",(int)100);
    guide_thickness=parameter_list.Get_Parameter("guide_thickness",(T)10);
    use_deforming_levelsets=parameter_list.Get_Parameter("use_deforming_levelsets",(bool)false);
    use_collisions_mass_modify=parameter_list.Get_Parameter("use_collisions_mass_modify",(bool)true);
    momentum_conserving_projection_iterations=parameter_list.Get_Parameter("momentum_conserving_projection_iterations",(int)3);
    use_momentum_conserving_before=parameter_list.Get_Parameter("use_momentum_conserving_before",(bool)true);
    use_non_momentum_conserving_before=parameter_list.Get_Parameter("use_non_momentum_conserving_before",(bool)false);
    use_momentum_conserving_after=parameter_list.Get_Parameter("use_momentum_conserving_after",(bool)false);
    use_non_momentum_conserving_after=parameter_list.Get_Parameter("use_non_momentum_conserving_after",(bool)true);
    use_progressive_collision_thickness=parameter_list.Get_Parameter("use_progressive_collision_thickness",(bool)true);
    project_second_node_to_surface=parameter_list.Get_Parameter("project_second_node_to_surface",(bool)false);
    head_friction=parameter_list.Get_Parameter("head_friction",(T)0);
    use_eulerian_level_set_interpolation=parameter_list.Get_Parameter("use_eulerian_level_set_interpolation",(bool)false);
    frame_rate=parameter_list.Get_Parameter("frame_rate",(T)30);

    T levelset_frame_rate=parameter_list.Get_Parameter("levelset_frame_rate",(T)60); // twice the sampling rate as 
    use_implicit=parameter_list.Get_Parameter("use_implicit",(bool)false);
    cameras=parameter_list.Get_Parameter("cameras",std::string(""));
    parameter_list.End_Parse();
    LOG::cout<<"Parameters were"<<std::endl;
    parameter_list.Write(LOG::cout);
    parameter_list.Write(output_directory+"/parameters");

    // frame rates and level set sampling rate
    levelset_frequency=(T)1/levelset_frame_rate;
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Limit_Solids_Dt(T& dt,const T time)
{
    if(use_deforming_levelsets) Update_Keyframed_Parameters_For_Time_Update(time);
    T max_binding_speed=sqrt(max_binding_speed_squared);
    T max_displacement=dt*max_binding_speed;
    const T max_displacement_allowed=(T).01;
    T old_dt=dt;
    if(max_displacement>max_displacement_allowed) dt=max_displacement_allowed/max_binding_speed;
    LOG::cout<<"Clamping to allowed displacement of "<<max_displacement_allowed<<" displacement was "<<max_displacement<<" with dt="<<old_dt<<" and now dt="<<dt<<std::endl;
}
//#####################################################################
// Function Compute_Binding_Velocities
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Compute_Binding_Velocities()
{
    binding_velocities.Resize(bindings1_to_apply.m);
    max_binding_speed_squared=0;
    for(int i=0;i<binding_velocities.m;i++){
        binding_velocities(i)=(bindings2_to_apply(i)-bindings1_to_apply(i))*levelset_frequency;
        max_binding_speed_squared=max(binding_velocities(i).Magnitude_Squared(),max_binding_speed_squared);}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Initialize_Bodies()
{
    //helper references
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    ARRAY<T> masses;
    
    // compute last frame
    if(use_deforming_levelsets){
        int levelset_frames;FILE_UTILITIES::Read_From_Text_File(STRING_UTILITIES::string_sprintf("%s/%s/motion/last_frame",data_directory.c_str(),sim_folder.c_str()),levelset_frames);
        last_frame=(int)((T)levelset_frames*levelset_frequency*(T)frame_rate)+round_number(start_time*frame_rate);
        last_frame-=1;}
    else{
        last_frame=1000;}

    // give mon hints
    LOG::cout<<"MONITOR begin_frame="<<this->first_frame<<std::endl;
    LOG::cout<<"MONITOR output_directory="<<(FILE_UTILITIES::Get_Working_Directory()+"/"+output_directory)<<std::endl;
    LOG::cout<<"MONITOR end_frame="<<last_frame<<std::endl;
    if(cameras!=""){
        LOG::cout<<"MONITOR gl="<<cameras<<std::endl;
        LOG::cout<<"MONITOR glkeys=("<<std::endl;}
    
    solids_parameters.triangle_collision_parameters.temporary_enable_collisions=true;
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.mass_modifier=use_collisions_mass_modify?this:0;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.triangle_collision_parameters.total_collision_loops=1;
    //solids_parameters.maximum_levelset_collision_projection_velocity=1.;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=30;
    solids_parameters.rigid_body_collision_parameters.contact_project_iterations=0;
    solids_parameters.rigid_body_collision_parameters.contact_iterations=1;
    solids_parameters.use_rigid_deformable_contact=false;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-5;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).9; // used for torus and other volumetrics
    solids_parameters.triangle_collision_parameters.collisions_final_repulsion_youngs_modulus=(T)30;
    solids_parameters.triangle_collision_parameters.repulsions_youngs_modulus=(T)30;
    solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=4;
    solid_body_collection.deformable_body_collection.triangle_collisions.Set_Attempts_For_Rigid_Collisions(true,0);
    solids_parameters.triangle_collision_parameters.check_initial_mesh_for_self_intersection=false;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.triangle_collision_parameters.collisions_disable_repulsions_based_on_proximity_factor=(T)1.5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.verbose_dt=true;
    solids_parameters.deformable_object_collision_parameters.disable_multiple_levelset_collisions=false;
    solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
    solids_parameters.use_post_cg_constraints=false;
    //solid_body_collection.print_residuals=true;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    //data_directory=parse_args.Get_String_Value("-directory");
    assert(!rigid_model.empty());assert(!data_directory.empty());assert(!sim_folder.empty());
    wind_start_time=(T)1;
    wind_stop_time=(T)5;

    if(restart) current_levelset=max(0,round_number(((float)restart_frame/(float)frame_rate-start_time)/levelset_frequency));

    //################################################################
    // Geometry Phase
    //################################################################
    ARRAY<int> fixed_nodes;
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/fixed_nodes",fixed_nodes);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/masses",masses);
    SEGMENTED_CURVE<TV>& edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& fixed_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& extra_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& bending_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& torsion_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& sim_guide_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& sim_guide_tet_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& guide_tet_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    guide_volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    sim_guide_volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);

    if(use_guide){
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/guide_edges.curve",guide_tet_edges.mesh);
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/guide_tets.gz",guide_volume->mesh);}
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/particles",particles);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/edges.curve",edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/fixed_edges.curve",fixed_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/extra_edges.curve",extra_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/bending_edges.curve",bending_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/torsion_edges.curve",torsion_edges.mesh);
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/tets.gz",volume->mesh);
    
    FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/project_mesh.gz",project_mesh);

    particles.Store_Velocity(true);
    for(int i=0;i<particles.array_collection->Size();i++) {active_particles.Append(i);}
    
    //fake deformable object for guide hairs
    offset=particles.array_collection->Size();
    if(!guide_sim_folder.empty()) {
        COLLISION_GEOMETRY_COLLECTION<TV> guide_list;
        guide_object1=new DEFORMABLE_BODY_COLLECTION<TV>(solid_body_collection.example_forces_and_velocities,guide_list);
        guide_object2=new DEFORMABLE_BODY_COLLECTION<TV>(solid_body_collection.example_forces_and_velocities,guide_list);
        guide_object1->Read(stream_type,guide_sim_folder+"/",0,-1,1,solids_parameters.write_from_every_process);
        guide_object2->Read(stream_type,guide_sim_folder+"/",1,-1,0,solids_parameters.write_from_every_process);
        SEGMENTED_CURVE<TV>& guide_edges=guide_object1->deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>(0);
        particles.array_collection->Add_Elements(guide_edges.particles.array_collection->Size());
        for(int i=0;i<guide_edges.mesh.elements.m;i++){sim_guide_edges.mesh.elements.Append(guide_edges.mesh.elements(i)+VECTOR<int,2>(offset,offset));}
        for(int i=0;i<guide_edges.particles.array_collection->Size();i++){particles.X(offset+i)=guide_edges.particles.X(i);particles.V(offset+i)=static_cast<DEFORMABLE_PARTICLES<TV>&>(guide_edges.particles).V(i);}
    }

    if(use_guide){
        for(int i=0;i<guide_tet_edges.mesh.elements.m;i++){
            if(i%6<=3&&i%6>0) sim_guide_tet_edges.mesh.elements.Append(guide_tet_edges.mesh.elements(i)+VECTOR<int,2>(0,offset));
            else sim_guide_tet_edges.mesh.elements.Append(guide_tet_edges.mesh.elements(i)+VECTOR<int,2>(offset,offset));}
        for(int i=1;i<guide_volume->mesh.elements.m;i++) sim_guide_volume->mesh.elements.Append(guide_volume->mesh.elements(i)+VECTOR<int,4>(0,offset,offset,offset));}

    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(edges,density,true);
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(extra_edges,density,true);

    fixed_edges.Update_Number_Nodes();
    edges.Update_Number_Nodes();
    extra_edges.Update_Number_Nodes();
    bending_edges.Update_Number_Nodes();
    torsion_edges.Update_Number_Nodes();
    sim_guide_edges.Update_Number_Nodes();
    sim_guide_tet_edges.Update_Number_Nodes();
    assert(masses.m==offset);
    for(int i=0;i<masses.m;i++) {assert(masses(i));particles.mass(i)=masses(i);}
    //for(int i=0;i<masses.m;i++) {assert(masses(i));particles.mass(i)=4.;}
    for(int i=offset+1;i<=particles.array_collection->Size();i++) particles.mass(i)=FLT_MAX;
    for(int i=1;i<fixed_nodes.m+1;i++) particles.mass(fixed_nodes(i))=FLT_MAX;

    // Find edge distances
    if(use_collisions_mass_modify || use_progressive_collision_thickness){
        LOG::cout<<"computing distances..."<<std::endl;
        distance_to_root.Resize(particles.array_collection->Size());
        distance_to_root.Fill(0);
        distance_to_root.Subset(fixed_nodes).Fill(1);
        edges.mesh.Initialize_Neighbor_Nodes();
        while(1){
            int marked=0;
            for(int p=0;p<distance_to_root.m;p++) if(distance_to_root(p)==0){
                ARRAY<int>& neighbor_nodes=(*edges.mesh.neighbor_nodes)(p);
                for(int j=0;j<neighbor_nodes.m;j++){
                    int neighbor_p=neighbor_nodes(j);
                    if(distance_to_root(neighbor_p)!=0){
                        distance_to_root(p)=distance_to_root(neighbor_p)+1;marked++;break;}}}
            if(!marked) break;}
        if(use_progressive_collision_thickness){
            collision_tolerances.Resize(particles.array_collection->Size());
            for(int p=0;p<distance_to_root.m;p++){
                collision_tolerances(p)=(T)2e-4*min((T)5,(T)distance_to_root(p));} // only do for the first 10 segments
            LOG::cout<<collision_tolerances<<std::endl;}
        solid_body_collection.deformable_body_collection.collisions.collision_tolerances=&collision_tolerances;}


    // Fix tetrahedra orientation
    int i=1;
    for(int t=0;t<volume->mesh.elements.m;t++){
        VECTOR<int,4>& nodes=volume->mesh.elements(t);
        if (TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) exchange(nodes[2],nodes[3]);
        i++;}
    for(int t=0;t<volume->mesh.elements.m;t++){
        VECTOR<int,4>& nodes=volume->mesh.elements(t);
        if (TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) PHYSBAM_FATAL_ERROR();}
    volume->Update_Number_Nodes();
    for(int t=0;t<sim_guide_volume->mesh.elements.m;t++){
        VECTOR<int,4>& nodes=sim_guide_volume->mesh.elements(t);
        if (TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) exchange(nodes[2],nodes[3]);
        i++;}
    for(int t=0;t<sim_guide_volume->mesh.elements.m;t++){
        VECTOR<int,4>& nodes=sim_guide_volume->mesh.elements(t);
        if (TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) PHYSBAM_FATAL_ERROR();}
    sim_guide_volume->Update_Number_Nodes();
    
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
    bool strain_limit=false;
    ARRAY<int> tet_node_list;
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,&active_particles,NULL));
    PHYSBAM_DEBUG_PRINT("Parameters",overdamping_fraction,edge_stiffness);
    //LINEAR_SPRINGS<TV>* guide_springs=Create_Edge_Springs(sim_guide_edges,edge_stiffness,overdamping_fraction);
    LINEAR_SPRINGS<TV>* edge_springs=Create_Edge_Springs(edges,edge_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    edge_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(edge_springs);
    LINEAR_SPRINGS<TV>* extra_edge_springs=Create_Edge_Springs(extra_edges,edge_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    extra_edge_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(extra_edge_springs);
    LINEAR_SPRINGS<TV>* bending_springs=Create_Edge_Springs(bending_edges,bending_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);

    LINEAR_SPRINGS<TV>* torsion_springs=0;
    bending_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(bending_springs);
    torsion_springs=Create_Edge_Springs(torsion_edges,torsion_stiffness,overdamping_fraction,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    torsion_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(torsion_springs);
    LINEAR_TET_SPRINGS<T> *tet_springs=Create_Tet_Springs(*volume,altitude_stiffness,overdamping_fraction,false,(T).1,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
    tet_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(tet_springs);
    if(use_guide){
        LINEAR_SPRINGS<TV>* guide_edge_springs=Create_Edge_Springs(sim_guide_tet_edges,guide_edge_stiffness,overdamping_fraction);
        extra_edge_springs->Clamp_Restlength(restlength_clamp);
        solid_body_collection.Add_Force(guide_edge_springs);
        LINEAR_TET_SPRINGS<T> *guide_tet_springs=Create_Tet_Springs(*sim_guide_volume,guide_altitude_stiffness,overdamping_fraction,false,(T).1,true,(T).1,true,(T)0,true);
        guide_tet_springs->Clamp_Restlength(restlength_clamp);
        solid_body_collection.Add_Force(guide_tet_springs);}
    SPARSE_UNION_FIND<> particle_connectivity(particles.array_collection->Size()+rigid_body_collection.rigid_body_particle.array_collection->Size());
    HAIR_ID next_segment_id(0);
    PHYSBAM_FATAL_ERROR();
#if 0
    particle_to_spring_id.Resize(particles.array_collection->Size());
    edge_springs->Add_Fragment_Connectivity(particle_connectivity);extra_edge_springs->Add_Fragment_Connectivity(particle_connectivity);
    if(torsion_springs) torsion_springs->Add_Fragment_Connectivity(particle_connectivity);bending_springs->Add_Fragment_Connectivity(particle_connectivity);//guide_springs->Add_Fragment_Connectivity(union_find);}
    for(int p=0;p<particles.array_collection->Size();p++){
        int root=particle_connectivity.Find(p);
        if(particle_to_spring_id(root)==HAIR_ID(0)){
            next_segment_id++;
            particle_to_spring_id(root)=next_segment_id;}
        particle_to_spring_id(p)=particle_to_spring_id(root);}
    number_of_hairs=next_segment_id;
#endif
    spring_id_to_particle.Resize(number_of_hairs);
    for(int i=0;i<fixed_nodes.m;i++) spring_id_to_particle(particle_to_spring_id(fixed_nodes(i)))=fixed_nodes(i);
    // adhesion forces
    if(use_adhesion){
        segment_adhesion=new SEGMENT_ADHESION<TV>(particles,edges.mesh,particle_to_spring_id,solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.intersecting_edge_edge_pairs);
        segment_adhesion->Set_Parameters(adhesion_stiffness,(T)10,adhesion_start_radius,adhesion_stop_radius,max_connections);
        solid_body_collection.Add_Force(segment_adhesion);
        if(restart) segment_adhesion->Read_State(stream_type,output_directory+STRING_UTILITIES::string_sprintf("/adhesion.%d",restart_frame));
        else segment_adhesion->Write_State(stream_type,output_directory+"/adhesion.0");}
    // drag forces
    ETHER_DRAG<GRID<TV> >* ether_drag=new ETHER_DRAG<GRID<TV> >(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,&active_particles,NULL);
    solid_body_collection.Add_Force(ether_drag);
    ether_drag->Use_Constant_Wind(20);
    PHYSBAM_ASSERT(!use_wind);
    //WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,volume->Get_Boundary_Object());
    //solid_body_collection.Add_Force(drag);
    //if(use_wind) {drag->Use_Linear_Normal_Viscosity(17);drag->Use_Constant_Wind(0,TV(0,(T).5,(T)-1.));}
    //else {drag->Use_Linear_Normal_Viscosity(1);drag->Use_Constant_Wind(0,TV());}

    //################################################################
    // Collision Bodies
    //################################################################
    //rigid bodies
    //deforming level sets
    if(use_deforming_levelsets){
        int id=rigid_body_collection.Add_Rigid_Body(stream_type,"",(T)1,false,false,false,false);
        implicit_rigid_body=&rigid_body_collection.Rigid_Body(id);
        implicit_rigid_body->is_static=true;
        //implicit_rigid_body->Set_Coefficient_Of_Friction((T)0);}
        implicit_rigid_body->Set_Coefficient_Of_Friction((T)head_friction);

        GRID<TV>& grid_1=*new GRID<TV>;
        ARRAY<T,VECTOR<int,3> >& phi_1=*new ARRAY<T,VECTOR<int,3> >;
        ARRAY<TV,VECTOR<int,3> >& velocity_1=*new ARRAY<TV,VECTOR<int,3> >;
        LEVELSET_3D<GRID<TV> > levelset_1(grid_1,phi_1);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s/motion/out.%d.phi",data_directory.c_str(),sim_folder.c_str(),current_levelset),levelset_1);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s/motion/velocities.%d",data_directory.c_str(),sim_folder.c_str(),current_levelset),velocity_1);
        LOG::cout<<"READING LEVELSET NUMBER: "<<current_levelset<<std::endl;
        LEVELSET_IMPLICIT_OBJECT<TV> *levelset_implicit_object_1=new LEVELSET_IMPLICIT_OBJECT<TV>(levelset_1.grid,levelset_1.phi);
        levelset_implicit_object_1->V=&velocity_1;
        GRID<TV> &grid_2=*new GRID<TV>;
        ARRAY<T,VECTOR<int,3> > &phi_2=*new ARRAY<T,VECTOR<int,3> >;
        ARRAY<TV,VECTOR<int,3> >& velocity_2=*new ARRAY<TV,VECTOR<int,3> >;
        LEVELSET_3D<GRID<TV> > levelset_2(grid_2,phi_2);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s/motion/out.%d.phi",data_directory.c_str(),sim_folder.c_str(),++current_levelset),levelset_2);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s/motion/velocities.%d",data_directory.c_str(),sim_folder.c_str(),current_levelset),velocity_2);
        LOG::cout<<"READING LEVELSET NUMBER: "<<current_levelset<<std::endl;
        LEVELSET_IMPLICIT_OBJECT<TV> *levelset_implicit_object_2=new LEVELSET_IMPLICIT_OBJECT<TV>(levelset_2.grid,levelset_2.phi);
        levelset_implicit_object_2->V=&velocity_2;
        if(use_eulerian_level_set_interpolation){
            IMPLICIT_OBJECT_COMBINED_EULERIAN<TV> &combined_eulerian=*new IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>(levelset_implicit_object_1,false,levelset_implicit_object_2,false);
            LOG::cout<<"combined eulerian object1="<<combined_eulerian.implicit_object1<<" object2="<<combined_eulerian.implicit_object2<<std::endl;
            combined_eulerian.Update_Box();
            combined_eulerian.alpha=0;
            combined_eulerian.dt=levelset_frequency;
            ((RIGID_BODY<TV>*)implicit_rigid_body)->Add_Structure(combined_eulerian);}
        else{
            IMPLICIT_OBJECT_COMBINED<TV> &combined=*new IMPLICIT_OBJECT_COMBINED<TV>(levelset_implicit_object_1,false,levelset_implicit_object_2,false);
            LOG::cout<<"object1="<<combined.implicit_object1<<" object2="<<combined.implicit_object2<<std::endl;
            combined.Update_Box();
            combined.alpha=0;
            ((RIGID_BODY<TV>*)implicit_rigid_body)->Add_Structure(combined);}}
    else{
        current_levelset++;
        head=&tests.Add_Rigid_Body(rigid_model,(T)1,(T)0,true);
        //head->is_static=true
        rigid_body_collection.rigid_body_particle.kinematic(head->particle_index)=true;
        init_frame=FRAME<TV>(head->X(),head->Rotation());
        implicit_rigid_body=head;}

    // ignore collisions on nodes that start inside
    for(int i=0;i<fixed_nodes.m;i++) deformable_body_collection.collisions.ignored_nodes.Append(fixed_nodes(i));
    if(restart)
        FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/ignored_nodes",deformable_body_collection.collisions.ignored_nodes);
    else{
        /*for(int i=0;i<edges.mesh.elements.m;i++) {
            const VECTOR<int,2> &nodes=edges.mesh.elements(i);
            if(implicit_rigid_body->Implicit_Geometry_Lazy_Inside(particles.X(nodes[0]))||implicit_rigid_body->Implicit_Geometry_Lazy_Inside(particles.X(nodes[1]))){
                deformable_body_collection.collisions.ignored_nodes.Append(nodes[0]);
                deformable_body_collection.collisions.ignored_nodes.Append(nodes[1]);}}
        Sort(deformable_body_collection.collisions.ignored_nodes);
        for(int i=deformable_body_collection.collisions.ignored_nodes.m;i>1;i--) if(deformable_body_collection.collisions.ignored_nodes(i)==deformable_body_collection.collisions.ignored_nodes(i-1)) deformable_body_collection.collisions.ignored_nodes.Remove_Index_Lazy(i);*/
        for(int i=0;i<particles.array_collection->Size();i++) 
            if(implicit_rigid_body->Implicit_Geometry_Lazy_Inside(particles.X(i))) deformable_body_collection.collisions.ignored_nodes.Append(i);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/ignored_nodes",deformable_body_collection.collisions.ignored_nodes);}

    comparator=new COLLISION_PAIR_COMPARATOR(&solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry,implicit_rigid_body);
    //comb
    if(test_number>6){
        ARRAY<int> tri_indicies;
        FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/"+sim_folder+"/interpolation",tri_indicies);
        for(int i=0;i<tri_indicies.m;i++){
            if(tri_indicies(i)<0) continue;
            interp_points.Append(implicit_rigid_body->World_Space_Point(implicit_rigid_body->simplicial_object->Get_Element(tri_indicies(i)).Center()));
            interp_normals.Append(implicit_rigid_body->simplicial_object->Get_Element(tri_indicies(i)).Normal());}
        RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("pick",(T).05,(T)0);
        init_pick=FRAME<TV>(sphere.X(),sphere.Rotation());
        sphere_id=sphere.particle_index;
        sphere.Is_Kinematic()=true;}
        
    //guide forces (head is needed)
    if (use_spring_guide&&!guide_sim_folder.empty()){
        ARRAY<int,HAIR_ID> roots(number_of_hairs);
        ARRAY<ARRAY<int>,HAIR_ID> hairs(next_segment_id);
        for(int p=0;p<particles.array_collection->Size();p++) hairs(particle_to_spring_id(p)).Append(p);
        for(HAIR_ID i(0);i<number_of_hairs;i++){
            T min_dist=-1;
            for (int p=1;p<=hairs(i).m;p++){
                T dist=implicit_rigid_body->implicit_object->Signed_Distance(particles.X(hairs(i)(p)));
                if (min_dist==-1||min_dist>dist) {min_dist=dist;roots(i)=hairs(i)(p);}}
            assert(min_dist!=-1);}
        guide_adhesion=new GUIDE_ADHESION<TV>(particles,edges.mesh,sim_guide_edges.mesh,particle_to_spring_id,roots);
        guide_adhesion->Set_Parameters(guide_stiffness,(T)10,guide_thickness,max_connections);
        guide_adhesion->Update_Springs(true);
        guide_adhesion->Write_State(stream_type,output_directory+"/adhesion.0");
        solid_body_collection.Add_Force(guide_adhesion);
    }
    //hairs
    deformable_body_collection.deformable_geometry.Add_Structure(&edges);
    if(use_guide){
        deformable_body_collection.deformable_geometry.Add_Structure(&sim_guide_edges);
        deformable_body_collection.deformable_geometry.Add_Structure(&sim_guide_tet_edges);}
    deformable_body_collection.deformable_geometry.Add_Structure(&fixed_edges);
    deformable_body_collection.deformable_geometry.Add_Structure(&extra_edges);
    deformable_body_collection.deformable_geometry.Add_Structure(&bending_edges);
    deformable_body_collection.deformable_geometry.Add_Structure(&torsion_edges);
    //Bindings
    ARRAY<TV> bindings1,bindings2; // Will be populated into ones to apply after MPI setup
    if(use_deforming_levelsets){
        FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/%s/fixed_positions.%d",data_directory.c_str(),sim_folder.c_str(),(current_levelset-1)),bindings1);
        FILE_UTILITIES::Read_From_File<T>(STRING_UTILITIES::string_sprintf("%s/%s/fixed_positions.%d",data_directory.c_str(),sim_folder.c_str(),current_levelset),bindings2);}
    else for(int i=0;i<fixed_nodes.m;i++){bindings1.Append(particles.X(fixed_nodes(i)));bindings2.Append(particles.X(fixed_nodes(i)));}
    // transform all points into world space
    for(int i=0;i<offset;i++) {particles.X(i)=implicit_rigid_body->World_Space_Point(particles.X(i));}

    // parallel
    if(solid_body_collection.deformable_body_collection.mpi_solids){
        solid_body_collection.deformable_body_collection.mpi_solids->KD_Tree_Partition_Subset(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,spring_id_to_particle,particles.X);
        ARRAY<PARTITION_ID>& partition_id_from_particle_index=solid_body_collection.deformable_body_collection.mpi_solids->partition_id_from_particle_index;
        ARRAY<ARRAY<int>,PARTITION_ID>& particles_of_partition=solid_body_collection.deformable_body_collection.mpi_solids->particles_of_partition;

        // give roots first good partition we run across if they don't have something already
        for(int p=0;p<particles.array_collection->Size();p++){
            int representative_particle=spring_id_to_particle(particle_to_spring_id(p));
            partition_id_from_particle_index(p)=partition_id_from_particle_index(representative_particle);}
        // update reverse map of all non root particles to new root (and connected component) asignment)
//        for(int p=0;p<particles.array_collection->Size();p++) partition_id_from_particle_index(p)=partition_id_from_particle_index(particle_connectivity.Find(p));
        // repopulate forward map
        for(PARTITION_ID i(0);i<particles_of_partition.Size();i++) particles_of_partition(i).Remove_All();
        for(int p=0;p<particles.array_collection->Size();p++) particles_of_partition(partition_id_from_particle_index(p)).Append(p);
        for(PARTITION_ID i(0);i<particles_of_partition.Size();i++) LOG::cout<<"Partition "<<i<<" has "<<particles_of_partition(i).Size()<<" particles"<<std::endl;
        //partition_fixed_nodes.Resize(solid_body_collection.deformable_body_collection.mpi_solids->Number_Of_Partitions());
        partition_spring_representative.Resize(solid_body_collection.deformable_body_collection.mpi_solids->Number_Of_Partitions());
        for(HAIR_ID hid(0);hid<spring_id_to_particle.Size();hid++) partition_spring_representative(partition_id_from_particle_index(spring_id_to_particle(hid))).Append(spring_id_to_particle(hid));

        // Restrict fixed nodes
        for(int i=fixed_nodes.m;i>=1;i--){
            //partition_fixed_nodes(partition_id_from_particle_index(fixed_nodes(i))).Append(fixed_nodes(i));
            if(partition_id_from_particle_index(fixed_nodes(i))==solid_body_collection.deformable_body_collection.mpi_solids->Partition()){
                fixed_nodes_indices.Append(i);
                fixed_nodes_to_apply.Append(fixed_nodes(i));
                bindings1_to_apply.Append(bindings1(i));
                bindings2_to_apply.Append(bindings2(i));}}
        
        // Restrict project_mesh to processors TODO: make this repartitionable
        for(int i=project_mesh.elements.m;i>=1;i--){
            const VECTOR<int,2>& nodes=project_mesh.elements(i);
            if(partition_id_from_particle_index(nodes(0))!=solid_body_collection.deformable_body_collection.mpi_solids->Partition() && partition_id_from_particle_index(nodes(0))!=solid_body_collection.deformable_body_collection.mpi_solids->Partition())
                project_mesh.elements.Remove_Index_Lazy(i);}
    }
    else{
        fixed_nodes_indices=IDENTITY_ARRAY<int>(fixed_nodes.m);
        fixed_nodes_to_apply=fixed_nodes;
        bindings1_to_apply=bindings1;
        bindings2_to_apply=bindings2;
    }

    // setup segment adhesion for partitions (be it single threaded  or not)
    // i.e. build local and boundary meshes and other mpi stuff for adhesion
    if(segment_adhesion) segment_adhesion->Update_Partitions(restart,solid_body_collection.deformable_body_collection.mpi_solids,output_directory);
        
    // initial projection rest lengths (needs to happen after MPI restriction)
    project_restlengths.Resize(project_mesh.elements.m);
    for(int i=0;i<project_mesh.elements.m;i++){
        const VECTOR<int,2>& nodes=project_mesh.elements(i);
        project_restlengths(i)=(deformable_body_collection.particles.X(nodes[0])-deformable_body_collection.particles.X(nodes[1])).Magnitude();}

    // compute binding velocities (needs to happen after MPI restriction)
    if(use_deforming_levelsets) Compute_Binding_Velocities();

    // collisions
    solid_body_collection.collision_body_list.Add_Bodies(*rigid_body_collection.rigid_geometry_collection.collision_body_list);
    solid_body_collection.deformable_body_collection.collisions.collision_structures.Append(&edges);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(&edges);
}
//#####################################################################
// Function Update_Levelsets_For_Time_Update
//#####################################################################
template<class T_input> template<class T_IMPLICIT_COMBINED> void HAIR_SIM_TESTS<T_input>::
Update_Keyframed_Parameters_For_Time_Update_Helper(const T time,T_IMPLICIT_COMBINED& combined)
{
    if(time<start_time) return;
    T pseudo_time=(time-start_time)/levelset_frequency;
    combined.alpha=(pseudo_time-(T)(int)pseudo_time);
    if (pseudo_time>=current_levelset){
        bindings1_to_apply=bindings2_to_apply;
        ARRAY<TV> bindings2;
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s/fixed_positions.%d",data_directory.c_str(),sim_folder.c_str(),++current_levelset),bindings2);
        bindings2_to_apply=bindings2.Subset(fixed_nodes_indices);
        Compute_Binding_Velocities();
        GRID<TV>& grid=*new GRID<TV>;
        ARRAY<T,VECTOR<int,3> >& phi=*new ARRAY<T,VECTOR<int,3> >;
        ARRAY<TV,VECTOR<int,3> >& velocity=*new ARRAY<TV,VECTOR<int,3> >;
        LEVELSET_3D<GRID<TV> > levelset(grid,phi);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s/motion/out.%d.phi",data_directory.c_str(),sim_folder.c_str(),current_levelset),levelset);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%s/motion/velocities.%d",data_directory.c_str(),sim_folder.c_str(),current_levelset),velocity);
        LOG::cout<<"PRE object1="<<combined.implicit_object1<<" object2="<<combined.implicit_object2<<std::endl;
        LOG::cout<<"READING LEVELSET NUMBER: "<<current_levelset<<std::endl;
        LEVELSET_IMPLICIT_OBJECT<TV> *levelset_implicit_object=new LEVELSET_IMPLICIT_OBJECT<TV>(levelset.grid,levelset.phi);
        levelset_implicit_object->V=&velocity;
        delete &((LEVELSET_IMPLICIT_OBJECT<TV>*)combined.implicit_object1)->levelset.grid;
        delete &((LEVELSET_IMPLICIT_OBJECT<TV>*)combined.implicit_object1)->levelset.phi;
        delete ((LEVELSET_IMPLICIT_OBJECT<TV>*)combined.implicit_object1)->V;
        delete combined.implicit_object1;
        combined.implicit_object1=combined.implicit_object2;
        combined.implicit_object2=levelset_implicit_object;
        combined.Update_Box();}
}
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Update_Keyframed_Parameters_For_Time_Update(const T time)
{
    if(time<start_time) return;
    if(use_eulerian_level_set_interpolation){
        IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>& implicit_combined_eulerian=
            ((IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>&)*((IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_BODY<TV> >*)implicit_rigid_body->implicit_object)->object_space_implicit_object);
        implicit_combined_eulerian.dt=levelset_frequency;
        Update_Keyframed_Parameters_For_Time_Update_Helper(time,implicit_combined_eulerian);}
    else
        Update_Keyframed_Parameters_For_Time_Update_Helper(time,
            ((IMPLICIT_OBJECT_COMBINED<TV>&)*((IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_BODY<TV> >*)implicit_rigid_body->implicit_object)->object_space_implicit_object));
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    for(int i=0;i<fixed_nodes_to_apply.m;i++) V(fixed_nodes_to_apply(i))=TV();
    if(guide_object1) for (int i=1;i<=guide_object1->particles.array_collection->Size();i++) {V(offset+i)=TV();}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(use_deforming_levelsets){
        Update_Keyframed_Parameters_For_Time_Update(velocity_time);
        for(int i=0;i<fixed_nodes_to_apply.m;i++) if(velocity_time<start_time) V(fixed_nodes_to_apply(i))=TV(); else V(fixed_nodes_to_apply(i))=binding_velocities(i);}
    else for(int i=0;i<fixed_nodes_to_apply.m;i++){
        TV point=bindings1_to_apply(i);
        V(fixed_nodes_to_apply(i))=implicit_rigid_body->Pointwise_Object_Velocity(implicit_rigid_body->World_Space_Point(point));}
    if(guide_object1){
        T frame=velocity_time*frame_rate;
        T alpha=frame-(T)(int)frame;
        if(frame>current_frame){
            guide_object1=guide_object2;
            guide_object2->Read(stream_type,guide_sim_folder+"/",++current_frame,-1,0,solids_parameters.write_from_every_process);}
        if(guide_object1) for (int i=1;i<=guide_object1->particles.array_collection->Size();i++) {V(offset+i)=(1-alpha)*guide_object1->particles.V(i)+alpha*guide_object2->particles.V(i);}}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    T pseudo_time=(time-start_time)/levelset_frequency;
    T alpha=(pseudo_time-(T)(int)pseudo_time);
    if(time<start_time) alpha=0;
    if(use_deforming_levelsets) Update_Keyframed_Parameters_For_Time_Update(time);
    for(int i=0;i<fixed_nodes_to_apply.m;i++){
        TV point;
        if(use_deforming_levelsets) point=(1-alpha)*bindings1_to_apply(i)+alpha*bindings2_to_apply(i);
        else point=bindings1_to_apply(i);
        X(fixed_nodes_to_apply(i))=implicit_rigid_body->World_Space_Point(point);}
    if(guide_object1){
        int frame=1;// TODO : broken fix before checkin
        T alpha=frame-(T)(int)frame;
        if(frame>current_frame){
            guide_object1=guide_object2;
            guide_object2->Read(stream_type,guide_sim_folder+"/",++current_frame,-1,0,solids_parameters.write_from_every_process);}
        if(guide_object1) for (int i=1;i<=guide_object1->particles.array_collection->Size();i++) {X(offset+i)=(1-alpha)*guide_object1->particles.X(i)+alpha*guide_object2->particles.X(i);}}
}
//#####################################################################
// Function Add_Connectivity
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Preprocess_Solids_Substep(const T time,const int substep) {
    if(segment_adhesion){segment_adhesion->Update_Collisions_List();segment_adhesion->Update_Springs(true);}
    //if(segment_adhesion){segment_adhesion->Update_Springs(true);}
    if(guide_adhesion) guide_adhesion->Update_Springs(false);
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Update_Time_Varying_Material_Properties(const T time)
{
    if(use_deforming_levelsets) Update_Keyframed_Parameters_For_Time_Update(time);
    //stop wind_drag
    //T wind_drag_off=(T)1.0;
    //T wind_drag_ramp_off=(T)1.1;
    T ether_drag_off=(T)1.5;
    T ether_drag_ramp_off=(T)1.6;
    PHYSBAM_ASSERT(!use_wind);
    //if(use_wind&&time>wind_drag_off){
    //    INTERPOLATION_CURVE<T,TV> wind_stop;
    //    wind_stop.Add_Control_Point(wind_drag_off,TV(0,(T).1,(T)-2.));
    //    wind_stop.Add_Control_Point(wind_drag_ramp_off,TV());
    //    WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
    //    drag.Use_Linear_Normal_Viscosity(0);drag.Use_Constant_Wind(0,TV());}
    if(time>ether_drag_off){
        INTERPOLATION_CURVE<T,T> wind_viscosity;
        wind_viscosity.Add_Control_Point(ether_drag_off,20);
        wind_viscosity.Add_Control_Point(ether_drag_ramp_off,0);
        ETHER_DRAG<GRID<TV> >& drag=solid_body_collection.template Find_Force<ETHER_DRAG<GRID<TV> >&>();
        drag.Use_Constant_Wind(wind_viscosity.Value(time));}
    if(use_drag&&time>start_time){
        ETHER_DRAG<GRID<TV> >& drag=solid_body_collection.template Find_Force<ETHER_DRAG<GRID<TV> >&>();
        drag.Use_Constant_Wind(drag_viscosity);}
    //start wind
    if(test_number%3!=1){
        if(time>wind_start_time){
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(17);drag.Use_Constant_Wind(0,TV(0,0,(T)-1.));}}
    //stop wind
    if(test_number%3==0){
        if(time>wind_stop_time){
            INTERPOLATION_CURVE<T,TV> wind_stop;
            wind_stop.Add_Control_Point(wind_stop_time-(T).3,TV(0,0,(T)-.25));
            wind_stop.Add_Control_Point(wind_stop_time,TV());
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(0);drag.Use_Constant_Wind(0,TV());}}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class T_input> bool HAIR_SIM_TESTS<T_input>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if(test_number>3 && test_number<=6 && id==head->particle_index){
        if(time>1&&time<10) twist.angular=TV((T)pi/2,0,0);}
    else if(test_number>6 && id==sphere_id){
        T duration = (T)2.;
        T start_time = (T)1.;
        INTERPOLATION_CURVE<T,TV> curve, normal;
        curve.Add_Control_Point(0,interp_points(0));
        normal.Add_Control_Point(0,interp_normals(0));
        for(int i=0;i<interp_points.m;i++) {
            curve.Add_Control_Point(start_time+duration*(T)i/interp_points.m,interp_points(i));
            normal.Add_Control_Point(start_time+duration*(T)i/interp_normals.m,interp_normals(i));}
        curve.Add_Control_Point(last_frame/frame_rate,interp_points(interp_points.m)+TV(0,0,(T)-100));
        normal.Add_Control_Point(last_frame/frame_rate,interp_normals(interp_normals.m));
        twist.linear=curve.Derivative(time);
        twist.angular=normal.Derivative(time);}
    else return false;
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if(test_number<=3 || test_number>6) frame=init_frame;
    if(test_number>3 && test_number<=6 && id==head->particle_index){
        if (time<=2||time>=5) frame=init_frame;
        else{
            FRAME<TV> transup, transdown, rotate;
            transup.t=TV(0,(T)0.2,0);
            transdown.t=TV(0,(T)-0.2,0);
            if (((int)(time*2))%2==0) rotate.r=ROTATION<TV>((time-(T)(int)time)*(T)pi,TV(1,0,0));
            else rotate.r=ROTATION<TV>((1-time+(T)(int)time)*(T)pi,TV(1,0,0));
            frame=init_frame*transdown*rotate*transup;
        }
    }
    if(test_number>6 && id==sphere_id){
        T duration = (T)2.;
        T start_time = (T)1.;
        FRAME<TV> rotation,transup,transhead;
        transup.t=TV(0,(T)(0.05*0.95),0);
        INTERPOLATION_CURVE<T,TV> curve, normal;
        curve.Add_Control_Point(0,interp_points(0)+TV(0,0,(T)10.));
        normal.Add_Control_Point(0,interp_normals(0));
        for(int i=0;i<interp_points.m;i++) {
            curve.Add_Control_Point(start_time+duration*(T)i/interp_points.m,interp_points(i));
            normal.Add_Control_Point(start_time+duration*(T)i/interp_normals.m,interp_normals(i));}
        curve.Add_Control_Point(last_frame/frame_rate,interp_points(interp_points.m)+TV(0,0,(T)-10.));
        normal.Add_Control_Point(last_frame/frame_rate,interp_normals(interp_normals.m));
        transhead.t=curve.Value(time);
        rotation.r=ROTATION<TV>::From_Rotated_Vector(TV(0,1,0),normal.Value(time));
        frame=transhead*rotation*transup;}
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Postprocess_Frame(const int frame)
{
    if (frame/frame_rate>=start_time) solids_parameters.triangle_collision_parameters.temporary_enable_collisions=solids_parameters.triangle_collision_parameters.perform_self_collision;
    //Write_Output_Files(frame);
    if(segment_adhesion) segment_adhesion->Write_State(stream_type,output_directory+STRING_UTILITIES::string_sprintf("/adhesion.%d",frame));
    if(guide_adhesion) guide_adhesion->Write_State(stream_type,output_directory+STRING_UTILITIES::string_sprintf("/adhesion.%d",frame));
}
//#####################################################################
// Function Add_External_Impulses_Before
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt)
{
    Add_External_Impulses_Helper(V,time,dt,use_momentum_conserving_before,use_non_momentum_conserving_before);
}
//#####################################################################
// Function Add_External_Impulses_Before
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt)
{
    Add_External_Impulses_Helper(V,time,dt,use_momentum_conserving_after,use_non_momentum_conserving_after);
}
//#####################################################################
// Function Add_External_Impulses_Helper
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Add_External_Impulses_Helper(ARRAY_VIEW<TV> V,const T time,const T dt,bool use_momentum_conserving,bool use_non_momentum_conserving)
{
    //if(time<start_time) return;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    if(write_substeps_level>-1){
        ARRAY<TV> positions_save(particles.X);
        for(int i=0;i<particles.array_collection->Size();i++) particles.X(i)=particles.X(i)+dt*V(i);
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
            if(project_second_node_to_surface && distance_to_root(child)==2){
                T phi;
                T tolerance=deformable_body_collection.collisions.collision_tolerances?(*deformable_body_collection.collisions.collision_tolerances)(child):deformable_body_collection.collisions.collision_tolerance;
                if(implicit_rigid_body->Implicit_Geometry_Lazy_Inside_And_Value(X_child_new_prime,phi,tolerance)){
                    X_child_new_prime-=phi*implicit_rigid_body->Implicit_Geometry_Normal(X_child_new_prime);}
            }
            V(child)=(X_child_new_prime-particles.X(child))/dt;}}

    if(write_substeps_level>-1){
        ARRAY<TV> positions_save(particles.X);
        for(int i=0;i<particles.array_collection->Size();i++) particles.X(i)=particles.X(i)+dt*V(i);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after impulse",2,2);
        particles.X=positions_save;}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_input> template<class T_IMPLICIT_COMBINED> void HAIR_SIM_TESTS<T_input>::
Write_Interpolated_Level_Set(const int frame,T_IMPLICIT_COMBINED& combined) const
{
    if(write_substeps_level){
        LEVELSET_IMPLICIT_OBJECT<TV> *levelset=(LEVELSET_IMPLICIT_OBJECT<TV>*)combined.implicit_object1;
        GRID<TV> grid(levelset->levelset.grid.counts.x,levelset->levelset.grid.counts.y,levelset->levelset.grid.counts.z,combined.box);
        ARRAY<T,VECTOR<int,3> > phi(grid.Domain_Indices());
        for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++) for(int ij=0;ij<grid.counts.z;ij++) phi(i,j,ij)=combined(grid.X(i,j,ij));
        LEVELSET_3D<GRID<TV> > interpolated_levelset(grid,phi);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/grid.%d",output_directory.c_str(),frame),grid);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/levelset.%d",output_directory.c_str(),frame),interpolated_levelset);}
    else{
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/grid.%d",output_directory.c_str(),frame),((LEVELSET_IMPLICIT_OBJECT<TV>*)combined.implicit_object1)->levelset.grid);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/levelset.%d",output_directory.c_str(),frame),((LEVELSET_IMPLICIT_OBJECT<TV>*)combined.implicit_object1)->levelset);}
}
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    if(segment_adhesion) segment_adhesion->Write_State(stream_type,output_directory+STRING_UTILITIES::string_sprintf("/adhesion.%d",frame));
    if(guide_adhesion) guide_adhesion->Write_State(stream_type,output_directory+STRING_UTILITIES::string_sprintf("/adhesion.%d",frame));
    if(use_deforming_levelsets){
        if(use_eulerian_level_set_interpolation) Write_Interpolated_Level_Set(frame,
            ((IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>&)*((IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_BODY<TV> >*)implicit_rigid_body->implicit_object)->object_space_implicit_object));
        else if(use_eulerian_level_set_interpolation) Write_Interpolated_Level_Set(frame,
            ((IMPLICIT_OBJECT_COMBINED<TV>&)*((IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_BODY<TV> >*)implicit_rigid_body->implicit_object)->object_space_implicit_object));
    }
}
//#####################################################################
// Function Point_Face_Mass
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,T& one_over_mass1,T& one_over_mass2,T& one_over_mass3,T& one_over_mass4){
    T factor = (T)0.5-sqr(attempt_ratio)/(T)2.;
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    const ARRAY<TV>& X=geometry.X_self_collision_free;
    VECTOR<T,3> face_embedded=weights[0]*X(nodes[1])+weights[1]*X(nodes[2])+weights[2]*X(nodes[3]);
    VECTOR<T,3> point=X(nodes[0]);
    if (implicit_rigid_body->implicit_object->Signed_Distance(face_embedded)<implicit_rigid_body->implicit_object->Signed_Distance(point))
    {one_over_mass2*=factor;one_over_mass3*=factor;one_over_mass4*=factor;}
    else one_over_mass1*=factor;
}
//#####################################################################
// Function Point_Face_Mass
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,VECTOR<T,4>& one_over_mass)
{Point_Face_Mass(attempt_ratio,nodes,weights,one_over_mass[0],one_over_mass[1],one_over_mass[2],one_over_mass[3]);}
//#####################################################################
// Function Point_Face_Mass
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Point_Face_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,3>& weights,ARRAY_VIEW<T>& one_over_mass){
    saved=one_over_mass.Subset(nodes);
    Point_Face_Mass(attempt_ratio,nodes,weights,one_over_mass(nodes[0]),one_over_mass(nodes[1]),one_over_mass(nodes[2]),one_over_mass(nodes[3]));
}
//#####################################################################
// Function Point_Face_Mass_Revert
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Point_Face_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass)
{Mass_Revert(nodes,one_over_mass);}
//#####################################################################
// Function Edge_Edge_Mass
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,T& one_over_mass1,T& one_over_mass2,T& one_over_mass3,T& one_over_mass4){
    //T factor = (T)0.5-sqr(attempt_ratio)/(T)2.;
    //TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry=solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry;
    //const ARRAY<TV>& X=geometry.X_self_collision_free;
    VECTOR<int,4> distances(distance_to_root.Subset(nodes));
    T distance1=(1-weights[0])*(T)distances[0]+weights[0]*(T)distances[1];
    T distance2=(1-weights[1])*(T)distances[2]+weights[1]*(T)distances[3];
    one_over_mass1*=distance1;one_over_mass2*=distance1;
    one_over_mass3*=distance2;one_over_mass4*=distance2;
}
//#####################################################################
// Function Edge_Edge_Mass
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,VECTOR<T,4>& one_over_mass)
{Edge_Edge_Mass(attempt_ratio,nodes,weights,one_over_mass[0],one_over_mass[1],one_over_mass[2],one_over_mass[3]);}
//#####################################################################
// Function Edge_Edge_Mass
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,4>& nodes,const VECTOR<T,2>& weights,ARRAY_VIEW<T>& one_over_mass){
    saved=one_over_mass.Subset(nodes);
    Edge_Edge_Mass(attempt_ratio,nodes,weights,one_over_mass(nodes[0]),one_over_mass(nodes[1]),one_over_mass(nodes[2]),one_over_mass(nodes[3]));
}
//#####################################################################
// Function Edge_Edge_Mass_Revert
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Edge_Edge_Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass)
{Mass_Revert(nodes,one_over_mass);}
//#####################################################################
// Function Mass_Revert
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Mass_Revert(const VECTOR<int,4>& nodes,ARRAY_VIEW<T>& one_over_mass)
{one_over_mass(nodes[0])=saved[0];one_over_mass(nodes[1])=saved[1];one_over_mass(nodes[2])=saved[2];one_over_mass(nodes[3])=saved[3];}
//#####################################################################
// Function Reorder_Pairs
//#####################################################################
template<class T_input> void HAIR_SIM_TESTS<T_input>::
Reorder_Pairs(ARRAY<VECTOR<int,4> >& edge_edge_pairs,ARRAY<VECTOR<int,4> >& point_face_pairs) {
    Sort(edge_edge_pairs,*comparator);
    Sort(point_face_pairs,*comparator);
}
//#####################################################################
template class HAIR_SIM_TESTS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HAIR_SIM_TESTS<double>;
#endif
