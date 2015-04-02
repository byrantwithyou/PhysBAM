//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Levine, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_BODY_COLLECTION
//#####################################################################
#include <Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Log/SCOPE.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/BINDING_SPRINGS.h>
#include <Deformables/Forces/DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Deformables/Forces/DEFORMABLES_FORCES.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLE_BODY_COLLECTION<TV>::
DEFORMABLE_BODY_COLLECTION(COLLISION_BODY_COLLECTION<TV>& collision_body_list)
    :particles(*new DEFORMABLE_PARTICLES<TV>),simulate(true),
    binding_list(*new BINDING_LIST<TV>(*this)),soft_bindings(*new SOFT_BINDINGS<TV>(binding_list)),mpi_solids(0),implicit_damping(true),
    print_diagnostics(false),print_residuals(false),print_energy(false),iterations_used_diagnostic(0),
    triangle_repulsions_and_collisions_geometry(*new TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>(*this)),
    triangle_repulsions(*new TRIANGLE_REPULSIONS<TV>(triangle_repulsions_and_collisions_geometry)),
    triangle_collisions(*new TRIANGLE_COLLISIONS<TV>(triangle_repulsions_and_collisions_geometry,triangle_repulsions.repulsion_thickness)),
    collisions(*new DEFORMABLE_OBJECT_COLLISIONS<TV>(particles,*this,structures,collision_body_list)),use_embedded_collisions(false),use_nonembedded_self_collision(false),
    check_stale(false)
{
    particles.Store_Velocity();
    particles.Store_Mass();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLE_BODY_COLLECTION<TV>::
~DEFORMABLE_BODY_COLLECTION()
{
    deformables_forces.Delete_Pointers_And_Clean_Memory();
    structures.Delete_Pointers_And_Clean_Memory();
    delete &soft_bindings;
    delete &binding_list;
    delete &collisions;
    delete &particles;
    delete &triangle_repulsions;
    delete &triangle_collisions;
    delete &triangle_repulsions_and_collisions_geometry;
}
//#####################################################################
// Function Initialize_Triangle_Collisions
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters)
{
    if(!triangle_collision_parameters.perform_self_collision || 
        (!triangle_repulsions_and_collisions_geometry.structures.m && !triangle_collision_parameters.initialize_collisions_without_objects))
        return;
    LOG::SCOPE scope("Initializing triangle collisions","Initializing triangle collisions");

    triangle_repulsions_and_collisions_geometry.Initialize(triangle_collision_parameters);
    if(triangle_collision_parameters.check_initial_mesh_for_self_intersection && triangle_repulsions_and_collisions_geometry.Check_For_Intersection(true)) PHYSBAM_FATAL_ERROR();
    triangle_repulsions.Initialize(triangle_collision_parameters);
    triangle_collisions.Initialize(triangle_collision_parameters);
}
//#####################################################################
// Function Update_Collision_Penalty_Forces_And_Derivatives
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Update_Collision_Penalty_Forces_And_Derivatives()
{
    for(int i=0;i<collision_penalty_forces.m;i++) collision_penalty_forces(i)->Update_Forces_And_Derivatives();
}
//#####################################################################
// Function Read_Dynamic_Variables
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,prefix+"/"+FILE_UTILITIES::Number_To_String(frame)+"/deformable_object_particles",particles);
    // if number==0, the particles format doesn't remember the set of attributes, so the following line makes restarts look more exact
    particles.Store_Velocity();
    std::string frame_string=FILE_UTILITIES::Number_To_String(frame);
    // read in bindings
    int local_frame=frame;
    std::string binding_state_list_name=prefix+"/common/bindings_list";
    std::string binding_state_name=prefix+"/"+frame_string+"/bindings";
    if(FILE_UTILITIES::File_Exists(binding_state_list_name)){
        if(!binding_list.frame_list){
            binding_list.frame_list=new ARRAY<int>;
            FILE_UTILITIES::Read_From_File(stream_type,binding_state_list_name,*binding_list.frame_list);}
        local_frame=(*binding_list.frame_list)(binding_list.frame_list->Binary_Search(frame));}
    if(binding_list.last_read!=local_frame && FILE_UTILITIES::File_Exists(binding_state_name)){
        FILE_UTILITIES::Read_From_File(stream_type,binding_state_name,binding_list);binding_list.last_read=local_frame;}
    local_frame=frame;
    std::string soft_binding_state_list_name=prefix+"/common/soft_bindings_list";
    std::string soft_binding_state_name=prefix+"/"+frame_string+"/soft_bindings";
    if(FILE_UTILITIES::File_Exists(soft_binding_state_list_name)){
        if(!soft_bindings.frame_list){soft_bindings.frame_list=new ARRAY<int>;FILE_UTILITIES::Read_From_File(stream_type,soft_binding_state_list_name,*soft_bindings.frame_list);}
        local_frame=(*soft_bindings.frame_list)(soft_bindings.frame_list->Binary_Search(frame));}
    if(soft_bindings.last_read!=local_frame && FILE_UTILITIES::File_Exists(soft_binding_state_name)){
        FILE_UTILITIES::Read_From_File(stream_type,soft_binding_state_name,soft_bindings);soft_bindings.last_read=local_frame;}
    // recompute auxiliary mass data (this data is destroyed when particles and read, and mass might have changed)
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Write_Dynamic_Variables
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const
{
    FILE_UTILITIES::Write_To_File(stream_type,prefix+"/"+FILE_UTILITIES::Number_To_String(frame)+"/deformable_object_particles",particles);
    std::string f=FILE_UTILITIES::Number_To_String(frame);
    if(binding_list.bindings.m>0 && !(check_stale && !binding_list.is_stale)){
        if(check_stale){
            if(!binding_list.frame_list) binding_list.frame_list=new ARRAY<int>;binding_list.frame_list->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"/common/bindings_list",*binding_list.frame_list);
            binding_list.is_stale=false;}
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"/"+f+"/bindings",binding_list);}
    if(soft_bindings.bindings.m>0 && !(check_stale && !soft_bindings.is_stale)){
        if(check_stale){
            if(!soft_bindings.frame_list) soft_bindings.frame_list=new ARRAY<int>;soft_bindings.frame_list->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"/common/soft_bindings_list",*soft_bindings.frame_list);
            soft_bindings.is_stale=false;}
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"/"+f+"/soft_bindings",soft_bindings);}
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,
    const bool read_from_every_process)
{
    particles.Store_Mass();
    if(include_static_variables) Read_Static_Variables(stream_type,static_prefix,static_frame);
    Read_Dynamic_Variables(stream_type,prefix,frame);
    if(mpi_solids && read_from_every_process) // make sure every process has all the correct data
        mpi_solids->Broadcast_Data(particles.X,particles.V);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,
    const bool write_from_every_process) const
{
    if(mpi_solids){
        DEFORMABLE_PARTICLES<TV>& const_particles=const_cast<DEFORMABLE_PARTICLES<TV>&>(particles);
        mpi_solids->Gather_Data(const_particles.X,const_particles.V);
        if(mpi_solids->rank && !write_from_every_process) return;}
    if(include_static_variables) Write_Static_Variables(stream_type,static_prefix,static_frame);
    Write_Dynamic_Variables(stream_type,prefix,frame);
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Update_Simulated_Particles()
{
    int particles_number=particles.Size();

    ARRAY<bool> particle_is_simulated(particles_number);for(int i=0;i<particles_number;i++) particle_is_simulated(i)=true;
    particle_is_simulated.Subset(particles.deletion_list).Fill(false);

    simulated_particles.Remove_All();
    dynamic_particles.Remove_All();

    if(mpi_solids){
        for(int p=0;p<mpi_solids->particles_of_partition(mpi_solids->Partition()).m;p++){
            int particle_index=mpi_solids->particles_of_partition(mpi_solids->Partition())(p);
            if(particle_is_simulated(particle_index)) simulated_particles.Append(particle_index);}
        
        binding_list.Clear_Hard_Bound_Particles(particle_is_simulated); // Prevent hard bound particles from being added to dynamic_particles
        
        for(int p=0;p<mpi_solids->particles_of_partition(mpi_solids->Partition()).m;p++){
            int particle_index=mpi_solids->particles_of_partition(mpi_solids->Partition())(p);
            if(particle_is_simulated(particle_index)) dynamic_particles.Append(particle_index);}

        SEGMENT_MESH dependency_mesh_force,dependency_mesh_binding;
        dependency_mesh_force.number_nodes=particles_number;
        dependency_mesh_force.Initialize_Incident_Elements();

        dependency_mesh_binding.number_nodes=particles_number;
        dependency_mesh_binding.Initialize_Incident_Elements();

        for(int i=0;i<deformables_forces.m;i++) deformables_forces(i)->Add_Dependencies(dependency_mesh_force);

        binding_list.Add_Dependencies(dependency_mesh_binding);
        soft_bindings.Add_Dependencies(dependency_mesh_binding);
        binding_list.Compute_Dependency_Closure_Based_On_Embedding(dependency_mesh_binding);

        dependency_mesh_force.Initialize_Neighbor_Nodes();
        dependency_mesh_binding.Initialize_Neighbor_Nodes();

        // setup mpi
        mpi_solids->Update(*this,dependency_mesh_force,dependency_mesh_binding);}
    else{
        for(int p=0;p<particles_number;p++) if(particle_is_simulated(p)) simulated_particles.Append(p);
        binding_list.Clear_Hard_Bound_Particles(particle_is_simulated); // Prevent hard bound particles from being added to dynamic_particles
        for(int p=0;p<particles_number;p++) if(particle_is_simulated(p)) dynamic_particles.Append(p);}

    ARRAY<bool> particle_is_simulated_actual(particles_number);
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> simulated_subset=particle_is_simulated_actual.Subset(simulated_particles);
    simulated_subset.Fill(true);
    for(int i=0;i<deformables_forces.m;i++) deformables_forces(i)->Update_Mpi(particle_is_simulated_actual,mpi_solids);

    if(mpi_solids){
        simulated_particles=dynamic_particles;
        for(PARTITION_ID n(0);n<mpi_solids->mpi_partition_binding.ghost_dynamic_particles.Size();n++)
            simulated_particles.Append_Elements(mpi_solids->mpi_partition_binding.ghost_dynamic_particles(n));
        for(PARTITION_ID n(0);n<mpi_solids->mpi_partition_force.ghost_dynamic_particles.Size();n++)
            simulated_particles.Append_Elements(mpi_solids->mpi_partition_force.ghost_dynamic_particles(n));
        simulated_particles.Prune_Duplicates();
        binding_list.Compute_Particle_Closure_Based_On_Embedding(simulated_particles);}

    collisions.Update_Simulated_Particles();
}
//#####################################################################
// Function Set_Mpi_Solids
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Set_Mpi_Solids(MPI_SOLIDS<TV>* mpi_solids_input)
{
    mpi_solids=mpi_solids_input;
    triangle_repulsions.mpi_solids=mpi_solids;
    triangle_collisions.mpi_solids=mpi_solids;
    triangle_repulsions_and_collisions_geometry.mpi_solids=mpi_solids;
}
//#####################################################################
// Function Update_CFL
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Update_CFL()
{
    bool cfl_valid=true;
    if(deformables_forces.m){
        for(int i=0;i<deformables_forces.m;i++){if(!deformables_forces(i)->CFL_Valid()){cfl_valid=false;break;}}}
    else cfl_valid=false;
    if(!cfl_valid){
        frequency.Resize(particles.Size(),false,false);
        INDIRECT_ARRAY<ARRAY<T_FREQUENCY_DEFORMABLE>,ARRAY<int>&> frequency_subset=frequency.Subset(simulated_particles);
        frequency_subset.Fill(T_FREQUENCY_DEFORMABLE());

        for(int i=0;i<deformables_forces.m;i++){deformables_forces(i)->Initialize_CFL(frequency);deformables_forces(i)->Validate_CFL();}
        cfl_elastic=FLT_MAX;cfl_damping=FLT_MAX;
        for(int i=0;i<simulated_particles.m;i++){int p=simulated_particles(i);
            cfl_elastic=min(cfl_elastic,Robust_Inverse(sqrt(frequency(p).elastic_squared)));
            cfl_damping=min(cfl_damping,Robust_Inverse(frequency(p).damping));}}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_BODY_COLLECTION<TV>::
CFL(const bool verbose)
{
    T dt_elastic_and_damping=CFL_Elastic_And_Damping(),dt_strain_rate=CFL_Strain_Rate();
    if(verbose){
        LOG::cout<<"dt_elastic_and_damping = "<<dt_elastic_and_damping<<std::endl;
        LOG::cout<<"dt_strain_rate = "<<dt_strain_rate<<std::endl;
        LOG::cout<<"min = "<<min(dt_elastic_and_damping,dt_strain_rate)<<std::endl;}
    return min(dt_elastic_and_damping,dt_strain_rate);
}
//#####################################################################
// Function CFL_Elastic_And_Damping
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_BODY_COLLECTION<TV>::
CFL_Elastic_And_Damping()
{
    T dt_elastic=CFL_Elastic();
    T dt_damping=FLT_MAX;if(!implicit_damping) dt_damping=CFL_Damping();
    T one_over_dt_full=1/dt_elastic+1/dt_damping;
    return Robust_Divide((T)1,one_over_dt_full);
}
//#####################################################################
// Function CFL_Elastic
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_BODY_COLLECTION<TV>::
CFL_Elastic()
{
    Update_CFL();
    return cfl_elastic;
}
//#####################################################################
// Function CFL_Damping
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_BODY_COLLECTION<TV>::
CFL_Damping()
{
    Update_CFL();
    return cfl_damping;
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_BODY_COLLECTION<TV>::
CFL_Strain_Rate()
{
    T dt_strain=FLT_MAX;
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->limit_time_step_by_strain_rate) dt_strain=min(dt_strain,deformables_forces(k)->CFL_Strain_Rate()); // otherwise not included
    return dt_strain;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    for(int k=0;k<deformables_forces.m;k++) deformables_forces(k)->Update_Position_Based_State(time,is_position_update);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F_full,const T time) const
{
    assert(F_full.Size()==particles.Size());
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->use_velocity_independent_forces) deformables_forces(k)->Add_Velocity_Independent_Forces(F_full,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
// can depend on position too
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T time) const
{
    assert(F_full.Size()==particles.Size());
    for(int k=0;k<deformables_forces.m;k++) if(deformables_forces(k)->use_velocity_dependent_forces) deformables_forces(k)->Add_Velocity_Dependent_Forces(V_full,F_full,time);
}
//#####################################################################
// Function Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T scale,const T time) const
{
    assert(V_full.Size()==particles.Size() && F_full.Size()==particles.Size());
    for(int k=0;k<deformables_forces.m;k++)
        if(deformables_forces(k)->use_implicit_velocity_independent_forces)
            deformables_forces(k)->Add_Implicit_Velocity_Independent_Forces(V_full,F_full,scale,time);
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int DEFORMABLE_BODY_COLLECTION<TV>::
Add_Force(DEFORMABLES_FORCES<TV>* force)
{
    deformables_forces.Append(force);
    force->Set_CFL_Number(cfl_number);
    return deformables_forces.m;
}
//#####################################################################
// Function Set_CFL_Number
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Set_CFL_Number(const T cfl_number_input)
{
    cfl_number=cfl_number_input;
    for(int i=0;i<deformables_forces.m;i++) deformables_forces(i)->Set_CFL_Number(cfl_number_input);
}
//#####################################################################
// Function Test_Energy
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Test_Energy(const T time)
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    RANDOM_NUMBERS<T> random;
    T e=(T)1e-5;
    ARRAY<TV> dX(particles.X.m);
    random.Fill_Uniform(dX,-e,e);
    ARRAY<TV> X2a(particles.X+dX);
    ARRAY_VIEW<TV> X1(X2a);
    for(int i=0;i<deformables_forces.m;i++){
        ARRAY<TV> F(particles.X.m);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        T PE0=deformables_forces(i)->Potential_Energy(time);
        particles.X.Exchange(X1);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        T PE1=deformables_forces(i)->Potential_Energy(time);
        particles.X.Exchange(X1);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        T W=F.Dot(dX)/2;
        T dPE=(PE0-PE1)/e,dW=W/e,rel=(dPE-dW)/max(abs(dW),(T)1e-20);
        LOG::cout<<"potential energy test d phi "<<dPE<<"  W "<<dW<<"   rel "<<rel<<"   "<<typeid(*deformables_forces(i)).name()<<std::endl;}
}
//#####################################################################
// Function Test_Force_Derivatives
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Test_Force_Derivatives(const T time)
{
    PHYSBAM_ASSERT(sizeof(T)==sizeof(double));
    RANDOM_NUMBERS<T> random;
    T e=(T)1e-5;
    ARRAY<TV> dX(particles.X.m);
    random.Fill_Uniform(dX,-e,e);
    ARRAY<TV> X2a(particles.X+dX);
    ARRAY_VIEW<TV> X1(X2a);
    for(int i=0;i<deformables_forces.m;i++){
        ARRAY<TV> F(particles.X.m),G(particles.X.m);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        deformables_forces(i)->Add_Implicit_Velocity_Independent_Forces(dX,G,(T).5,time);
        F*=-(T)1;
        particles.X.Exchange(X1);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        deformables_forces(i)->Add_Velocity_Independent_Forces(F,time);
        deformables_forces(i)->Add_Implicit_Velocity_Independent_Forces(dX,G,(T).5,time);
        particles.X.Exchange(X1);
        deformables_forces(i)->Update_Position_Based_State(time,true);
        T MF=sqrt(F.Magnitude_Squared());
        T MG=sqrt(G.Magnitude_Squared());
        T MD=sqrt((F-G).Magnitude_Squared());
        LOG::cout<<"force derivative error "<<MD<<" vs "<<MF<<"   "<<MG<<"   rel "<<MD/max((T)1e-30,MF,MG)<<"    "<<typeid(*deformables_forces(i)).name()<<std::endl;}
}
//#####################################################################
// Function Read_Static_Variables
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Read_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame)
{
    std::string f=frame==-1?"common":FILE_UTILITIES::Number_To_String(frame);
    std::istream* input_raw=FILE_UTILITIES::Safe_Open_Input(prefix+"/"+f+"/deformable_object_structures");
    TYPED_ISTREAM input(*input_raw,stream_type);
    int m;Read_Binary(input,m);
    // TODO: merge this functionality with dynamic lists to allow for more flexibility
    if(!structures.m){ // // create and read all structures from scratch
        structures.Resize(m);
        for(int k=0;k<structures.m;k++) structures(k)=STRUCTURE<TV>::Create_Structure(input,particles);}
    else if(structures.m<=m){
        int old_number_of_structures=structures.m;structures.Resize(m);
        for(int k=0;k<old_number_of_structures;k++) structures(k)->Read_Structure(input);
        for(int k=old_number_of_structures;k<m;k++) structures(k)=STRUCTURE<TV>::Create_Structure(input,particles);}
    else{
        LOG::cout<<"Current number of structures ("<<structures.m<<") is greater than number in file ("<<m<<").";
        PHYSBAM_FATAL_ERROR();}
    delete input_raw;
}
//#####################################################################
// Function Write_Static_Variables
//#####################################################################
template<class TV> void DEFORMABLE_BODY_COLLECTION<TV>::
Write_Static_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const
{
    std::string f=frame==-1?"common":FILE_UTILITIES::Number_To_String(frame);
    std::ostream* output_raw=FILE_UTILITIES::Safe_Open_Output(prefix+"/"+f+"/deformable_object_structures");
    TYPED_OSTREAM output(*output_raw,stream_type);
    Write_Binary(output,structures.m);
    for(int k=0;k<structures.m;k++) structures(k)->Write_Structure(output);
    delete output_raw;
}
//#####################################################################
namespace PhysBAM{
template class DEFORMABLE_BODY_COLLECTION<VECTOR<float,1> >;
template class DEFORMABLE_BODY_COLLECTION<VECTOR<float,2> >;
template class DEFORMABLE_BODY_COLLECTION<VECTOR<float,3> >;
template class DEFORMABLE_BODY_COLLECTION<VECTOR<double,1> >;
template class DEFORMABLE_BODY_COLLECTION<VECTOR<double,2> >;
template class DEFORMABLE_BODY_COLLECTION<VECTOR<double,3> >;
}
