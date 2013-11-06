//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Log/LOG.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Utilities/DEBUG_CAST.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <climits>
namespace PhysBAM{
template<class TV>
struct ALLOCATE_BODY_HELPER:public ALLOCATE_HELPER<TV>
{
    RIGID_BODY_COLLECTION<TV>& collection;
    ALLOCATE_BODY_HELPER(RIGID_BODY_COLLECTION<TV>& collection_input):collection(collection_input) {}
    RIGID_BODY<TV>* Create(int index=0) PHYSBAM_OVERRIDE {return new RIGID_BODY<TV>(collection,true,index);}
    virtual ~ALLOCATE_BODY_HELPER(){}
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_COLLECTION<TV>::
RIGID_BODY_COLLECTION(COLLISION_BODY_COLLECTION<TV>* collision_body_list_input)
    :rigid_body_particles(*new RIGID_BODY_PARTICLES<TV>()),collision_body_list(collision_body_list_input),structure_list(*new STRUCTURE_LIST<TV,int>),always_create_structure(false),
    structure_hash(*new HASHTABLE<std::string,int>),frame_list_key(0),frame_list_active(0),check_stale(false),last_read_key(-1),last_read_active(-1),
    allocate_helper(new ALLOCATE_BODY_HELPER<TV>(*this)),owns_collision_body_list(false),
    articulated_rigid_body(*new ARTICULATED_RIGID_BODY<TV>(*this)),
    rigid_body_cluster_bindings(*new RIGID_BODY_CLUSTER_BINDINGS<TV>(*this,articulated_rigid_body)),
    dynamic_rigid_body_particles(0),print_diagnostics(false),print_residuals(false),print_energy(false),iterations_used_diagnostic(0)
{
    if(!allocate_helper) allocate_helper=new ALLOCATE_BODY_HELPER<TV>(*this);
    if(!collision_body_list){
        collision_body_list=new COLLISION_BODY_COLLECTION<TV>;
        owns_collision_body_list=true;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_COLLECTION<TV>::
~RIGID_BODY_COLLECTION()
{
    rigids_forces.Delete_Pointers_And_Clean_Memory();
    delete &articulated_rigid_body;
    delete &rigid_body_cluster_bindings;
    delete &rigid_body_particles;
    delete &structure_hash;delete allocate_helper;delete &structure_list;
    if(owns_collision_body_list) delete collision_body_list;
}
//#####################################################################
// Function Rigid_Body
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGID_BODY_COLLECTION<TV>::
Rigid_Body(const int particle_index)
{
    return *rigid_body_particles.rigid_body(particle_index);
}
//#####################################################################
// Function Rigid_Body
//#####################################################################
template<class TV> const RIGID_BODY<TV>& RIGID_BODY_COLLECTION<TV>::
Rigid_Body(const int particle_index) const
{
    return *rigid_body_particles.rigid_body(particle_index);
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
// structures already added to their respective lists
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body(RIGID_BODY<TV>* rigid_body,const int simplicial_boundary_id,const int implicit_object_id,const int simplicial_interior_id)
{
    int id=rigid_body->particle_index;
    if(simplicial_boundary_id>=0) rigid_body_particles.structure_ids(id)(0)=simplicial_boundary_id;
    if(implicit_object_id>=0) rigid_body_particles.structure_ids(id)(1)=implicit_object_id;
    if(simplicial_interior_id>=0) rigid_body_particles.structure_ids(id)(2)=simplicial_interior_id;
    for(int i=0;i<rigid_body_particles.structure_ids(id).m;i++)
        if(rigid_body_particles.structure_ids(id)(i)>=0 && !structure_list.Element(rigid_body_particles.structure_ids(id)(i)))
            PHYSBAM_FATAL_ERROR();
    return id;
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
// adds body's segmented curve and implicit curve to their respective lists
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body_And_Geometry(RIGID_BODY<TV>* rigid_body)
{
    int id=rigid_body->particle_index;
    assert(rigid_body->structures.m<=3);
    rigid_body_particles.structure_ids(id)=VECTOR<int,3>(-1,-1,-1);
    for(int s=0;s<rigid_body->structures.m;s++) rigid_body_particles.structure_ids(id)(s)=structure_list.Add_Element(rigid_body->structures(s));
    return id;
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body(const STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,const bool read_simplicial_boundary,const bool read_implicit_object,
    const bool read_simplicial_interior,const bool read_rgd_file)
{
    return Add_Rigid_Body(stream_type,false,basename,scaling_factor,read_simplicial_boundary,read_implicit_object,read_simplicial_interior,read_rgd_file);
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body(const STREAM_TYPE stream_type,const bool thin_shell,const std::string& basename,const T scaling_factor,const bool read_simplicial_boundary,const bool read_implicit_object,
    const bool read_simplicial_interior,const bool read_rgd_file)
{
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(*this,true);
    rigid_body->thin_shell=thin_shell;

    // rigid body
    std::string rgd=TV::dimension==2?"rgd2d":"rgd";
    if(read_rgd_file){
        try{FILE_UTILITIES::Read_From_File(stream_type,basename+"."+rgd,rigid_body->Mass(),rigid_body->Inertia_Tensor(),rigid_body->Frame());}
        catch(FILESYSTEM_ERROR&){LOG::cout<<"Note: No "<<rgd<<" file for "<<basename<<" (using default values)"<<std::endl;}}
    if(scaling_factor!=1) rigid_body->Rescale(scaling_factor);
    rigid_body->Update_Angular_Velocity();

    int id=Add_Rigid_Body(rigid_body,stream_type,basename,scaling_factor,read_simplicial_boundary,read_implicit_object,read_simplicial_interior,read_rgd_file);
    return id;
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body(RIGID_BODY<TV>* rigid_body,STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,
    const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file)
{
    int id=rigid_body->particle_index;

    // structures
    ARRAY<int> structure_ids;
    TV structure_center=rigid_body_particles.frame(id).t;
    if(TV::dimension==2){
        if(read_simplicial_boundary && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".curve2d",scaling_factor,structure_center))
            LOG::cout<<"Note: No curve2d file for "<<basename<<std::endl;
        if(read_implicit_object && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".phi2d",scaling_factor,structure_center))
            LOG::cout<<"Note: No phi2d file for "<<basename<<std::endl;
        if(read_simplicial_interior && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".tri2d",scaling_factor,structure_center))
            LOG::cout<<"Note: No tri2d file for "<<basename<<std::endl;}
    else{
        if(read_simplicial_boundary && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".tri",scaling_factor,structure_center))
            LOG::cout<<"Note: No tri file for "<<basename<<std::endl;
        if(read_implicit_object && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".phi",scaling_factor,structure_center) && 
            !Find_Or_Read_Structure(stream_type,structure_ids,basename+".oct",scaling_factor,structure_center))
            LOG::cout<<"Note: No phi or oct file for "<<basename<<std::endl;
        if(read_simplicial_interior && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".tet",scaling_factor,structure_center))
            LOG::cout<<"Note: No tet file for "<<basename<<std::endl;}
    assert(structure_ids.m<=3);
    rigid_body_particles.structure_ids(id)=VECTOR<int,3>(-1,-1,-1);
    for(int i=0;i<structure_ids.m;i++){
        if(structure_ids(i)>=0){
            rigid_body_particles.structure_ids(id)(i)=structure_ids(i);
            Rigid_Body(id).Add_Structure(*structure_list.Element(structure_ids(i)));}}

    return id;
}
//#####################################################################
// Function Update_Angular_Velocity
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Velocity()
{
    for(int p=0;p<rigid_body_particles.Size();p++) if(Is_Active(p)) Rigid_Body(p).Update_Angular_Velocity();
}
//#####################################################################
// Function Update_Angular_Momentum
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Momentum()
{
    for(int p=0;p<rigid_body_particles.Size();p++) if(Is_Active(p)) Rigid_Body(p).Update_Angular_Momentum();
}
//#####################################################################
// Function Update_Angular_Velocity
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Velocity(const ARRAY<int>& particle_indices)
{
    for(int i=0;i<particle_indices.m;i++) Rigid_Body(particle_indices(i)).Update_Angular_Velocity();
}
//#####################################################################
// Function Update_Angular_Momentum
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Momentum(const ARRAY<int>& particle_indices)
{
    for(int i=0;i<particle_indices.m;i++) Rigid_Body(particle_indices(i)).Update_Angular_Momentum();
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Simulated_Particles()
{
    int rigid_particles_number=rigid_body_particles.Size();

    ARRAY<bool> particle_is_simulated(rigid_particles_number);
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> simulated_subset=particle_is_simulated.Subset(rigid_body_particles.deletion_list);
    simulated_subset.Fill(false);
    for(int i=0;i<rigid_particles_number;i++)
        if(Is_Active(i) && Rigid_Body(i).Is_Simulated()) // TODO: Can't everything be defaulted to true?
            particle_is_simulated(i)=true;
    for(int i=0;i<rigid_particles_number;i++)
        if(!Is_Active(i) || !Rigid_Body(i).Is_Simulated())
            particle_is_simulated(i)=false;

    simulated_rigid_body_particles.Remove_All();
    dynamic_rigid_body_particles.Remove_All();

    for(int p=0;p<rigid_particles_number;p++) if(particle_is_simulated(p)) simulated_rigid_body_particles.Append(p);

    rigid_body_cluster_bindings.Clear_Hard_Bound_Particles(particle_is_simulated);

    for(int p=0;p<rigid_particles_number;p++) if(particle_is_simulated(p)) dynamic_rigid_body_particles.Append(p);

    static_rigid_bodies.Remove_All();kinematic_rigid_bodies.Remove_All();static_and_kinematic_rigid_bodies.Remove_All();
    for(int p=0;p<rigid_particles_number;p++) if(Is_Active(p)){RIGID_BODY<TV>& rigid_body=Rigid_Body(p);
        if(rigid_body.is_static){static_rigid_bodies.Append(p);static_and_kinematic_rigid_bodies.Append(p);}
        if(rigid_body_particles.kinematic(p)){kinematic_rigid_bodies.Append(p);static_and_kinematic_rigid_bodies.Append(p);}}

    ARRAY<bool> rigid_particle_is_simulated(rigid_particles_number);
    rigid_particle_is_simulated.Subset(simulated_rigid_body_particles).Fill(true);
    for(int i=0;i<rigids_forces.m;i++) rigids_forces(i)->Update_Mpi(rigid_particle_is_simulated);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=0;k<rigids_forces.m;k++)
        if(rigids_forces(k)->use_velocity_independent_forces) rigids_forces(k)->Add_Velocity_Independent_Forces(rigid_F_full,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
// can depend on position too
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=0;k<rigids_forces.m;k++)
        if(rigids_forces(k)->use_velocity_dependent_forces) rigids_forces(k)->Add_Velocity_Dependent_Forces(rigid_V_full,rigid_F_full,time);
}
//#####################################################################
// Function Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T scale,const T time) const
{
    assert(rigid_F_full.Size()==rigid_body_particles.Size());
    for(int k=0;k<rigids_forces.m;k++)
        if(rigids_forces(k)->use_implicit_velocity_independent_forces)
            rigids_forces(k)->Add_Implicit_Velocity_Independent_Forces(rigid_V_full,rigid_F_full,scale,time);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Position_Based_State(const T time)
{
    for(int k=0;k<rigids_forces.m;k++) rigids_forces(k)->Update_Position_Based_State(time);
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const
{
    potential_energy=0;
    kinetic_energy=0;
    for(int i=0;i<rigids_forces.m;i++) potential_energy+=rigids_forces(i)->Potential_Energy(time);
    for(int i=0;i<dynamic_rigid_body_particles.m;i++){int p=dynamic_rigid_body_particles(i);
        kinetic_energy+=Rigid_Body(p).Kinetic_Energy();}
}
//#####################################################################
// Function Print_Energy
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Print_Energy(const T time,const int step) const
{
    if(print_energy){
        T potential_energy=0,kinetic_energy=0;
        Compute_Energy(time,kinetic_energy,potential_energy);
        LOG::cout<<"total energy = "<<(potential_energy+kinetic_energy)<<"    (KE = "<<kinetic_energy<<"   PE = "<<potential_energy<<")  Step "<<step<<std::endl;}
}
//#####################################################################
// Function CFL_Rigid
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY_COLLECTION<TV>::
CFL_Rigid(const RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters,const bool verbose_dt)
{
    static T static_min_bounding_box_width=FLT_MAX;
    T min_bounding_box_width=FLT_MAX;
    for(int i=0;i<rigid_body_particles.Size();i++) if(Is_Active(i)){
            const RANGE<TV>& box=Rigid_Body(i).Object_Space_Bounding_Box();
            TV edge_lengths=box.Edge_Lengths();min_bounding_box_width=min(min_bounding_box_width,edge_lengths.Min());}
    if(min_bounding_box_width!=static_min_bounding_box_width){
        static_min_bounding_box_width=min_bounding_box_width;
        LOG::Stat("minimum rigid body bounding box width",min_bounding_box_width);}

    T max_distance_per_time_step=rigid_body_evolution_parameters.max_rigid_body_linear_movement_fraction_per_time_step*min_bounding_box_width;
    T dt=FLT_MAX;
    bool no_active_bodies=true;
    for(int p=0;p<rigid_body_particles.Size();p++) if(Is_Active(p)){
        dt=min(dt,Rigid_Body(p).CFL(max_distance_per_time_step,rigid_body_evolution_parameters.max_rigid_body_rotation_per_time_step,verbose_dt));
        no_active_bodies=false;}
    if(no_active_bodies) return FLT_MAX; // don't apply rigid dt bounds if there aren't any active rigid bodies
    dt=Robust_Multiply(rigid_body_evolution_parameters.rigid_cfl,dt);
    T dt_clamped=clamp(dt,rigid_body_evolution_parameters.rigid_minimum_dt,rigid_body_evolution_parameters.rigid_maximum_dt);
    if(dt_clamped>dt && verbose_dt) LOG::cout<<"Warning: taking larger time step ("<<dt_clamped<<") than CFL dt ("<<dt<<")"<<std::endl;
    return dt_clamped;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Force(RIGIDS_FORCES<TV>* force)
{
    rigids_forces.Append(force);
    force->Set_CFL_Number((T).5);
    return rigids_forces.m;
}
//#####################################################################
// Function Update_Kinematic_Particles
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Kinematic_Particles()
{
    static_rigid_bodies.Remove_All();
    kinematic_rigid_bodies.Remove_All();
    for(int p=0;p<rigid_body_particles.Size();p++)
        if(Is_Active(p)){
            RIGID_BODY<TV>& rigid_body=Rigid_Body(p);
            if(rigid_body.is_static) static_rigid_bodies.Append(p);
            else kinematic_rigid_bodies.Append(p);}
}
//#####################################################################
// Function Register_Analytic_Replacement_Structure
//#####################################################################
template<class TV> bool RIGID_BODY_COLLECTION<TV>::
Register_Analytic_Replacement_Structure(const std::string& filename,const T scaling_factor,STRUCTURE<TV>* structure)
{
    std::string hashname=STRING_UTILITIES::string_sprintf("%s@%.6f",filename.c_str(),scaling_factor); // mangle hash name
    if(structure_hash.Contains(hashname)) return false;
    int id=structure?structure_list.Add_Element(structure):0;
    structure_hash.Insert(hashname,id);
    return true;
}
template<class T> void
Wrap_Structure_Helper(STRUCTURE<VECTOR<T,1> >*& structure,const VECTOR<T,1>& center)
{}
template<class TV> void
Wrap_Structure_Helper(STRUCTURE<TV>*& structure,const TV& center)
{   // TODO(jontg): This shouldn't be necessary
}
//#####################################################################
// Function Find_Or_Read_Structure
//#####################################################################
template<class TV> bool RIGID_BODY_COLLECTION<TV>::
Find_Or_Read_Structure(const STREAM_TYPE stream_type,ARRAY<int>& structure_ids,const std::string& filename,const T scaling_factor,const TV& center)
{
    int id;
    if(!FILE_UTILITIES::File_Exists(filename)) return false;
    std::string hashname=STRING_UTILITIES::string_sprintf("%s@%.6f",filename.c_str(),scaling_factor); // mangle hash name
    if(!always_create_structure&&structure_hash.Get(hashname,id)){ // already read in
        if(!structure_list.Is_Active(id)) PHYSBAM_FATAL_ERROR();} // // only works if the referenced geometry is still in memory
    else{ // read in for the first time
        STRUCTURE<TV>* structure=0;
        if(!stream_type.use_doubles)
            structure=STRUCTURE<TV>::template Create_From_File<float>(filename);
        else
            structure=STRUCTURE<TV>::template Create_From_File<double>(filename);
        if(scaling_factor!=1){
            Wrap_Structure_Helper(structure,center);
            structure->Rescale(scaling_factor);}
        id=structure_list.Add_Element(structure);
        if(!always_create_structure) structure_hash.Insert(hashname,id);}
    structure_ids.Append(id);
    return true;
}
//#####################################################################
// Function Destroy_Unreferenced_Geometry
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Destroy_Unreferenced_Geometry() 
{
    ARRAY<bool> referenced(structure_list.Number_Of_Active_Elements());
    for(int i=0;i<rigid_body_particles.Size();i++) for(int j=0;j<rigid_body_particles.structure_ids(i).m;j++) if(rigid_body_particles.structure_ids(i)(j)>=0)
        referenced(structure_list.Element_Index(rigid_body_particles.structure_ids(i)(j)))=true;
    for(int i=structure_list.Number_Of_Active_Elements()-1;i>=0;i--) if(!referenced(i)) structure_list.Remove_Element(structure_list.Active_Element_Id(i));
}
//#####################################################################
// Function Destroy_Unreferenced_Geometry
//#####################################################################
template<class TV> RIGID_BODY<TV>* RIGID_BODY_COLLECTION<TV>::
New_Body(int index)
{
    return allocate_helper->Create(index);
}
//#####################################################################
// Function Update_Level_Set_Transforms
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Level_Set_Transforms()
{
    for(int i=0;i<rigid_body_particles.Size();i++)
        if(rigid_body_particles.rigid_body(i))
            if(rigid_body_particles.rigid_body(i)->implicit_object)
                rigid_body_particles.rigid_body(i)->implicit_object->transform=&rigid_body_particles.frame(i);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame,ARRAY<int>* needs_init,ARRAY<int>* needs_destroy)
{
    if(stream_type.use_doubles) structure_list.template Read<double>(directory,"rigid_body_structure_",frame);
    else structure_list.template Read<float>(directory,"rigid_body_structure_",frame);
    ARRAY<RIGID_BODY<TV>*> bodies(rigid_body_particles.rigid_body);
    rigid_body_particles.rigid_body.Fill(0);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_particles",directory.c_str(),frame),rigid_body_particles);
    while(rigid_body_particles.rigid_body.m<bodies.m) delete bodies.Pop();
    rigid_body_particles.rigid_body.Prefix(bodies.m)=bodies;

    ARRAY<int> needs_init_default;
    if(needs_init) needs_init->Remove_All();
    if(needs_destroy) needs_destroy->Remove_All();
    char version;int next_id=0;ARRAY<int> active_ids;int local_frame=frame;
    std::string active_list_name=STRING_UTILITIES::string_sprintf("%s/common/rigid_body_active_ids_list",directory.c_str(),frame);
    std::string active_name=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_active_ids",directory.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(active_list_name)){
        if(!frame_list_active){frame_list_active=new ARRAY<int>;FILE_UTILITIES::Read_From_File(stream_type,active_list_name,*frame_list_active);}
        local_frame=(*frame_list_active)(frame_list_active->Binary_Search(frame));}
    if(last_read_active!=local_frame && FILE_UTILITIES::File_Exists(active_name)){
        FILE_UTILITIES::Read_From_File(stream_type,active_name,version,next_id,active_ids);
        last_read_active=local_frame;PHYSBAM_ASSERT(version==1);
        if(needs_destroy) for(int i=next_id;i<rigid_body_particles.Size();i++) if(!rigid_body_particles.rigid_body(i)) needs_destroy->Append(i);
        rigid_body_particles.Resize(next_id);}
    else{
        for(int id=0;id<rigid_body_particles.Size();id++) if(Is_Active(id)) active_ids.Append(id);
        next_id=rigid_body_particles.Size();}
    if(rigid_body_particles.rigid_body.Subset(active_ids).Contains(0)){ // don't need to re-read these things if we will not be initializing any newly-active bodies
        std::string key_file_list=STRING_UTILITIES::string_sprintf("%s/common/rigid_body_key_list",directory.c_str());
        std::string key_file=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_key",directory.c_str(),frame);
        char version;
        if(FILE_UTILITIES::File_Exists(key_file_list)){
            if(!frame_list_key){frame_list_key=new ARRAY<int>;FILE_UTILITIES::Read_From_File(stream_type,key_file_list,*frame_list_key);}
            local_frame=(*frame_list_key)(frame_list_key->Binary_Search(frame));}
        if(last_read_key!=local_frame && FILE_UTILITIES::File_Exists(key_file)){
            FILE_UTILITIES::Read_From_File(stream_type,key_file,version,rigid_body_particles.structure_ids);
            last_read_active=local_frame;PHYSBAM_ASSERT(version==2 || version==3);}

        try{
            std::istream* input=FILE_UTILITIES::Safe_Open_Input(directory+"/common/rigid_body_names",false);
            int num;*input>>num;input->ignore(INT_MAX,'\n');
            rigid_body_names.Resize(num);
            for(int i=0;i<rigid_body_names.Size();i++) std::getline(*input,rigid_body_names(i));
            delete input;}
        catch(FILESYSTEM_ERROR&){
            LOG::cerr<<"Did not find rigid body names."<<std::endl;
            rigid_body_names.Clean_Memory();}}

    ARRAY<bool> exists(next_id);exists.Subset(active_ids).Fill(true);
    for(int i=0;i<exists.m;i++) if(!exists(i) && Is_Active(i)) Deactivate_Body(i);
    if(active_ids.m>0){
        for(int i=0;i<active_ids.m;i++){int p=active_ids(i);
            // initialize new rigid body with given id, and initialize geometry
            if(!rigid_body_particles.rigid_body(p)){
                if(needs_init) needs_init->Append(p);
                RIGID_BODY<TV>* rigid_body=New_Body(p);
                if(p<rigid_body_names.Size()) rigid_body->name=rigid_body_names(p);
                for(int s=0;s<rigid_body_particles.structure_ids(p).m;s++)
                    if(rigid_body_particles.structure_ids(p)(s)>=0)
                        rigid_body->Add_Structure(*structure_list.Element(rigid_body_particles.structure_ids(p)(s)));}
            if(Is_Active(p) && rigid_body_particles.structure_ids(p)==VECTOR<int,3>(-1,-1,-1)) Deactivate_Body(p);}}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const
{
    articulated_rigid_body.Write(stream_type,directory,frame);
    if(stream_type.use_doubles) structure_list.template Write<double>(directory,"rigid_body_structure_",frame);
    else structure_list.template Write<float>(directory,"rigid_body_structure_",frame);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_particles",directory.c_str(),frame),rigid_body_particles);

    // update names
    rigid_body_names.Resize(rigid_body_particles.Size());
    ARRAY<int> active_ids;
    for(int id=0;id<rigid_body_particles.Size();id++) if(Is_Active(id)){active_ids.Append(id);rigid_body_names(id)=Rigid_Body(id).name;}
    if(active_ids.m>0 && !(check_stale && is_stale_active)){
        if(check_stale){
            if(!frame_list_active) frame_list_active=new ARRAY<int>;frame_list_active->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/rigid_body_active_ids_list",directory.c_str()),*frame_list_active);
            is_stale_active=false;}
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_active_ids",directory.c_str(),frame),(char)1,rigid_body_particles.Size(),active_ids);}
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(directory+"/common/rigid_body_names",false);
    *output<<rigid_body_names.Size()<<std::endl;
    for(int i=0;i<rigid_body_names.Size();i++) *output<<rigid_body_names(i)<<std::endl;
    delete output;
    char version=3;
    if(rigid_body_particles.structure_ids.m>0 && !(check_stale && is_stale_key)){
        if(check_stale){
            if(!frame_list_key) frame_list_key=new ARRAY<int>;frame_list_key->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/rigid_body_key_list",directory.c_str()),*frame_list_key);
            is_stale_key=false;}
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_key",directory.c_str(),frame),version,rigid_body_particles.structure_ids);}
}
//#####################################################################
template class RIGID_BODY_COLLECTION<VECTOR<float,1> >;
template class RIGID_BODY_COLLECTION<VECTOR<float,2> >;
template class RIGID_BODY_COLLECTION<VECTOR<float,3> >;
template class RIGID_BODY_COLLECTION<VECTOR<double,1> >;
template class RIGID_BODY_COLLECTION<VECTOR<double,2> >;
template class RIGID_BODY_COLLECTION<VECTOR<double,3> >;
}
