//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_COLLECTION
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <climits>
namespace PhysBAM{
template<class TV>
struct ALLOCATE_GEOMETRY_HELPER:public ALLOCATE_HELPER<TV>
{
    RIGID_GEOMETRY_COLLECTION<TV>& collection;
    ALLOCATE_GEOMETRY_HELPER(RIGID_GEOMETRY_COLLECTION<TV>& collection_input): collection(collection_input) {}
    RIGID_GEOMETRY<TV>* Create(int index=0) PHYSBAM_OVERRIDE {PHYSBAM_FATAL_ERROR("Missing parameter \"create_collision_geometry\" in RIGID_GEOMETRY<TV> constructor call below.");return 0;}
    //RIGID_GEOMETRY<TV>* Create(int index=0) PHYSBAM_OVERRIDE {return new RIGID_GEOMETRY<TV>(collection,index);}
    virtual ~ALLOCATE_GEOMETRY_HELPER(){}
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_GEOMETRY_COLLECTION<TV>::
RIGID_GEOMETRY_COLLECTION(RIGID_GEOMETRY_PARTICLES<TV>& particles_input,RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>* rigid_geometry_example_velocities_input,
    COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input,ALLOCATE_HELPER<TV>* allocate_helper_input)
    :particles(particles_input),collision_body_list(collision_body_list_input),structure_list(*new STRUCTURE_LIST<TV,int>),always_create_structure(false),
    structure_hash(*new HASHTABLE<std::string,int>),frame_list_key(0),frame_list_active(0),check_stale(false),last_read_key(-1),last_read_active(-1),
    allocate_helper(allocate_helper_input),rigid_geometry_example_velocities(rigid_geometry_example_velocities_input),owns_particles(false),owns_collision_body_list(false)
{
    if(!allocate_helper) allocate_helper=new ALLOCATE_GEOMETRY_HELPER<TV>(*this);
    if(!collision_body_list){
        collision_body_list=new COLLISION_GEOMETRY_COLLECTION<TV>();
        owns_collision_body_list=true;}
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_GEOMETRY_COLLECTION<TV>::
RIGID_GEOMETRY_COLLECTION(RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>* rigid_geometry_example_velocities_input,COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input,
    ALLOCATE_HELPER<TV>* allocate_helper_input)
    :particles(*new RIGID_GEOMETRY_PARTICLES<TV>()),collision_body_list(collision_body_list_input),structure_list(*new STRUCTURE_LIST<TV,int>),always_create_structure(false),
    structure_hash(*new HASHTABLE<std::string,int>),frame_list_key(0),frame_list_active(0),check_stale(false),last_read_key(-1),last_read_active(-1),
    allocate_helper(allocate_helper_input),rigid_geometry_example_velocities(rigid_geometry_example_velocities_input),owns_particles(true),owns_collision_body_list(false)
{
    if(!allocate_helper) allocate_helper=new ALLOCATE_GEOMETRY_HELPER<TV>(*this);
    if(!collision_body_list){
        collision_body_list=new COLLISION_GEOMETRY_COLLECTION<TV>();
        owns_collision_body_list=true;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_GEOMETRY_COLLECTION<TV>::
~RIGID_GEOMETRY_COLLECTION()
{
    if(owns_particles) delete &particles;
    delete &structure_hash;delete allocate_helper;delete &structure_list;
    if(owns_collision_body_list) delete collision_body_list;
}
//#####################################################################
// Function Exists
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
Exists(const int particle) const
{
    return particle>=0 && particle<particles.Size() && particles.rigid_geometry(particle);
}
//#####################################################################
// Function Exists
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
Is_Active(const int particle) const
{
    return Exists(particle) && particles.rigid_geometry(particle)->particle_index>=0;
}
//#####################################################################
// Function Update_Kinematic_Particles
//#####################################################################
template<class TV> void RIGID_GEOMETRY_COLLECTION<TV>::
Update_Kinematic_Particles()
{
    static_rigid_geometry.Remove_All();kinematic_rigid_geometry.Remove_All();
    for(int p=0;p<particles.Size();p++) if(Is_Active(p)){RIGID_GEOMETRY<TV>& rigid_geometry=Rigid_Geometry(p);
        if(rigid_geometry.is_static) static_rigid_geometry.Append(p); else kinematic_rigid_geometry.Append(p);}
}
//#####################################################################
// Function Add_Rigid_Geometry
//#####################################################################
template<class TV> int RIGID_GEOMETRY_COLLECTION<TV>::
Add_Rigid_Geometry(STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file)
{
    RIGID_GEOMETRY<TV>* rigid_geometry=new RIGID_GEOMETRY<TV>(*this,true);
    return Add_Rigid_Geometry(rigid_geometry,stream_type,basename,scaling_factor,read_simplicial_boundary,read_implicit_object,read_simplicial_interior,read_rgd_file);
}
//#####################################################################
// Function Add_Rigid_Geometry
//#####################################################################
template<class TV> int RIGID_GEOMETRY_COLLECTION<TV>::
Add_Rigid_Geometry(RIGID_GEOMETRY<TV>* rigid_geometry,STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,
    const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file)
{
    int id=rigid_geometry->particle_index;

    // structures
    ARRAY<int> structure_ids;
    TV structure_center=particles.frame(id).t;
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
    particles.structure_ids(id)=VECTOR<int,3>(-1,-1,-1);
    for(int i=0;i<structure_ids.m;i++){
        if(structure_ids(i)>=0){
            particles.structure_ids(id)(i)=structure_ids(i);
            Rigid_Geometry(id).Add_Structure(*structure_list.Element(structure_ids(i)));}}

    return id;
}
//#####################################################################
// Function Register_Analytic_Replacement_Structure
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
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
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
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
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else
            structure=STRUCTURE<TV>::template Create_From_File<double>(filename);
#endif
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
template<class TV> void RIGID_GEOMETRY_COLLECTION<TV>::
Destroy_Unreferenced_Geometry() 
{
    ARRAY<bool> referenced(structure_list.Number_Of_Active_Elements());
    for(int i=0;i<particles.Size();i++) for(int j=0;j<particles.structure_ids(i).m;j++) if(particles.structure_ids(i)(j)>=0)
        referenced(structure_list.Element_Index(particles.structure_ids(i)(j)))=true;
    for(int i=structure_list.Number_Of_Active_Elements()-1;i>=0;i--) if(!referenced(i)) structure_list.Remove_Element(structure_list.Active_Element_Id(i));
}
//#####################################################################
// Function Destroy_Unreferenced_Geometry
//#####################################################################
template<class TV> RIGID_GEOMETRY<TV>* RIGID_GEOMETRY_COLLECTION<TV>::
New_Body(int index)
{
    return allocate_helper->Create(index);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void RIGID_GEOMETRY_COLLECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame,ARRAY<int>* needs_init,ARRAY<int>* needs_destroy)
{
    if(stream_type.use_doubles) structure_list.template Read<double>(directory,"rigid_body_structure_",frame);
    else structure_list.template Read<float>(directory,"rigid_body_structure_",frame);
    ARRAY<RIGID_GEOMETRY<TV>*> bodies(particles.rigid_geometry);
    particles.rigid_geometry.Fill(0);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",directory.c_str(),frame),particles);
    while(particles.rigid_geometry.m<bodies.m) delete bodies.Pop();
    particles.rigid_geometry.Prefix(bodies.m)=bodies;

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
        if(needs_destroy) for(int i=next_id;i<particles.Size();i++) if(!particles.rigid_geometry(i)) needs_destroy->Append(i);
        particles.Resize(next_id);}
    else{
        for(int id=0;id<particles.Size();id++) if(Is_Active(id)) active_ids.Append(id);
        next_id=particles.Size();}
    if(particles.rigid_geometry.Subset(active_ids).Contains(0)){ // don't need to re-read these things if we will not be initializing any newly-active bodies
        std::string key_file_list=STRING_UTILITIES::string_sprintf("%s/common/rigid_body_key_list",directory.c_str());
        std::string key_file=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_key",directory.c_str(),frame);
        char version;
        if(FILE_UTILITIES::File_Exists(key_file_list)){
            if(!frame_list_key){frame_list_key=new ARRAY<int>;FILE_UTILITIES::Read_From_File(stream_type,key_file_list,*frame_list_key);}
            local_frame=(*frame_list_key)(frame_list_key->Binary_Search(frame));}
        if(last_read_key!=local_frame && FILE_UTILITIES::File_Exists(key_file)){
            FILE_UTILITIES::Read_From_File(stream_type,key_file,version,particles.structure_ids);
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
    for(int i=0;i<exists.m;i++) if(!exists(i) && Is_Active(i)) Deactivate_Geometry(i);
    if(active_ids.m>0){
        for(int i=0;i<active_ids.m;i++){int p=active_ids(i);
            // initialize new rigid body with given id, and initialize geometry
            if(!particles.rigid_geometry(p)){
                if(needs_init) needs_init->Append(p);
                RIGID_GEOMETRY<TV>* rigid_geometry=New_Body(p);
                if(p<rigid_body_names.Size()) rigid_geometry->name=rigid_body_names(p);
                for(int s=0;s<particles.structure_ids(p).m;s++)
                    if(particles.structure_ids(p)(s)>=0)
                        rigid_geometry->Add_Structure(*structure_list.Element(particles.structure_ids(p)(s)));}
            if(Is_Active(p) && particles.structure_ids(p)==VECTOR<int,3>(-1,-1,-1)) Deactivate_Geometry(p);}}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void RIGID_GEOMETRY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const
{
    if(stream_type.use_doubles) structure_list.template Write<double>(directory,"rigid_body_structure_",frame);
    else structure_list.template Write<float>(directory,"rigid_body_structure_",frame);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",directory.c_str(),frame),particles);

    // update names
    rigid_body_names.Resize(particles.Size());
    ARRAY<int> active_ids;
    for(int id=0;id<particles.Size();id++) if(Is_Active(id)){active_ids.Append(id);rigid_body_names(id)=Rigid_Geometry(id).name;}
    if(active_ids.m>0 && !(check_stale && is_stale_active)){
        if(check_stale){
            if(!frame_list_active) frame_list_active=new ARRAY<int>;frame_list_active->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/rigid_body_active_ids_list",directory.c_str()),*frame_list_active);
            is_stale_active=false;}
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_active_ids",directory.c_str(),frame),(char)1,particles.Size(),active_ids);}
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(directory+"/common/rigid_body_names",false);
    *output<<rigid_body_names.Size()<<std::endl;
    for(int i=0;i<rigid_body_names.Size();i++) *output<<rigid_body_names(i)<<std::endl;
    delete output;
    char version=3;
    if(particles.structure_ids.m>0 && !(check_stale && is_stale_key)){
        if(check_stale){
            if(!frame_list_key) frame_list_key=new ARRAY<int>;frame_list_key->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/rigid_body_key_list",directory.c_str()),*frame_list_key);
            is_stale_key=false;}
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_key",directory.c_str(),frame),version,particles.structure_ids);}
}
//#####################################################################
template class RIGID_GEOMETRY_COLLECTION<VECTOR<float,1> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<float,2> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GEOMETRY_COLLECTION<VECTOR<double,1> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<double,2> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<double,3> >;
#endif
}
