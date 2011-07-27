//#####################################################################
// Copyright 2005, Jared Go, Ranjitha Kumar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <fstream>
#include <sstream>
#include "../../Public_Library/Articulated_Rigid_Bodies/BEND_JOINT.h"
#include "ARTICULATION_VISUALIZATION.h"
#include "OPENGL_COMPONENT_JOINT_3D.h"
#include "OPENGL_COMPONENT_MUSCLE_3D.h"
#include "TOOL.h"
//#include <hash_map>
#include <string.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ARTICULATION_VISUALIZATION<T>::
ARTICULATION_VISUALIZATION()
    : ANIMATED_VISUALIZATION(), currmuscle(0), triangulated_surface_component(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> ARTICULATION_VISUALIZATION<T>::
~ARTICULATION_VISUALIZATION()
{
    triangulated_surface_component.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Add args
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);
    
    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Set_Extra_Arguments(-1, "<tri_file>");
}
//#####################################################################

// Parse Arguments
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);
}
//#####################################################################
// Load a list
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Load_List(std::string filename, ARRAY<std::string> &list)
{
    std::istream* is = FILE_UTILITIES::Safe_Open_Input(filename,false);
    while (!is->eof()) {
        std::string name; *is>>name;
        list.Append(name);  
    }
}
//#####################################################################
// Load bones from file (first step)
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Generate_Bones(std::string listfile)
{
    std::istream* is = FILE_UTILITIES::Safe_Open_Input(listfile,false);
    while (!is->eof()) 
    {
        // Load these
        std::string name; std::string file; *is>>name>>file;                       
        std::cout<< "Loading " << name << "...\n";

        // Transform file if dash
        if (file == "-") file = name;
        std::string filepath=bone_path + file;

        // Create surface
        OPENGL_COMPONENT_TRIANGULATED_SURFACE<T> *bone = new OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>(filepath + ".tri.gz");            
        bone->opengl_triangulated_surface.front_material = OPENGL_MATERIAL(OPENGL_COLOR(0.8f, 0.82f, 0.72f, 1.0f));
        //bone->selectable = true;        
        Add_Component(bone,name);

        // Add to hash based on name
        bone_map.Set(name,bone);        

        // Set component frame
        std::string rgd_filename=filepath+".rgd.gz";
        if(FILE_UTILITIES::File_Exists(rgd_filename)){

            RIGID_BODY_3D<T> rigid_body;
            FILE_UTILITIES::Read_From_File<T>(rgd_filename,rigid_body);                        
            bone->opengl_triangulated_surface.frame = new FRAME_3D<float>(rigid_body.frame);            
        }
    }
    delete is;
}
//#####################################################################
// Load muscles from file (second step)
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Generate_Muscles(std::string listfile)
{
    Load_List(listfile, muscle_list);

    for(int i=1;i<=muscle_list.m;i++) 
    {
        // Name
        std::string filename=muscle_list(i);
        std::string filepath=muscle_path + filename + ".viapts";

        std::cout<<"Loading "<<filepath<<"...\n";
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(filepath,false);

        // Read in params
        int num_points;T x,y,z;std::string bone_name,attached,segment_type;
        *input>>num_points;
        
        // Create muscle component for both sides (?)
        OPENGL_COMPONENT_MUSCLE_3D<T> *muscle = new OPENGL_COMPONENT_MUSCLE_3D<T>(this);        
        muscle->parameters.Resize(num_points);
        // Read in this information for each point
        for(int i=1;i<=num_points;i++)
        {      
            muscle->parameters(i).Resize(4);
            if (i==1) {
                *input>>x>>y>>z>>bone_name>>attached;
                segment_type = "";
                for(int p=1;p<=4;p++) muscle->parameters(i)(p)=0.0;
            }
            else {
                *input>>x>>y>>z>>bone_name>>attached>>segment_type;
                if(segment_type=="analytic") 
                    *input>>muscle->parameters(i)(1)>>muscle->parameters(i)(2)>>muscle->parameters(i)(3)>>muscle->parameters(i)(4);
            }
                
            if (!bone_map.Get(bone_name)) {
                if (!bone_map.Get(bone_name + "_right")) {
                    std::cout<<"Unable to find bone "<<bone_name<<", skipping muscle.\n";
                    continue;
                } else {
                    bone_name += "_right";
                }                
            } else { std::cout<<"Attaching to bone "<<bone_name<<".\n"; }

            OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>* surf;
            if (bone_map.Get(bone_name,surf)) {                   
                FRAME_3D<float> *frameptr = surf->opengl_triangulated_surface.frame;
                muscle->AddPointInFrame(VECTOR<T,3>(x,y,z),FRAME_3D<T>(*frameptr),bone_name,attached,segment_type);                
            }            
        }

        // Add this
        muscle->name = filename;
        muscle->selectable = true;
        muscle_map.Set(filename, muscle);
        Add_Component(muscle, filename);

        // Close stream
        delete input;
    }
}
//#####################################################################
// Load articulation points (joints) from file
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Generate_Joints(std::string listfile)
{
    std::istream* is = FILE_UTILITIES::Safe_Open_Input(listfile,false);
    int count; *is>>count;
    for (int i = 0; i < count; i++) 
    {
        std::string parent, child;
        JOINT_3D<T>* joint = Parse_Joint(*is,parent,child);
        if (joint == NULL) continue;

        OPENGL_COMPONENT_JOINT_3D<T>* jcomp = new OPENGL_COMPONENT_JOINT_3D<T>(this,joint);
        jcomp->name = joint->name;
        bone_map.Get(parent,jcomp->parent_surf);
        bone_map.Get(child,jcomp->child_surf);
        Add_Component(jcomp, joint->name);
    }    
    delete is;
}
//#####################################################################
// Parse a single joint from the file
//#####################################################################
template<class T> JOINT_3D<T>* ARTICULATION_VISUALIZATION<T>::
Parse_Joint(std::istream &s, std::string &parent, std::string &child)
{
    // Joint params
    std::string name,type,bone_child,bone_parent,joint_center_frame;
    VECTOR<T,3> center,axis;

    // Read in parameters
    std::string token;  
    s>>token; if (token != "beginjoint") { std::cerr << "Not a joint, exiting.\n"; return NULL; }
    for (s>>token; token != "endjoint"; s>>token)
    {
        if (token == "joint_name") { s>>name; }
        else if (token == "joint_type") { s>>type; }
        else if (token == "bone_child") { s>>bone_child; }
        else if (token == "bone_parent") { s>>bone_parent; }
        else if (token == "joint_center") { s>>center; }
        else if (token == "joint_frame") { s>>joint_center_frame; }
        else if (token == "axis") { s>>axis; }
        else { std::cerr << "Unable to process directive "<<token<<", exiting.\n"; }
    }

    // Create joint of appropriate type, and parse type-dependent args
    JOINT_3D<T>* joint = NULL;
    if (type == "bend") {joint = new BEND_JOINT_3D<T>();}
    else if (type == "point") {joint = new POINT_JOINT_3D<T>();}

    // Parse general args, name, etc.
    joint->name = name;

    // Reflect hack
    bone_parent += "_right";
    bone_child += "_right";

    // Set output
    parent = bone_parent;
    child = bone_child;

    // Get frames and check names
    OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>* parent_surface;
    OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>* child_surface;
    if (!bone_map.Get(bone_parent,parent_surface)) { std::cerr << "Unable to find parent "<<bone_parent; return NULL; }
    if (!bone_map.Get(bone_child,child_surface)) { std::cerr << "Unable to find parent "<<bone_child; return NULL; }

    // Get frames    
    FRAME_3D<T> parent_frame = FRAME_3D<T>(*(parent_surface->opengl_triangulated_surface.frame));
    FRAME_3D<T> child_frame = FRAME_3D<T>(*(child_surface->opengl_triangulated_surface.frame));

    // Set frames
    FRAME_3D<T> joint_frame=FRAME_3D<T>(center);
    if(type=="bend") joint_frame.r=QUATERNION<T>::Rotation_Quaternion(VECTOR<T,3>(1,0,0),axis);

    // Set parent to child transforms    
    if(joint_center_frame=="parent"){
        joint->Set_Joint_To_Parent_Frame(joint_frame);
        joint->Set_Joint_To_Child_Frame(child_frame.Inverse()*parent_frame*joint_frame);
    }
    else if(joint_center_frame=="child"){
        joint->Set_Joint_To_Parent_Frame(parent_frame.Inverse()*child_frame*joint_frame);
        joint->Set_Joint_To_Child_Frame(joint_frame);
    }
    else if(joint_center_frame=="world"){
        joint->Set_Joint_To_Parent_Frame(parent_frame.Inverse()*joint_frame);
        joint->Set_Joint_To_Child_Frame(child_frame.Inverse()*joint_frame);
    }
    else std::cerr<<"Unrecognized joint_center_frame: "<<joint_center_frame<<std::endl;

    // Return the joint
    return joint;
}
//#####################################################################
// Load articulation points (joints) from file
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Save_Bones(std::string listfile)
{    
}
//#####################################################################
// Load articulation points (joints) from file
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Save_Muscles()
{   
    // Iterate over all muscles and write modified data back out
    for(int i=1;i<=muscle_list.m;i++) 
    {
        // get name
        std::string filename=muscle_list(i);
        std::string filepath=muscle_path + filename + ".viapts";

        // save
        std::cout<<"Saving "<<filepath<<"...\n";
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filepath,false);

        // get the muscle        
        OPENGL_COMPONENT_MUSCLE_3D<T> *m=0; muscle_map.Get(filename,m);
        
        // write points now        
        *output<<m->points.m<<std::endl;

        // write list now
        for (int i = 1; i <= m->points.m; i++) 
        {
            // convert bone name back to one side hack
            std::string bname = m->bones(i);
            std::string::size_type off = bname.find("_right");

            // shorten name
            if (off != std::string::npos) { bname = bname.substr(0,off); }

            // write out
            *output<<m->points(i)<<" "<<bname<<" "<<m->attach(i)<<" "<<m->segment_type(i);
            if(m->segment_type(i)=="analytic") for(int p=1;p<=4;p++) *output<<" "<<m->parameters(i)(p);
            *output<<std::endl;
        }

        // close
        output->flush();        
        delete output;
    }  
}
//#####################################################################
// Load articulation points (joints) from file
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Save_Joints(std::string listfile)
{    
}
//#####################################################################
// Initialize_Components
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Initialize_Components()
{
    ANIMATED_VISUALIZATION::Initialize_Components();

    // Parse data directory
    bone_path = "../../Public_Data/Rigid_Bodies/New_Visible_Human_Bones/";
    muscle_path = "../../Public_Data/SIMM_Data/muscle_viapoints/";

    // Create bone opengl components
    Generate_Bones("../../Public_Data/SIMM_Data/bones.pbl.txt");
    Generate_Muscles("../../Public_Data/SIMM_Data/muscles.pml.txt");
    Generate_Joints("../../Public_Data/SIMM_Data/hand_joints");
    
    // set custom tool
    this->opengl_world.Set_External_Mouse_Handler(new TOOL<T>(opengl_world));

}
//#####################################################################
// Handle_Object_Selected
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Handle_Object_Selected(OPENGL_COMPONENT_MUSCLE_3D<T>*m,int viapt)
{
    // On selection, change to appropriate mode, set hilights, etc.
    std::cout<<"Mouse handler setting for "<<m->name<<" point "<<viapt<<".\n";
    std::cout<<"This is the SEGMENT: "<<m->segment<<" \n";
    //unhighlight previous muscle and then set current muscle to be selected and curr_muscle
    if (currmuscle != 0) { currmuscle->selected = false; }
    m->selected = true;
    currmuscle = m;
}
//#####################################################################
// Add_OpenGL_Initialization
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION::Add_OpenGL_Initialization();
}
//#####################################################################
// Add_Key_Bindings
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Add_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Add_Key_Bindings();
    opengl_world.Bind_Key('s', Save_Muscles_Command_CB("Save muscles back to disk"));
    opengl_world.Bind_Key('t', Change_Segment_Type_Command_CB("Change muscle segment type"));
    //opengl_world.Bind_Key(

}
//#####################################################################
// Update_OpenGL_Strings
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Save_Muscles_Command()
{    
    std::cout<<"Saving muscles!\n";
    Save_Muscles();
}


template<class T> bool ARTICULATION_VISUALIZATION<T>::
Valid_Segment_Type(int segment_type) {
    return (segment_type>0 && segment_type<=MUSCLE_SEGMENT<T,VECTOR<T,3> >::NUM_MUSCLE_SEGMENT_TYPE);
}
template<class T> void ARTICULATION_VISUALIZATION<T>::
Set_Segment_Type(int segment_type){
    if (segment_type==1) {
        currmuscle->segment_type(currmuscle->segment)="linear";
    }
    else if(segment_type==2) {
        currmuscle->segment_type(currmuscle->segment)="analytic";
        char buf[16];
        sprintf(buf,"%f",currmuscle->parameters(currmuscle->segment)(1));
        opengl_world.Prompt_User("Enter Curve Thickness: ", Change_Curve_Thickness_Prompt_CB(),buf);
    }
}
template<class T> void ARTICULATION_VISUALIZATION<T>::
Change_Segment_Type_Prompt()
{
    if(!opengl_world.prompt_response.empty()){
        int segment_type=atoi(opengl_world.prompt_response.c_str());
        if(Valid_Segment_Type(segment_type)) Set_Segment_Type(segment_type);}
}
template<class T> void ARTICULATION_VISUALIZATION<T>::
Change_Segment_Type_Command()
{    
    std::cout<<"Changing segment_type!\n";
    opengl_world.Prompt_User("Enter Muscle Type: ", Change_Segment_Type_Prompt_CB(),"");
}

template<class T> void ARTICULATION_VISUALIZATION<T>::
Change_Curve_Thickness_Prompt()
{
    if(!opengl_world.prompt_response.empty()){
        T curve_thickness=atof(opengl_world.prompt_response.c_str());
        std::cout<<curve_thickness<<std::endl;
        if(curve_thickness>0) currmuscle->parameters(currmuscle->segment)(1)=curve_thickness;
        char buf[16];
        sprintf(buf,"%f",currmuscle->parameters(currmuscle->segment)(2));
        opengl_world.Prompt_User("Enter Curve Offset Thickness: ", Change_Curve_Offset_Thickness_Prompt_CB(),buf);
    }
}
template<class T> void ARTICULATION_VISUALIZATION<T>::
Change_Curve_Offset_Thickness_Prompt()
{
    if(!opengl_world.prompt_response.empty()){
        T curve_offset_thickness=atof(opengl_world.prompt_response.c_str());
        std::cout<<curve_offset_thickness<<std::endl;
        if(curve_offset_thickness>0) currmuscle->parameters(currmuscle->segment)(2)=curve_offset_thickness;
        char buf[16];
        sprintf(buf,"%f",currmuscle->parameters(currmuscle->segment)(3));
        opengl_world.Prompt_User("Enter Tendon Fraction 1: ", Change_Tendon_Fraction_1_Prompt_CB(),buf);
    }
}
template<class T> void ARTICULATION_VISUALIZATION<T>::
Change_Tendon_Fraction_1_Prompt()
{
    if(!opengl_world.prompt_response.empty()){
        T tendon_fraction_1=atof(opengl_world.prompt_response.c_str());
        std::cout<<tendon_fraction_1<<std::endl;
        if(tendon_fraction_1>0) currmuscle->parameters(currmuscle->segment)(3)=tendon_fraction_1;
        char buf[16];
        sprintf(buf,"%f",currmuscle->parameters(currmuscle->segment)(4));
        opengl_world.Prompt_User("Enter Tendon Fraction 2: ", Change_Tendon_Fraction_2_Prompt_CB(),buf);
    }
}
template<class T> void ARTICULATION_VISUALIZATION<T>::
Change_Tendon_Fraction_2_Prompt()
{
    if(!opengl_world.prompt_response.empty()){
        T tendon_fraction_2=atof(opengl_world.prompt_response.c_str());
        std::cout<<tendon_fraction_2<<std::endl;
        if(tendon_fraction_2>0) currmuscle->parameters(currmuscle->segment)(4)=tendon_fraction_2;}
}

//#####################################################################
// Update_OpenGL_Strings
//#####################################################################
template<class T> void ARTICULATION_VISUALIZATION<T>::
Update_OpenGL_Strings()
{
    ANIMATED_VISUALIZATION::Update_OpenGL_Strings();

    std::ostringstream output_stream;
    if(current_selection){
        output_stream << "Selected " << current_selection->object->name << std::endl;
        current_selection->object->Print_Selection_Info(output_stream,current_selection);
    }
    opengl_world.Add_String(output_stream.str());
}
//#####################################################################
template class ARTICULATION_VISUALIZATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ARTICULATION_VISUALIZATION<double>;
#endif
