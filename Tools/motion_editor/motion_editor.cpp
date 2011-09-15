//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/BASIC_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_ARCBALL.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MOUSE_HANDLER.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRANSLATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Dynamics/Motion/LINEAR_BLEND_SKINNING.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Dynamics_Components/OPENGL_COMPONENT_BODY_MOTION_SEQUENCE.h>
using namespace PhysBAM;

#define ROTATE 1
#define TRANSLATE 2
#define SCALE 3

template<class T>
class ATTACHMENT_VISUALIZATION:public BASIC_VISUALIZATION,public OPENGL_MOUSE_HANDLER
{
    using BASIC_VISUALIZATION::opengl_world;

    typedef VECTOR<T,3> TV;
public:
    std::string filename;
    std::string rigid_body_filename;
    std::string mesh_filename;
    PARTICLES<TV> particles;
    TRIANGLE_MESH human_tri_mesh;
    TETRAHEDRON_MESH human_tet_mesh;
    TRIANGULATED_SURFACE<T> human_surface;
    TETRAHEDRALIZED_VOLUME<T> human_volume;
    BODY_MOTION_SEQUENCE<T> *body_motion;
    LINEAR_BLEND_SKINNING<T> *blend_skinning;
    RIGID_BODY_COLLECTION<TV> *rigid_body_collection;
    int saved_id;
    OPENGL_TRIANGULATED_SURFACE<T> *opengl_triangulated_surface;
    OPENGL_TETRAHEDRALIZED_VOLUME<T> *opengl_tetrahedralized_volume;
    OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T> *opengl_body_motion_component;
    OPENGL_COMPONENT_BASIC<OPENGL_TRIANGULATED_SURFACE<T> > *opengl_triangulated_surface_component;
    OPENGL_COMPONENT_BASIC<OPENGL_TETRAHEDRALIZED_VOLUME<T> > *opengl_tetrahedralized_volume_component;
    OPENGL_COMPONENT_BASIC<OPENGL_ARCBALL<T> > *opengl_arcball_component;
    OPENGL_COMPONENT_BASIC<OPENGL_TRANSLATION<T> > *opengl_translation_component,*opengl_scale_component;
    OPENGL_ARCBALL<T> *rotation;
    OPENGL_TRANSLATION<T> *translation,*scale;
    bool paint_mode;
    bool dragging;
    bool draw_rot,draw_trans,draw_scale;
    T length,old_mag,motion_frame_rate;
    std::string motion_output;

    int Detect_Rigid_Intersection(int x,int y)
    {
        RAY<TV> ray=opengl_world.Ray_Through_Normalized_Image_Coordinate(opengl_world.Convert_Mouse_Coordinates(x,y));
        for(int id(1);id<=rigid_body_collection->rigid_body_particle.array_collection->Size();id++) if(rigid_body_collection->Rigid_Body(id).Simplex_Intersection(ray)) return id;
        return int(0);
    }

    void Update_OpenGL_Strings()
    {
        opengl_world.Clear_Strings();
        std::ostringstream output_stream;
        output_stream<<"Paint the Fence: "<<(paint_mode?"On":"Off")<<std::endl;
        if(saved_id!=int(0)){
            output_stream<<"Selected Bone: "<<body_motion->names(opengl_body_motion_component->id_to_index.Get(saved_id))<<std::endl;
            output_stream<<"Index: "<<opengl_body_motion_component->id_to_index.Get(saved_id)<<std::endl;}
        output_stream<<"Registration Frame: "<<body_motion->saved_frame<<std::endl;
        opengl_world.Add_String(output_stream.str());
    }

    void Update_OpenGL_Current_Frame()
    {
        opengl_body_motion_component->Reinitialize(true);
    }

    void Handle_Click(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed)
    {
        if(state==GLUT_DOWN && paint_mode){dragging=true;
            if(draw_rot) rotation->Begin_Drag(opengl_world.Convert_Mouse_Coordinates(x,y));
            if(draw_scale) scale->Begin_Drag(opengl_world.Convert_Mouse_Coordinates(x,y));
            if(draw_trans) translation->Begin_Drag(opengl_world.Convert_Mouse_Coordinates(x,y));}
        else if(state==GLUT_UP){dragging=false;
            if(draw_rot) rotation->End_Drag(opengl_world.Convert_Mouse_Coordinates(x,y));
            if(draw_scale){
                int bone_index=opengl_body_motion_component->id_to_index.Get(saved_id);
                scale->End_Drag(opengl_world.Convert_Mouse_Coordinates(x,y));scale->Reinitialize();
                body_motion->base_position(bone_index).length=body_motion->trajectories(bone_index)(opengl_body_motion_component->Get_Frame()+1).length;
                body_motion->Update_Base_Position();
                opengl_body_motion_component->ui_scale=1;}
            if(draw_trans) translation->End_Drag(opengl_world.Convert_Mouse_Coordinates(x,y));}
        if(state==GLUT_DOWN && ctrl_pressed){
            if(saved_id!=int(0)){Update_Current_Frame();}
            saved_id=Detect_Rigid_Intersection(x,y);
            if(draw_rot){draw_rot=false;opengl_arcball_component->Set_Draw(draw_rot);if(saved_id!=int(0)) Rotation();}
            if(draw_scale){draw_scale=false;opengl_scale_component->Set_Draw(draw_scale);if(saved_id!=int(0)) Scale();}
            if(draw_trans){draw_trans=false;opengl_translation_component->Set_Draw(draw_trans);if(saved_id!=int(0)) Translation();}
            Update_OpenGL_Strings();}
    }

    void Handle_Drag(int x,int y)
    {
        if(dragging&&saved_id!=int(0)){
            int body_index=opengl_body_motion_component->id_to_index.Get(saved_id);
            if(draw_rot){
                rotation->Update(opengl_world.Convert_Mouse_Coordinates(x,y));
                FRAME<TV> new_transform=FRAME<TV>(TV(),rotation->qNow);
                body_motion->ui_position(body_index).targeted_rotation=new_transform*body_motion->ui_position(body_index).rotation;
                LOG::cout<<"Rotation: "<<rotation->qNow<<std::endl;
                Update_Current_Transforms(new_transform,ROTATE);}
            if(draw_trans){
                translation->Update(opengl_world.Convert_Mouse_Coordinates(x,y));
                FRAME<TV> new_transform=FRAME<TV>(translation->qNow,ROTATION<TV>());
                body_motion->ui_position(body_index).targeted_translation=new_transform*body_motion->ui_position(body_index).translation;
                LOG::cout<<"Translation: "<<translation->qNow<<std::endl;
                Update_Current_Transforms(new_transform,TRANSLATE);}
            if(draw_scale){
                assert(body_motion->base_position(body_index).length);
                scale->Update(opengl_world.Convert_Mouse_Coordinates(x,y));
                T signed_magnitude=scale->qNow.x+scale->qNow.y+scale->qNow.z;
                T scale_factor=(T)1+signed_magnitude/body_motion->base_position(body_index).length;
                body_motion->ui_position(body_index).length=scale_factor;
                LOG::cout<<"Scale: "<<signed_magnitude<<std::endl;
                FRAME<TV> new_transform=FRAME<TV>(opengl_body_motion_component->opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(saved_id).Frame()*TV(1,0,0)*signed_magnitude,ROTATION<TV>());
                //body_motion->Update_Base_Transforms(body_index,opengl_body_motion_component->Get_Frame());
                Update_Current_Transforms(new_transform,SCALE);}
        }
    }

    ATTACHMENT_VISUALIZATION(std::string& filename_input,const std::string& rigid_body_filename_input,const std::string& mesh_filename_input)
        :filename(filename_input),rigid_body_filename(rigid_body_filename_input),mesh_filename(mesh_filename_input),human_surface(human_tri_mesh,particles),human_volume(human_tet_mesh,particles),
        blend_skinning(0),saved_id(int(0)),opengl_triangulated_surface_component(0),opengl_tetrahedralized_volume_component(0),paint_mode(false),dragging(false),draw_rot(false),draw_trans(false),
        draw_scale(false),old_mag((T)0),motion_frame_rate((T)120)
    {
    }

    void Initialize_Components_And_Key_Bindings() PHYSBAM_OVERRIDE
    {
        opengl_body_motion_component=new OPENGL_COMPONENT_BODY_MOTION_SEQUENCE<T>(filename,false,rigid_body_filename,5);
        opengl_body_motion_component->Set_Draw(true);
        Add_Component(opengl_body_motion_component,"Body Motion Sequence",'\0',BASIC_VISUALIZATION::OWNED);

        FRAME<TV> trans,rot;
        //trans.t=TV(0,-1*opengl_body_motion_component->opengl_component_rigid_bodies.rigid_body_list(int(1))->axis_aligned_bounding_box.min_corner.y,0);
        trans.t=TV(-1*opengl_body_motion_component->opengl_component_rigid_body_collection.rigid_body_collection.Rigid_Body(1).axis_aligned_bounding_box.min_corner.x,0,0);
        //rot.r=ROTATION<TV>((T)pi/2.,TV(0,0,-1));
        opengl_body_motion_component->rigid_body_base_transform=rot*trans;
        opengl_body_motion_component->Reinitialize(true);

        if(mesh_filename.find(".tri")!=std::string::npos){
            FILE_UTILITIES::Read_From_File<T>(mesh_filename,human_surface);
            human_surface.Discard_Valence_Zero_Particles_And_Renumber();
            //FILE_UTILITIES::Write_To_File<T>(mesh_input+"_2",human_surface);
            opengl_triangulated_surface=new OPENGL_TRIANGULATED_SURFACE<T>(human_surface);
            human_surface.Rescale(1.04);
            opengl_triangulated_surface_component=new OPENGL_COMPONENT_BASIC<OPENGL_TRIANGULATED_SURFACE<T> >(*opengl_triangulated_surface);
            //opengl_triangulated_surface_component->Set_Draw(true);
            Add_Component(opengl_triangulated_surface_component,"tri surface",',',BASIC_VISUALIZATION::OWNED);}
        else if(mesh_filename.find_last_of(".tet")!=std::string::npos){ 
            FILE_UTILITIES::Read_From_File<T>(mesh_filename,human_volume);
            opengl_tetrahedralized_volume=new OPENGL_TETRAHEDRALIZED_VOLUME<T>(&human_volume.mesh,&human_volume.particles,OPENGL_MATERIAL::Metal(OPENGL_COLOR::Magenta(1,1)),true);
            opengl_tetrahedralized_volume_component=new OPENGL_COMPONENT_BASIC<OPENGL_TETRAHEDRALIZED_VOLUME<T> >(*opengl_tetrahedralized_volume);
            Add_Component(opengl_tetrahedralized_volume_component,"tet volume",',',BASIC_VISUALIZATION::OWNED);}
        
        scale=new OPENGL_TRANSLATION<T>();
        opengl_scale_component=new OPENGL_COMPONENT_BASIC<OPENGL_TRANSLATION<T> >(*scale);
        scale=&opengl_scale_component->object;
        scale->translate_opengl=false;
        Add_Component(opengl_scale_component,"Scale",'\0',BASIC_VISUALIZATION::OWNED);
        opengl_scale_component->Set_Draw(draw_scale);
        
        translation=new OPENGL_TRANSLATION<T>();
        opengl_translation_component=new OPENGL_COMPONENT_BASIC<OPENGL_TRANSLATION<T> >(*translation);
        translation=&opengl_translation_component->object;
        Add_Component(opengl_translation_component,"Translation",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN );
        opengl_translation_component->Set_Draw(draw_trans);

        rotation=new OPENGL_ARCBALL<T>(opengl_world);
        rotation->use_sphere_center=true;rotation->use_object_space=true;
        opengl_arcball_component=new OPENGL_COMPONENT_BASIC<OPENGL_ARCBALL<T> >(*rotation);
        rotation=&opengl_arcball_component->object;
        
        Add_Component(opengl_arcball_component,"Rotation",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN );
        opengl_arcball_component->Set_Draw(draw_rot);
    
        motion_output=filename.replace(filename.find(".mtn"),4,"_adjusted.mtn");

        body_motion=&opengl_body_motion_component->body_motion;
        rigid_body_collection=&opengl_body_motion_component->opengl_component_rigid_body_collection.rigid_body_collection;

        opengl_world.Set_Key_Binding_Category("Motion Editor");
        opengl_world.Append_Bind_Key('p',Toggle_Paint_Mode_CB("Toggle Paint Mode"));
        opengl_world.Append_Bind_Key('s',Set_Registration_Frame_CB("Set Registration Frame"));
        opengl_world.Append_Bind_Key('u',Update_All_Frames_CB("Update All Frames"));
        opengl_world.Append_Bind_Key('g',Goto_Frame_CB("Goto Frame"));
        opengl_world.Append_Bind_Key('r',Rotation_CB("Rotation"));
        opengl_world.Append_Bind_Key('t',Translation_CB("Translation"));
        opengl_world.Append_Bind_Key('e',Scale_CB("Scale"));
        opengl_world.Append_Bind_Key('z',Undo_CB("Undo"));
        opengl_world.Append_Bind_Key('l',Linear_Blend_Skinning_CB("Linear Blend Skinning"));
        opengl_world.Append_Bind_Key(GLUT_KEY_LEFT,Previous_Frame_CB("Previous Frame"));
        opengl_world.Append_Bind_Key(GLUT_KEY_RIGHT,Next_Frame_CB("Next Frame"));
        opengl_world.Append_Bind_Key('w',Write_Motion_Data_CB("Write Motion Data"));
        opengl_world.Set_Key_Binding_Category("Mesh Manipulation");
    }

    ~ATTACHMENT_VISUALIZATION()
    {}

    void Add_Key_Bindings()
    {
    if(opengl_triangulated_surface_component) opengl_world.Append_Bind_Key(',',opengl_triangulated_surface_component->Toggle_Draw_CB());
    if(opengl_tetrahedralized_volume_component) opengl_world.Append_Bind_Key(',',opengl_tetrahedralized_volume_component->Toggle_Draw_CB());}

    void Update_Skinning(int frame)
    {blend_skinning->Update_Particles(particles.X,(T)opengl_body_motion_component->Get_Frame()/motion_frame_rate);}

    void Reinitialize()
    {if(saved_id==int(0)) return;if(draw_rot) Rotation();if(draw_trans) Translation();if(draw_scale) Scale();}

    void Set_Registration_Frame()
    {opengl_body_motion_component->Save_Frame();Update_OpenGL_Strings();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Set_Registration_Frame);

    void Rotation()
    {if(saved_id==int(0)||draw_rot) return;
    int bone_index=opengl_body_motion_component->id_to_index.Get(saved_id);
    if(draw_trans){
        draw_trans=false;opengl_translation_component->Set_Draw(draw_trans);
        body_motion->Update_Base_Position();
        body_motion->ui_position(bone_index).translation=body_motion->ui_position(bone_index).targeted_translation;}
    if(draw_scale){
        draw_scale=false;opengl_scale_component->Set_Draw(draw_scale);
        body_motion->Update_Base_Position();
        body_motion->ui_position(bone_index).translation=body_motion->ui_position(bone_index).targeted_translation;}
    body_motion->Update_Transforms(bone_index,opengl_body_motion_component->Get_Frame());
    RIGID_BODY<TV> *rigid_body=&rigid_body_collection->Rigid_Body(saved_id);
    rigid_body->Update_Bounding_Box();
    rotation->Reinitialize();
    rotation->sphere=SPHERE<TV>(rigid_body->Oriented_Bounding_Box().corner,(T).03);
    draw_rot=true;opengl_arcball_component->Set_Draw(draw_rot);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Rotation);

    void Translation()
    {if(saved_id==int(0)||draw_trans) return;
    if(draw_rot){
        int bone_index=opengl_body_motion_component->id_to_index.Get(saved_id);
        draw_rot=false;opengl_arcball_component->Set_Draw(draw_rot);
        body_motion->Update_Base_Position();
        body_motion->ui_position(bone_index).rotation=body_motion->ui_position(bone_index).targeted_rotation;}
    if(draw_scale){
        int bone_index=opengl_body_motion_component->id_to_index.Get(saved_id);
        draw_scale=false;opengl_scale_component->Set_Draw(draw_scale);
        body_motion->Update_Base_Position();
        body_motion->ui_position(bone_index).translation=body_motion->ui_position(bone_index).targeted_translation;}
    RIGID_BODY<TV> *rigid_body=&rigid_body_collection->Rigid_Body(saved_id);
    rigid_body->Update_Bounding_Box();
    translation->Reinitialize();
    translation->Initialize_Ray(rigid_body->Oriented_Bounding_Box().corner,(T).03);
    draw_trans=true;opengl_translation_component->Set_Draw(draw_trans);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Translation);

    void Scale()
    {if(saved_id==int(0)||draw_scale) return;
    int bone_index=opengl_body_motion_component->id_to_index.Get(saved_id);
    if(draw_rot){
        draw_rot=false;opengl_arcball_component->Set_Draw(draw_rot);
        body_motion->Update_Base_Position();
        body_motion->ui_position(bone_index).rotation=body_motion->ui_position(bone_index).targeted_rotation;}
    if(draw_trans){
        draw_trans=false;opengl_translation_component->Set_Draw(draw_trans);
        body_motion->Update_Base_Position();
        body_motion->ui_position(bone_index).translation=body_motion->ui_position(bone_index).targeted_translation;}
    RIGID_BODY<TV> *rigid_body=&rigid_body_collection->Rigid_Body(saved_id);
    rigid_body->Update_Bounding_Box();
    body_motion->base_position(bone_index).length=body_motion->trajectories(bone_index)(opengl_body_motion_component->Get_Frame()+1).length;
    scale->Reinitialize();
    scale->Initialize_Ray(rigid_body->Oriented_Bounding_Box().corner,(T).03);
    draw_scale=true;opengl_scale_component->Set_Draw(draw_scale);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Scale);

    void Undo()
    {int body_index=opengl_body_motion_component->id_to_index.Get(saved_id);
    body_motion->ui_position(body_index).targeted_rotation=FRAME<TV>();
    body_motion->ui_position(body_index).targeted_translation=FRAME<TV>();
    FRAME<TV> frame;Update_Current_Transforms(frame,draw_trans?TRANSLATE:(draw_rot?ROTATE:SCALE));}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Undo);

    void Toggle_Paint_Mode()
    {if(paint_mode){paint_mode=false;opengl_world.Set_External_Mouse_Handler();}
    else{paint_mode=true;opengl_world.Set_External_Mouse_Handler(this);}
    Update_OpenGL_Strings();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Toggle_Paint_Mode);

    void Goto_Frame_Response()
    {if(!opengl_world.prompt_response.empty()){
        int new_frame;std::stringstream sstream(opengl_world.prompt_response);sstream>>new_frame;
        if(new_frame>=body_motion->time_grid.domain.min_corner.x&&new_frame<=body_motion->time_grid.domain.max_corner.x) opengl_body_motion_component->Set_Frame(new_frame);
        else LOG::cout<<new_frame<<" is not a valid frame.";}}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Goto_Frame_Response);

    void Goto_Frame()
    {Update_All_Frames();Reinitialize();opengl_world.Prompt_User("Enter Frame Number: ",Goto_Frame_Response_CB());}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Goto_Frame);

    void Previous_Frame()
    {Update_All_Frames();Reinitialize();opengl_body_motion_component->Prev_Frame();if(blend_skinning) Update_Skinning(opengl_body_motion_component->Get_Frame());}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Previous_Frame);
    
    void Next_Frame()
    {Update_All_Frames();Reinitialize();opengl_body_motion_component->Next_Frame();if(blend_skinning) Update_Skinning(opengl_body_motion_component->Get_Frame());}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Next_Frame);
    
    void Update_All_Frames()
    {int frame=opengl_body_motion_component->Get_Frame();
    Update_Current_Frame();
    body_motion->Propogate_Transforms(frame);
    body_motion->Propogate_Lengths(frame);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Update_All_Frames); 

    void Update_Current_Frame()
    {if(saved_id!=int(0)){body_motion->Update_Trajectories(opengl_body_motion_component->id_to_index.Get(saved_id));body_motion->Update_Base_Position();}}

    void Update_Current_Transforms(FRAME<TV> &frame,int type)
    {if(type==SCALE) body_motion->Update_Lengths(opengl_body_motion_component->id_to_index.Get(saved_id),opengl_body_motion_component->Get_Frame());
    else body_motion->Update_Targeted_Transforms(opengl_body_motion_component->id_to_index.Get(saved_id),opengl_body_motion_component->Get_Frame());
    if(type==ROTATE) body_motion->Update_Children_Rotation(opengl_body_motion_component->id_to_index.Get(saved_id),opengl_body_motion_component->Get_Frame(),frame);
    else if(type==TRANSLATE) body_motion->Update_Children_Translation(opengl_body_motion_component->id_to_index.Get(saved_id),opengl_body_motion_component->Get_Frame(),frame);
    else body_motion->Update_Children_Scale(opengl_body_motion_component->id_to_index.Get(saved_id),opengl_body_motion_component->Get_Frame(),frame);
    Update_OpenGL_Current_Frame();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Update_Current_Frame);

    void Linear_Blend_Skinning()
    {if(!opengl_triangulated_surface_component) return;
    blend_skinning=new LINEAR_BLEND_SKINNING<T>(human_surface,opengl_body_motion_component->body_motion,motion_frame_rate);
    blend_skinning->Calculate_Weights(0,1);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Linear_Blend_Skinning);

    void Write_Motion_Data()
    {Update_All_Frames();
    for(int i=1;i<=body_motion->trajectories.m;i++) for(int j=1;j<=body_motion->trajectories(i).counts.x;j++)
        body_motion->trajectories(i)(j).transform=body_motion->trajectories(i)(j).targeted_translation*body_motion->trajectories(i)(j).targeted_rotation;
    FILE_UTILITIES::Write_To_File<T>(motion_output,*body_motion);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Write_Motion_Data);
};
//#####################################################################
// Function main
//#####################################################################
int main(int argc,char *argv[])
{
    std::string htr,rigid,human;
    if(argc<4) {
        LOG::cout << std::endl << "Usage: ./motion_editor htr_file rigid_body_mesh [human_mesh]" << std::endl;
        exit(1);
    }
    htr=argv[1];rigid=argv[2];human=argv[3];
    ATTACHMENT_VISUALIZATION<float> visualization(htr,rigid,human);
    // This going to mess with the args?
    int dummy_argc=1;const char *dummy_argv[]={"paint"};
    visualization.Initialize_And_Run(dummy_argc,(char**)dummy_argv);

    return 0;
}
