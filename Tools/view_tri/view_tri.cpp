#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BOX_HIERARCHY_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace PhysBAM;

template<class T> 
class VISUALIZATION : public ANIMATED_VISUALIZATION
{
public:
    VISUALIZATION();
    ~VISUALIZATION();

protected:
    virtual void Add_Arguments(PARSE_ARGS &parse_args);
    virtual void Parse_Arguments(PARSE_ARGS &parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Add_OpenGL_Initialization();
    virtual void Update_OpenGL_Strings();

    // callbacks
    void Command_Prompt();
    void Command_Prompt_Response();
    void Toggle_Highlight_Boundary();
    DEFINE_CALLBACK_CREATOR(VISUALIZATION, Command_Prompt);
    DEFINE_CALLBACK_CREATOR(VISUALIZATION, Command_Prompt_Response);
    DEFINE_CALLBACK_CREATOR(VISUALIZATION, Toggle_Highlight_Boundary);

private:
    // components
    ARRAY<OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>*> triangulated_surface_component;
    ARRAY<OPENGL_COMPONENT_BASIC<OPENGL_BOX_HIERARCHY_3D<T> >*> particle_hierarchy_component;
    ARRAY<std::string> filenames;
    ARRAY<FRAME<VECTOR<T,3> > > frames;
    bool use_rgd_frames;
};

// ------------------------------------------------------------------

template<class T> VISUALIZATION<T>::
VISUALIZATION()
    : ANIMATED_VISUALIZATION(), triangulated_surface_component(0)
{
}

template<class T> VISUALIZATION<T>::
~VISUALIZATION()
{}

template<class T> void VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-use_rgd");
    parse_args.Set_Extra_Arguments(-1, "<tri_file>");
}

template<class T> void VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);

    animation_enabled=false;
    use_rgd_frames=parse_args.Is_Value_Set("-use_rgd");
    if (!parse_args.Num_Extra_Args()) parse_args.Print_Usage(true);
    else{for(int i=1;i<=parse_args.Num_Extra_Args();i++){
        animation_enabled=animation_enabled||FILE_UTILITIES::Is_Animated(parse_args.Extra_Arg(i));
        filenames.Append(parse_args.Extra_Arg(i));}}
}

template<class T> void VISUALIZATION<T>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Initialize_Components_And_Key_Bindings();

    for(int i=0;i<filenames.m;i++){std::string filename=filenames(i);
        if (!FILE_UTILITIES::Frame_File_Exists(filename, start_frame)){
            std::cerr << "Can't open " << FILE_UTILITIES::Get_Frame_Filename(filename, start_frame) << std::endl;
            exit(1);}

        triangulated_surface_component.Append(new OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>(filename));
        frames.Append(FRAME<VECTOR<T,3> >());
        triangulated_surface_component.Last()->selectable = true;
        Add_Component(triangulated_surface_component.Last(),"Triangulated surface",'t',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE);

        if(use_rgd_frames && !FILE_UTILITIES::Is_Animated(filename)){
            std::string rgd_filename=FILE_UTILITIES::Get_Basename(filename)+".rgd";
            if(FILE_UTILITIES::File_Exists(rgd_filename)){
                RIGID_BODY_PARTICLES<VECTOR<T,3> > rigid_particle;
                RIGID_BODY<VECTOR<T,3> > rigid_body(rigid_particle);FILE_UTILITIES::Read_From_File<T>(rgd_filename,rigid_body);
                std::cout << "Got rgd frame for " << filename << " from " << rgd_filename << std::endl;
                frames.Last()=rigid_body.Frame();
                *triangulated_surface_component.Last()->opengl_triangulated_surface.frame=FRAME<VECTOR<float,3> >(rigid_body.Frame());}}
    }

#if 0
    triangulated_surface->Initialize_Particle_Hierarchy();
    OPENGL_BOX_HIERARCHY_3D<T> *opengl_particle_hierarchy=new OPENGL_BOX_HIERARCHY_3D<T>(triangulated_surface->particle_hierarchy);
    opengl_particle_hierarchy->min_height=1;
    opengl_particle_hierarchy->max_height=100;
    particle_hierarchy_component = new OPENGL_COMPONENT_BASIC<OPENGL_BOX_HIERARCHY_3D<T> >(*opengl_particle_hierarchy);
    Add_Component(particle_hierarchy_component,"particle hierarchy");
#endif

    opengl_world.Bind_Key('~', Command_Prompt_CB());
    opengl_world.Bind_Key('b', Toggle_Highlight_Boundary_CB("Toggle boundary segments"));

    // initialize selection priority
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX)=100;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT)=90;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE)=80;

}

template<class T> void VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION::Add_OpenGL_Initialization();
}

template<class T> void VISUALIZATION<T>::
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

template<class T> void VISUALIZATION<T>::
Command_Prompt_Response()
{
    if (opengl_world.prompt_response!="")
    {
        std::string command;
        std::istringstream sstream(opengl_world.prompt_response);
        sstream >> command;
        if (command == "s")
        {
            // Select
            std::string type;
            int index;
            sstream >> type >> index;
            if (type == "v") 
                Set_Current_Selection(triangulated_surface_component(1)->opengl_triangulated_surface.Get_Vertex_Selection(index));
            else if (type == "s")
                Set_Current_Selection(triangulated_surface_component(1)->opengl_triangulated_surface.Get_Segment_Selection(index));
            else if (type == "t")
                Set_Current_Selection(triangulated_surface_component(1)->opengl_triangulated_surface.Get_Triangle_Selection(index));
            Selection_Callback();
        }
        else if (command == "h")
        {
            int min_height, max_height;
            sstream >> min_height >> max_height;
            for(int i=0;i<particle_hierarchy_component.m;i++) if(particle_hierarchy_component(i))
            {
                particle_hierarchy_component(i)->object.min_height = min_height;
                particle_hierarchy_component(i)->object.max_height = max_height;
            }
        }
    }
}

template<class T> void VISUALIZATION<T>::
Command_Prompt()
{
    opengl_world.Prompt_User("Command: ", Command_Prompt_Response_CB());
}

template<class T> void VISUALIZATION<T>::
Toggle_Highlight_Boundary()
{
    for(int i=0;i<triangulated_surface_component.m;i++) if(triangulated_surface_component(i))
        triangulated_surface_component(i)->opengl_triangulated_surface.highlight_boundary = !triangulated_surface_component(i)->opengl_triangulated_surface.highlight_boundary;
}

// ==========================================================================

int main(int argc, char *argv[])
{
    bool type_double = false;   // float by default
    if (PARSE_ARGS::Find_And_Remove("-float", argc, argv))
        type_double = false;
    if (PARSE_ARGS::Find_And_Remove("-double", argc, argv))
        type_double = true;

    ANIMATED_VISUALIZATION *visualization = 0;
    if(!type_double) visualization=new VISUALIZATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else visualization=new VISUALIZATION<double>;
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    visualization->Initialize_And_Run(argc, argv);

    return 0;
}
