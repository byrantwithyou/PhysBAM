#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <fstream>
#include <sstream>
#include "../../Personal_Libraries/Eran_Library/PARSE_UTILS.h"
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
    virtual void Initialize_Components();
    virtual void Add_OpenGL_Initialization();
    virtual void Add_Key_Bindings();
    virtual void Update_OpenGL_Strings();

    // callbacks
    void Command_Prompt();
    void Command_Prompt_Response();
    DEFINE_CALLBACK_CREATOR(VISUALIZATION, Command_Prompt);
    DEFINE_CALLBACK_CREATOR(VISUALIZATION, Command_Prompt_Response);

private:
    // components
    TRIANGULATED_AREA<T> *triangulated_area;
    OPENGL_COMPONENT_BASIC<OPENGL_TRIANGULATED_AREA<T> > *triangulated_area_component;
    std::string     filename;
};

// ------------------------------------------------------------------

template<class T> VISUALIZATION<T>::
VISUALIZATION()
    : ANIMATED_VISUALIZATION(), triangulated_area(0), triangulated_area_component(0)
{
}

template<class T> VISUALIZATION<T>::
~VISUALIZATION()
{
    if (triangulated_area_component) delete &triangulated_area_component->object;
    delete triangulated_area_component;
    delete triangulated_area;
}

template<class T> void VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Set_Extra_Arguments(1, "<phi_file>");
}

template<class T> void VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);

    if (parse_args.Num_Extra_Args() != 1) parse_args.Print_Usage(true);
    else filename = parse_args.Extra_Arg(1);

    animation_enabled = FILE_UTILITIES::Is_Animated(filename);
}

template<class T> void VISUALIZATION<T>::
Initialize_Components()
{
    ANIMATED_VISUALIZATION::Initialize_Components();

    FILE_UTILITIES::Create_From_File<T>(filename,triangulated_area);
    triangulated_area->triangle_mesh.Initialize_Segment_Mesh();
    OPENGL_TRIANGULATED_AREA<T> *opengl_triangulated_area = new OPENGL_TRIANGULATED_AREA<T>(*triangulated_area);
    triangulated_area_component = new OPENGL_COMPONENT_BASIC<OPENGL_TRIANGULATED_AREA<T> >(*opengl_triangulated_area);
    triangulated_area_component->selectable = true;
    Add_Component(triangulated_area_component, "Triangulated area");

    // initialize selection priority
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX)=100;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT)=90;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE)=80;
}

template<class T> void VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION::Add_OpenGL_Initialization();

    opengl_world.Set_2D_Mode();
}

template<class T> void VISUALIZATION<T>::
Add_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Add_Key_Bindings();

    opengl_world.Bind_Key('~', Command_Prompt_CB());
}

template<class T> void VISUALIZATION<T>::
Update_OpenGL_Strings()
{
    ANIMATED_VISUALIZATION::Update_OpenGL_Strings();

    std::ostringstream buffer;

    if (current_selection && current_selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX)
    {
        OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T> *real_selection = (OPENGL_SELECTION_TRIANGULATED_AREA_VERTEX<T> *)current_selection;
        int index = real_selection->index;
        buffer << "Selected vertex " << index << std::endl;
        buffer << triangulated_area->particles.X(index) << std::endl;
    }
    else if (current_selection && current_selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT)
    {
        OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T> *real_selection = (OPENGL_SELECTION_TRIANGULATED_AREA_SEGMENT<T> *)current_selection;
        int index = real_selection->index;
        buffer << "Selected segment " << index << std::endl;
        if (triangulated_area->triangle_mesh.segment_mesh)
        {
            int node1,node2;
            triangulated_area->triangle_mesh.segment_mesh->segments.Get(index,node1,node2);
            buffer << node1 << ": " << triangulated_area->particles.X(node1) << std::endl;
            buffer << node2 << ": " << triangulated_area->particles.X(node2) << std::endl;
        }
        else
        {
            buffer << "segment_mesh not initialized!" << std::endl;
        }
    }
    else if (current_selection && current_selection->type == OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE)
    {
        OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T> *real_selection = (OPENGL_SELECTION_TRIANGULATED_AREA_TRIANGLE<T> *)current_selection;
        int index = real_selection->index;
        buffer << "Selected triangle " << index << std::endl;
        int node1,node2,node3;
        triangulated_area->triangle_mesh.triangles.Get(index,node1,node2,node3);
        buffer << node1 << ": " << triangulated_area->particles.X(node1) << std::endl;
        buffer << node2 << ": " << triangulated_area->particles.X(node2) << std::endl;
        buffer << node3 << ": " << triangulated_area->particles.X(node3) << std::endl;
    }

    opengl_world.Add_String(buffer.str().c_str());
}

template<class T> void VISUALIZATION<T>::
Command_Prompt_Response()
{
    if (OPENGL_WORLD::prompt_response && strcmp(OPENGL_WORLD::prompt_response,""))
    {
        std::string command;
        std::istringstream sstream(OPENGL_WORLD::prompt_response);
        sstream >> command;
        if (command == "s")
        {
            // Select
            std::string type;
            int index;
            sstream >> type >> index;
            if (type == "v") 
                Set_Current_Selection(triangulated_area_component->object.Get_Vertex_Selection(index));
            else if (type == "s")
                Set_Current_Selection(triangulated_area_component->object.Get_Segment_Selection(index));
            else if (type == "t")
                Set_Current_Selection(triangulated_area_component->object.Get_Triangle_Selection(index));
            Selection_Callback();
        }
    }
}

template<class T> void VISUALIZATION<T>::
Command_Prompt()
{
    OPENGL_WORLD::Prompt_User("Command: ", Command_Prompt_Response_CB());
}

// ==========================================================================

int main(int argc, char *argv[])
{
    bool type_double = false;   // float by default
    if (PARSE_UTILS::Find_And_Remove("-float", argc, argv))
        type_double = false;
    if (PARSE_UTILS::Find_And_Remove("-double", argc, argv))
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
