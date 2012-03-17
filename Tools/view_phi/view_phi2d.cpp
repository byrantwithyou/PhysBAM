#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_2D.h>
#include <fstream>
#include <sstream>
#include <string>

using namespace PhysBAM;
using namespace std;

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

private:
    bool Read_Grid();

    OPENGL_COMPONENT_LEVELSET_2D<T> *levelset_component;
    OPENGL_COMPONENT_BASIC<OPENGL_GRID_2D<T> > *grid_component;

    GRID_2D<T>  grid;
    std::string levelset_filename;
};

static void *the_visualization = 0;

// ------------------------------------------------------------------

template<class T> VISUALIZATION<T>::
VISUALIZATION()
    : ANIMATED_VISUALIZATION(), levelset_component(0), grid_component(0)
{
    the_visualization = this;
    add_axes = false;
}

template<class T> VISUALIZATION<T>::
~VISUALIZATION()
{
    delete levelset_component;
    delete &grid_component->object;
    delete grid_component;
}

template<class T> void VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Option_Argument("-float");    // parsed in main()
    parse_args.Add_Option_Argument("-double");

    parse_args.Set_Extra_Arguments(1, "<phi2d_file>");
}

template<class T> void VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);

    if (parse_args.Num_Extra_Args() != 1) parse_args.Print_Usage(true);
    else levelset_filename = parse_args.Extra_Arg(1);

    if (!strstr(levelset_filename.c_str(), "%d")) animation_enabled = false;
}

template<class T> bool VISUALIZATION<T>::
Read_Grid()
{
    std::string filename = FILE_UTILITIES::Get_Frame_Filename(levelset_filename, start_frame);
    ARRAYS<VECTOR<T,2> > phi;
    LEVELSET_2D<T> levelset(grid, phi);
    std::istream* levelset_file=FILE_UTILITIES::Safe_Open_Input(filename);
    if (!levelset_file || !*levelset_file) return false;
    levelset.template Read<T>(*levelset_file);
    delete levelset_file;
    std::cout << "Got " << ((grid.MAC_offset==0)?"regular":"MAC") << " grid: " << grid.m << "x" << grid.n << std::endl;
    return true;
}

template<class T> void VISUALIZATION<T>::
Initialize_Components()
{
    ANIMATED_VISUALIZATION::Initialize_Components();

    if (!FILE_UTILITIES::Frame_File_Exists(levelset_filename, start_frame))
    {
        std::cerr << "Can't open " << FILE_UTILITIES::Get_Frame_Filename(levelset_filename, start_frame) << std::endl;
        exit(1);
    }

    if (!Read_Grid()) {
        std::cerr << "Could not read grid" << std::endl;
        exit(1);
    }

    levelset_component = new OPENGL_COMPONENT_LEVELSET_2D<T>(levelset_filename.c_str());
    levelset_component->opengl_levelset->outside_color = OPENGL_COLOR::Red(0.6);
    levelset_component->opengl_levelset->inside_color = OPENGL_COLOR::Blue(0.6);
    levelset_component->opengl_levelset->color_mode = OPENGL_LEVELSET_2D<T>::COLOR_GRADIENT;
    levelset_component->opengl_levelset->Update();
    Add_Component(levelset_component);

    OPENGL_GRID_2D<T> *opengl_grid = new OPENGL_GRID_2D<T>(grid, OPENGL_COLOR::Gray(0));
    grid_component = new OPENGL_COMPONENT_BASIC<OPENGL_GRID_2D<T> >(*opengl_grid);
    grid_component->selectable = true;
    Add_Component(grid_component);

    Selection_Priority(OPENGL_SELECTION::GRID_NODE_2D)=100;
    Selection_Priority(OPENGL_SELECTION::GRID_CELL_2D)=90;
}

template<class T> void VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION::Add_OpenGL_Initialization();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    opengl_world.Set_2D_Mode();
}

template<class T> void VISUALIZATION<T>::
Add_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Add_Key_Bindings();

    if (levelset_component)
    {
        opengl_world.Bind_Key('l', levelset_component->Toggle_Draw_CB());
        opengl_world.Bind_Key('L', levelset_component->Toggle_Color_Mode_CB());
        opengl_world.Bind_Key('`', levelset_component->Toggle_Smooth_CB());
        opengl_world.Bind_Key(';', levelset_component->Toggle_Draw_Mode_CB());
    }

    opengl_world.Bind_Key('1', grid_component->Toggle_Draw_CB());
}

template<class T> void VISUALIZATION<T>::
Update_OpenGL_Strings()
{
    ANIMATED_VISUALIZATION::Update_OpenGL_Strings();

    std::ostringstream output_stream;

    const LEVELSET_2D<T> *levelset=(levelset_component)?&levelset_component->opengl_levelset->levelset:0;

    if (current_selection && current_selection->type == OPENGL_SELECTION::GRID_CELL_2D)
    {
        VECTOR_2D<int> selected_index = ((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;

        output_stream << "Selected cell " << selected_index << " (" << grid.Center(selected_index) << ")" << std::endl;
        if (levelset && levelset->grid.MAC_offset==0.5)
            output_stream << "phi = " << levelset->phi(selected_index) << std::endl;
    }
    else if (current_selection && current_selection->type == OPENGL_SELECTION::GRID_NODE_2D)
    {
        VECTOR_2D<int> selected_index = ((OPENGL_SELECTION_GRID_NODE_2D<T>*)current_selection)->index;

        output_stream << "Selected node " << selected_index << " (" << grid.Node(selected_index) << ")" << std::endl;
        if (levelset && levelset->grid.MAC_offset==0)
            output_stream << "phi = " << levelset->phi(selected_index) << std::endl;
    }

    opengl_world.Add_String(output_stream.str().c_str());
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
    if (!type_double)    visualization = new VISUALIZATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else visualization=new VISUALIZATION<double>();
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    visualization->Initialize_And_Run(argc, argv);

    return 0;
}
