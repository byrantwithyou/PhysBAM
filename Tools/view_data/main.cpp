#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST_CONTAINER.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE_MANAGER.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_HEIGHTFIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace PhysBAM;
using namespace std;

template<class T>
struct GRID_INFO_3D
{
    std::string name;
    GRID_3D<T> grid;
    GRID_3D<T> mac_grid;
    OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> > *grid_component;
    OPENGL_SLICE_MANAGER slice_manager;

    GRID_INFO_3D(OPENGL_WORLD &world)
        :grid_component(0), slice_manager(world)
    {}
};

template<class T>
struct GRID_INFO_2D
{
    std::string name;
    GRID_2D<T> grid;
    GRID_2D<T> mac_grid;

    GRID_INFO_2D(OPENGL_WORLD &world)
    {}
};

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

    void Add_Grid_2D(PARAMETER_LIST &parameter_list);
    void Add_Grid_3D(PARAMETER_LIST &parameter_list);
    void Add_Triangulated_Surface(PARAMETER_LIST &parameter_list);
    void Add_Levelset_3D(PARAMETER_LIST &parameter_list);
    void Add_Heightfield_2D(PARAMETER_LIST &parameter_list);
    template<class T2> void Add_Scalar_Field_3D(PARAMETER_LIST &parameter_list, bool mac_grid);
    template<class T2> void Add_Face_Scalar_Field_3D(PARAMETER_LIST &parameter_list);

private:
    GRID_INFO_3D<T> *Find_Grid_Info_3D_By_Name(const std::string &grid_name, const std::string &referenced_by="") const;
    GRID_INFO_2D<T> *Find_Grid_Info_2D_By_Name(const std::string &grid_name, const std::string &referenced_by="") const;

    ARRAY<GRID_INFO_3D<T>*> grid_info_list_3d;
    ARRAY<GRID_INFO_2D<T>*> grid_info_list_2d;

    std::string     parameter_list_container_filename;
};

// ------------------------------------------------------------------

template<class T> VISUALIZATION<T>::
VISUALIZATION()
    : ANIMATED_VISUALIZATION()
{
}

template<class T> VISUALIZATION<T>::
~VISUALIZATION()
{
}

template<class T> void VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Set_Extra_Arguments(1, "<param_file>");
}

template<class T> void VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);

    if (parse_args.Num_Extra_Args() != 1) parse_args.Print_Usage(true);
    else parameter_list_container_filename = parse_args.Extra_Arg(1);
}

template<class T> void VISUALIZATION<T>::
Initialize_Components()
{
    ANIMATED_VISUALIZATION::Initialize_Components();

    std::string identifier;
    PARAMETER_LIST parameter_list;
    PARAMETER_LIST_CONTAINER container(parameter_list_container_filename);

    while(container.Get_Next_Parameter_List(identifier,parameter_list))
    {
        if(identifier == "Grid_2D") Add_Grid_2D(parameter_list);
        else if(identifier == "Grid_3D") Add_Grid_3D(parameter_list);
        else if(identifier == "Triangulated_Surface") Add_Triangulated_Surface(parameter_list);
        else if(identifier == "Levelset_3D") Add_Levelset_3D(parameter_list);
        else if(identifier == "Heightfield_2D") Add_Heightfield_2D(parameter_list);
        else if(identifier == "Center_Scalar_Field_3D") {
            std::string value_type = parameter_list.template Get_Parameter<std::string>("value_type","T");
            if(value_type == "T") Add_Scalar_Field_3D<T>(parameter_list,true);
            else if(value_type == "bool") Add_Scalar_Field_3D<bool>(parameter_list,true);
            else if(value_type == "int") Add_Scalar_Field_3D<int>(parameter_list,true);
            else {
                std::cerr << "Error: Unrecognized value_type '" << value_type << "'" << std::endl;
                exit(1);
            }
        }
        else if(identifier == "Node_Scalar_Field_3D") {
            std::string value_type = parameter_list.template Get_Parameter<std::string>("value_type","T");
            if(value_type == "T") Add_Scalar_Field_3D<T>(parameter_list,false);
            else if(value_type == "bool") Add_Scalar_Field_3D<bool>(parameter_list,false);
            else if(value_type == "int") Add_Scalar_Field_3D<int>(parameter_list,false);
            else {
                std::cerr << "Error: Unrecognized value_type '" << value_type << "'" << std::endl;
                exit(1);
            }
        }
        else if(identifier == "Face_Scalar_Field_3D") {
            std::string value_type = parameter_list.template Get_Parameter<std::string>("value_type","T");
            if(value_type == "T") Add_Face_Scalar_Field_3D<T>(parameter_list);
            else if(value_type == "bool") Add_Face_Scalar_Field_3D<bool>(parameter_list);
            else if(value_type == "int") Add_Face_Scalar_Field_3D<int>(parameter_list);
            else {
                std::cerr << "Error: Unrecognized value_type '" << value_type << "'" << std::endl;
                exit(1);
            }
        }
        else {
            std::cerr << "Error: Unrecognized type '" << identifier << "'" << std::endl;
            exit(1);
        }
    }

    // initialize selection priority (highest on top)
    Selection_Priority(OPENGL_SELECTION::COMPONENT_PARTICLES_3D)=100;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX)=90;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT)=89;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE)=88;
    Selection_Priority(OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX)=85;
    Selection_Priority(OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON)=84;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D)=80;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_DEFORMABLE_OBJECT_LIST_3D)=80;
    Selection_Priority(OPENGL_SELECTION::GRID_NODE_3D)=70;
    Selection_Priority(OPENGL_SELECTION::OCTREE_NODE)=70;
    Selection_Priority(OPENGL_SELECTION::OCTREE_FACE)=65;
    Selection_Priority(OPENGL_SELECTION::GRID_CELL_3D)=60;
    Selection_Priority(OPENGL_SELECTION::OCTREE_CELL)=60;
}

template<class T> void VISUALIZATION<T>::
Add_Grid_2D(PARAMETER_LIST &parameter_list)
{
    if(!parameter_list.Is_Defined("name")) {
        std::cerr << "Grid_2D requires a 'name' parameter" << std::endl;
        exit(1);
    }

    GRID_INFO_2D<T> *grid_info = new GRID_INFO_2D<T>(opengl_world);
    grid_info->name = parameter_list.template Get_Parameter<std::string>("name");
    
    if(parameter_list.Is_Defined("filename")) {
        std::string filename = parameter_list.template Get_Parameter<std::string>("filename");
        filename = FILE_UTILITIES::Get_Frame_Filename(filename, start_frame);
        FILE_UTILITIES::Read_From_File<T>(filename,grid_info->grid);
    } else if(parameter_list.Is_Defined("levelset_filename")) {
        std::string filename = parameter_list.template Get_Parameter<std::string>("levelset_filename");
        filename = FILE_UTILITIES::Get_Frame_Filename(filename, start_frame);
        ARRAYS<VECTOR<T,2> > phi;
        LEVELSET_2D<T> levelset(grid_info->grid, phi);
        FILE_UTILITIES::Read_From_File<T>(filename,levelset);
    } else if(parameter_list.Is_Defined("m") && parameter_list.Is_Defined("n")) {
        int m = parameter_list.template Get_Parameter<int>("m");
        int n = parameter_list.template Get_Parameter<int>("n");
        T xmin = parameter_list.template Get_Parameter<T>("xmin", 0);
        T xmax = parameter_list.template Get_Parameter<T>("xmax", 1);
        T ymin = parameter_list.template Get_Parameter<T>("ymin", 0);
        T ymax = parameter_list.template Get_Parameter<T>("ymax", 1);
        std::string type = parameter_list.template Get_Parameter<std::string>("type", "regular");
        grid_info->grid=GRID_2D<T>(m,n,xmin,xmax,ymin,ymax);
        if (type != "regular") grid_info->grid.Set_MAC_Grid();
    } else {
        std::cerr << "Error: Grid_2D requires a 'filename', 'levelset_filename', or 'm','n','xmin','xmax','ymin','ymax','type' attributes" << std::endl;
        exit(1);
    }

    grid_info->mac_grid = grid_info->grid.Get_MAC_Grid();
    std::cout << "Added Grid_2D '" << grid_info->name << "': " << grid_info->grid 
              << " (" << ((grid_info->grid.MAC_offset==0.5)?"MAC":"regular") << ")" << std::endl;

    grid_info_list_2d.Append(grid_info);
}

template<class T> void VISUALIZATION<T>::
Add_Grid_3D(PARAMETER_LIST &parameter_list)
{
    if(!parameter_list.Is_Defined("name")) {
        std::cerr << "Grid_3D requires a 'name' parameter" << std::endl;
        exit(1);
    }

    GRID_INFO_3D<T> *grid_info = new GRID_INFO_3D<T>(opengl_world);
    grid_info->name = parameter_list.template Get_Parameter<std::string>("name");
    
    if(parameter_list.Is_Defined("filename")) {
        std::string filename = parameter_list.template Get_Parameter<std::string>("filename");
        filename = FILE_UTILITIES::Get_Frame_Filename(filename, start_frame);
        FILE_UTILITIES::Read_From_File<T>(filename,grid_info->grid);
    } else if(parameter_list.Is_Defined("levelset_filename")) {
        std::string filename = parameter_list.template Get_Parameter<std::string>("levelset_filename");
        filename = FILE_UTILITIES::Get_Frame_Filename(filename, start_frame);
        ARRAYS<VECTOR<T,3> > phi;
        LEVELSET_3D<T> levelset(grid_info->grid, phi);
        FILE_UTILITIES::Read_From_File<T>(filename,levelset);
    } else if(parameter_list.Is_Defined("m") && parameter_list.Is_Defined("n") && parameter_list.Is_Defined("mn")) {
        int m = parameter_list.template Get_Parameter<int>("m");
        int n = parameter_list.template Get_Parameter<int>("n");
        int mn = parameter_list.template Get_Parameter<int>("mn");
        T xmin = parameter_list.template Get_Parameter<T>("xmin", 0);
        T xmax = parameter_list.template Get_Parameter<T>("xmax", 1);
        T ymin = parameter_list.template Get_Parameter<T>("ymin", 0);
        T ymax = parameter_list.template Get_Parameter<T>("ymax", 1);
        T zmin = parameter_list.template Get_Parameter<T>("zmin", 0);
        T zmax = parameter_list.template Get_Parameter<T>("zmax", 1);
        std::string type = parameter_list.template Get_Parameter<std::string>("type", "regular");
        grid_info->grid=GRID_3D<T>(m,n,mn,xmin,xmax,ymin,ymax,zmin,zmax);
        if (type != "regular") grid_info->grid.Set_MAC_Grid();
    } else {
        std::cerr << "Error: Grid_3D requires a 'filename', 'levelset_filename', or 'm','n','mn','xmin','xmax','ymin','ymax','zmin','zmax','type' attributes" << std::endl;
        exit(1);
    }

    grid_info->mac_grid = grid_info->grid.Get_MAC_Grid();
    std::cout << "Added Grid_3D '" << grid_info->name << "': " << grid_info->grid 
              << " (" << ((grid_info->grid.MAC_offset==0.5)?"MAC":"regular") << ")" << std::endl;

    OPENGL_UNIFORM_SLICE* slice=new OPENGL_UNIFORM_SLICE(opengl_world);slice->Initialize(GRID_3D<float>(grid_info->grid));
    grid_info->slice_manager.slice=slice;
    OPENGL_GRID_3D<T> *opengl_grid = new OPENGL_GRID_3D<T>(grid_info->grid, OPENGL_COLOR::Gray(0.5));
    grid_info->grid_component = new OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> >(*opengl_grid);
    grid_info->grid_component->Set_Draw(true);
    grid_info->grid_component->selectable = true;
    Add_Component(grid_info->grid_component, grid_info->name);
    grid_info->slice_manager.Add_Object(grid_info->grid_component);

    grid_info_list_3d.Append(grid_info);
}

template<class T> void VISUALIZATION<T>::
Add_Triangulated_Surface(PARAMETER_LIST &parameter_list)
{
    std::string name = parameter_list.template Get_Parameter<std::string>("name", "noname");

    std::string grid_name = parameter_list.template Get_Parameter<std::string>("grid", "");
    GRID_INFO_3D<T> *grid_info = 0;
    if(!grid_name.empty()) grid_info = Find_Grid_Info_3D_By_Name(grid_name, "Triangulated_Surface '" + name + "'");

    if(!parameter_list.Is_Defined("filename")) {
        std::cerr << "Error: Triangulated_Surface requires a 'filename' parameter" << std::endl;
        exit(1);
    }

    std::string filename = parameter_list.template Get_Parameter<std::string>("filename");
    if(!FILE_UTILITIES::Frame_File_Exists(filename, start_frame)) {
        std::cerr << "Skipping Triangulated_Surface '" << name << "'" << std::endl;
        return;
    }

    OPENGL_COMPONENT_TRIANGULATED_SURFACE<T> *component = new OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>(filename);

    if(parameter_list.Is_Defined("color")) {
        OPENGL_COLOR color = parameter_list.template Get_Parameter<OPENGL_COLOR>("color");
        component->opengl_triangulated_surface.Set_Front_Material(OPENGL_MATERIAL::Plastic(color));
    }

    Add_Component(component, name);
    if(grid_info) grid_info->slice_manager.Add_Object(component);
}

template<class T> void VISUALIZATION<T>::
Add_Levelset_3D(PARAMETER_LIST &parameter_list)
{
    std::string name = parameter_list.template Get_Parameter<std::string>("name", "noname");

    std::string grid_name = parameter_list.template Get_Parameter<std::string>("grid", "");
    GRID_INFO_3D<T> *grid_info = 0;
    if(!grid_name.empty()) grid_info = Find_Grid_Info_3D_By_Name(grid_name, "Levelset_3D '" + name + "'");

    if(!parameter_list.Is_Defined("filename")) {
        std::cerr << "Error: Levelset_3D requires a 'filename' parameter" << std::endl;
        exit(1);
    }

    std::string filename = parameter_list.template Get_Parameter<std::string>("filename");
    if(!FILE_UTILITIES::Frame_File_Exists(filename, start_frame)) {
        std::cerr << "Skipping Levelset_3D '" << name << "'" << std::endl;
        return;
    }

    OPENGL_COMPONENT_LEVELSET_3D<T> *component = new OPENGL_COMPONENT_LEVELSET_3D<T>(filename.c_str(),"",false,false);
    component->Toggle_Slice_Color_Mode();

    Add_Component(component, name);
    if(grid_info) grid_info->slice_manager.Add_Object(component);
}

template<class T> void VISUALIZATION<T>::
Add_Heightfield_2D(PARAMETER_LIST &parameter_list)
{
    std::string name = parameter_list.template Get_Parameter<std::string>("name", "noname");
    std::string full_name = "Heightfield_2D + '" + name + "'";

    std::string grid_name = parameter_list.template Get_Parameter<std::string>("grid", "");
    GRID_INFO_2D<T> *grid_info = 0;
    grid_info = Find_Grid_Info_2D_By_Name(grid_name, full_name);

    if(!parameter_list.Is_Defined("filename")) {
        std::cerr << "Error: "<< full_name << " requires a 'filename' parameter" << std::endl;
        exit(1);
    }

    std::string filename = parameter_list.template Get_Parameter<std::string>("filename");
    if(!FILE_UTILITIES::Frame_File_Exists(filename, start_frame)) {
        std::cerr << "Skipping " << full_name << std::endl;
        return;
    }

    OPENGL_COMPONENT_HEIGHTFIELD_2D<T> *component = new OPENGL_COMPONENT_HEIGHTFIELD_2D<T>(grid_info->grid, filename.c_str());

    if(parameter_list.Is_Defined("scale")) {
        T scale = parameter_list.template Get_Parameter<T>("scale", 1);
        component->Set_Scale(scale);
    }

    Add_Component(component, name);
}

template<class T> template<class T2> void VISUALIZATION<T>::
Add_Scalar_Field_3D(PARAMETER_LIST &parameter_list, bool mac_grid)
{
    std::string name = parameter_list.template Get_Parameter<std::string>("name", "noname");

    if(!parameter_list.Is_Defined("grid")) {
        std::cerr << "Error: Scalar_Field_3D requires a 'grid' parameter" << std::endl;
        exit(1);
    }

    std::string grid_name = parameter_list.template Get_Parameter<std::string>("grid");
    GRID_INFO_3D<T> *grid_info = Find_Grid_Info_3D_By_Name(grid_name, "Scalar_Field_3D '" + name + "'");

    if(!parameter_list.Is_Defined("filename")) {
        std::cerr << "Error: Scalar_Field_3D (Grid_Based) requires a 'filename' parameter" << std::endl;
        exit(1);
    }
    std::string filename = parameter_list.template Get_Parameter<std::string>("filename");
    if(!FILE_UTILITIES::Frame_File_Exists(filename, start_frame)) {
        std::cerr << "Skipping Scalar_Field_3D '" << name << "'" << std::endl;
        return;
    }

    std::string mode = parameter_list.template Get_Parameter<std::string>("mode", "Texture");
    typename OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_MODE draw_mode;
    if (mode == "Texture") draw_mode = OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_TEXTURE;
    else draw_mode = OPENGL_SCALAR_FIELD_3D<T,T2>::DRAW_POINTS;
    const GRID_3D<T> &grid = (mac_grid) ? grid_info->mac_grid : grid_info->grid;
    OPENGL_INDEXED_COLOR_MAP *indexed_color_map = OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map();
    indexed_color_map->Set_Index_Mode(OPENGL_INDEXED_COLOR_MAP::PERIODIC);
    OPENGL_COLOR_MAP<T2> *color_map = (OPENGL_COLOR_MAP<T2>*)(void *)indexed_color_map;
    OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2> *component = 
        new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,T2>(grid,filename,
//                    OPENGL_COLOR_RAMP<T2>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0), OPENGL_COLOR::Gray(1,1)), draw_mode);
                color_map, draw_mode);
    component->opengl_scalar_field.Update();
    Add_Component(component, name);
    grid_info->slice_manager.Add_Object(component);
}

template<class T> template<class T2> void VISUALIZATION<T>::
Add_Face_Scalar_Field_3D(PARAMETER_LIST &parameter_list)
{
    std::string name = parameter_list.template Get_Parameter<std::string>("name", "noname");

    if(!parameter_list.Is_Defined("grid")) {
        std::cerr << "Error: Scalar_Field_3D requires a 'grid' parameter" << std::endl;
        exit(1);
    }

    std::string grid_name = parameter_list.template Get_Parameter<std::string>("grid");
    GRID_INFO_3D<T> *grid_info = Find_Grid_Info_3D_By_Name(grid_name, "Face_Scalar_Field_3D '" + name + "'");
    const GRID_3D<T> &grid = grid_info->grid;

    OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2> *component = 0;
    if(parameter_list.Is_Defined("filename")) {
        std::string filename = parameter_list.template Get_Parameter<std::string>("filename");
        if(!FILE_UTILITIES::Frame_File_Exists(filename, start_frame)) {
            std::cerr << "Skipping Scalar_Field_3D '" << name << "'" << std::endl;
            return;
        }
        component = new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>(grid,filename,new OPENGL_CONSTANT_COLOR_MAP<T2>(OPENGL_COLOR::Cyan()));
    } 
    else if(parameter_list.Is_Defined("filename_x") &&
            parameter_list.Is_Defined("filename_y") &&
            parameter_list.Is_Defined("filename_z")) {
        std::string filename_x = parameter_list.template Get_Parameter<std::string>("filename_x");
        std::string filename_y = parameter_list.template Get_Parameter<std::string>("filename_y");
        std::string filename_z = parameter_list.template Get_Parameter<std::string>("filename_z");

        if(!FILE_UTILITIES::Frame_File_Exists(filename_x, start_frame) ||
           !FILE_UTILITIES::Frame_File_Exists(filename_y, start_frame) ||
           !FILE_UTILITIES::Frame_File_Exists(filename_z, start_frame)) {
            std::cerr << "Skipping Scalar_Field_3D '" << name << "'" << std::endl;
            return;
        }
        component = new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T2>(grid,filename_x,filename_y,filename_z,new OPENGL_CONSTANT_COLOR_MAP<T2>(OPENGL_COLOR::Cyan()));
    }

    if (component) {
        Add_Component(component, name);
        grid_info->slice_manager.Add_Object(component);
    }
}

template<class T> GRID_INFO_2D<T> *VISUALIZATION<T>::
Find_Grid_Info_2D_By_Name(const std::string &grid_name, const std::string &referenced_by) const
{
    GRID_INFO_2D<T> *grid_info=0;
    int index;
    for(index=0;index<grid_info_list_2d.m;index++) if(grid_info_list_2d(index)->name == grid_name) break;
    if(index == grid_info_list_2d.m+1) {
        std::cerr << "Error: Could not find Grid_2D '" << grid_name << "'";
        if(!referenced_by.empty()) std::cerr << " referenced by " << referenced_by;
        std::cerr << std::endl;
        exit(1);
    } else grid_info=grid_info_list_2d(index);
    return grid_info;
}

template<class T> GRID_INFO_3D<T> *VISUALIZATION<T>::
Find_Grid_Info_3D_By_Name(const std::string &grid_name, const std::string &referenced_by) const
{
    GRID_INFO_3D<T> *grid_info=0;
    int index;
    for(index=0;index<grid_info_list_3d.m;index++) if(grid_info_list_3d(index)->name == grid_name) break;
    if(index == grid_info_list_3d.m+1) {
        std::cerr << "Error: Could not find Grid_3D '" << grid_name << "'";
        if(!referenced_by.empty()) std::cerr << " referenced by " << referenced_by;
        std::cerr << std::endl;
        exit(1);
    } else grid_info=grid_info_list_3d(index);
    return grid_info;
}

template<class T> void VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION::Add_OpenGL_Initialization();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

template<class T> void VISUALIZATION<T>::
Add_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Add_Key_Bindings();

//    const int control=-'a'+1; // only works reliably for lowercase letters
    // TOOD: different keys for different elisce managers
    for(int i=0;i<grid_info_list_3d.m;i++) {
        opengl_world.Bind_Key('\b', grid_info_list_3d(i)->slice_manager.Toggle_Slice_Mode_CB("Toggle slice mode"));
        opengl_world.Bind_Key('\\', grid_info_list_3d(i)->slice_manager.Toggle_Slice_Axis_CB("Toggle slice axis"));
        opengl_world.Bind_Key(']', grid_info_list_3d(i)->slice_manager.Increment_Slice_CB("Increment slice"));
        opengl_world.Bind_Key('[', grid_info_list_3d(i)->slice_manager.Decrement_Slice_CB("Decrement slice"));
    }

    for(int i=1;i<=PhysBAM::min(component_list.m,10);i++) {
        char key=(i==10)?'0':'0'+i;
        opengl_world.Bind_Key(key, component_list(i)->Toggle_Draw_CB());
    }
}

template<class T> void VISUALIZATION<T>::
Update_OpenGL_Strings()
{
    ANIMATED_VISUALIZATION::Update_OpenGL_Strings();

    std::ostringstream output_stream;

    if(current_selection && current_selection->type == OPENGL_SELECTION::GRID_NODE_3D){
        VECTOR_3D<int> index = ((OPENGL_SELECTION_GRID_NODE_3D<T>*)current_selection)->index;
        GRID_3D<T>& grid=((OPENGL_GRID_3D<T>*)current_selection->object)->grid;
        output_stream << "Selected node " << index << " (" << grid.Node(index) << ")" << std::endl;
    }
    else if(current_selection && current_selection->type == OPENGL_SELECTION::GRID_CELL_3D){
        VECTOR_3D<int> index = ((OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection)->index;
        GRID_3D<T>& grid=((OPENGL_GRID_3D<T>*)current_selection->object)->grid;
        output_stream << "Selected cell " << index << " (" << grid.Center(index) << ")" << std::endl;
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
    if (type_double)    visualization = new VISUALIZATION<double>;
    else                visualization = new VISUALIZATION<float>;
    visualization->Initialize_And_Run(argc, argv);

    return 0;
}
