//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_HEIGHTFIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RGB_COLORS.h>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>

using namespace PhysBAM;
using namespace std;

template<class T>
class VISUALIZATION : public ANIMATED_VISUALIZATION
{
    typedef VECTOR<T,2> TV;
public:
    VISUALIZATION();
    ~VISUALIZATION();

protected:
    virtual void Add_Arguments(PARSE_ARGS &parse_args);
    virtual void Parse_Arguments(PARSE_ARGS &parse_args);
    virtual void Initialize_Components_And_Key_Bindings();

private:
    bool Read_Grid();

    // Components
    OPENGL_COMPONENT_HEIGHTFIELD_2D<T> *density_component;
    OPENGL_COMPONENT_HEIGHTFIELD_2D<T> *pressure_component;
    OPENGL_COMPONENT_BASIC<OPENGL_TRIANGULATED_SURFACE<T> > *world_geometry_component;

    GRID<TV> grid;

    // Options
    std::string basedir;
};

template<class T> VISUALIZATION<T>::
VISUALIZATION()
    : ANIMATED_VISUALIZATION(),density_component(0),pressure_component(0),world_geometry_component(0)
{}

template<class T> VISUALIZATION<T>::
~VISUALIZATION()
{
    delete density_component;
    delete pressure_component;
    delete world_geometry_component;
}

template<class T> void VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    basedir=".";   // default basedir
    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);
    parse_args.Set_Extra_Arguments(-1, "[<basedir>]");
}

template<class T> void VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);
    if (parse_args.Num_Extra_Args() > 1)
        parse_args.Print_Usage(true);
    else if (parse_args.Num_Extra_Args() == 1)
        basedir=parse_args.Extra_Arg(0);
    last_frame_filename = std::string(basedir)+"/common/last_frame";
    
    if (!parse_args.Is_Value_Set("-camera_script")) 
        camera_script_filename = basedir + std::string("/camera_script");
    opengl_window_title = std::string("2D Euler Visualization: ") + FILE_UTILITIES::Real_Path(basedir);
}

template<class T> bool VISUALIZATION<T>::
Read_Grid()
{
    std::string filename;
    filename=basedir+"/common/grid";
    if(!FILE_UTILITIES::File_Exists(filename)) return false;
    else{FILE_UTILITIES::Read_From_File<T>(filename,grid);return true;}
}

template<class T> void VISUALIZATION<T>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Initialize_Components_And_Key_Bindings();
    opengl_world.Unbind_Keys("12");
    if (!Read_Grid()){
        std::cerr << "Error reading grid" << std::endl;
        exit(-1);}
    int height_m_start=1, height_m_end=grid.counts.x, height_n_start=1, height_n_end=grid.counts.y;
    std::cout << "Using domain: " << height_m_start << " " << height_m_end << " " << height_n_start << " " << height_n_end << std::endl;

    std::string density_filename = basedir+"/%d/density";
    density_component = new OPENGL_COMPONENT_HEIGHTFIELD_2D<T>(grid,density_filename,"","",height_m_start,height_m_end,height_n_start,height_n_end);
    density_component->Toggle_Draw_Velocities();
    density_component->selectable=true;
    density_component->Reinitialize(true);
    Add_Component(density_component, "density",'1',0);
    std::string pressure_filename = basedir+"/%d/pressure";
    pressure_component = new OPENGL_COMPONENT_HEIGHTFIELD_2D<T>(grid,pressure_filename,"","",height_m_start,height_m_end,height_n_start,height_n_end);
    pressure_component->Toggle_Draw_Velocities();
    pressure_component->Toggle_Draw();
    pressure_component->selectable=true;
    pressure_component->Reinitialize(true);
    Add_Component(pressure_component, "pressure",'2',0);

    // Initialize priority ordering for selections
    Selection_Priority(OPENGL_SELECTION::COMPONENT_HEIGHTFIELD_2D);
}

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

    delete visualization;
    return 0;
}
