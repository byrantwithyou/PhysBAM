//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Nick Rasmussen, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include <Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEBUG_PARTICLES_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D.h>
#include <Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <Dynamics/Particles/SPH_PARTICLES.h>
#include <fstream>
#include <sstream>
#include <string.h>
using namespace PhysBAM;
using namespace std;
//#####################################################################
// OPENGL_2D_VISUALIZATION
//#####################################################################
template<class T>
class OPENGL_2D_VISUALIZATION:public ANIMATED_VISUALIZATION<T>
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
public:
    using ANIMATED_VISUALIZATION<T>::add_axes;using ANIMATED_VISUALIZATION<T>::start_frame;
    using ANIMATED_VISUALIZATION<T>::opengl_world;using ANIMATED_VISUALIZATION<T>::opengl_window_title;
    using ANIMATED_VISUALIZATION<T>::last_frame_filename;using ANIMATED_VISUALIZATION<T>::frame;
    using ANIMATED_VISUALIZATION<T>::frame_increment;using ANIMATED_VISUALIZATION<T>::frame_title;
    using ANIMATED_VISUALIZATION<T>::Update_OpenGL_Strings;using ANIMATED_VISUALIZATION<T>::camera_script_filename;
    using ANIMATED_VISUALIZATION<T>::Add_Component;using ANIMATED_VISUALIZATION<T>::Selection_Priority;
    using ANIMATED_VISUALIZATION<T>::frame_rate;using ANIMATED_VISUALIZATION<T>::draw_all_objects_cb;
    using ANIMATED_VISUALIZATION<T>::Set_Current_Selection;using ANIMATED_VISUALIZATION<T>::stream_type;
    OPENGL_2D_VISUALIZATION(STREAM_TYPE stream_type);
    ~OPENGL_2D_VISUALIZATION();

protected:
    virtual void Add_Arguments(PARSE_ARGS& parse_args);
    virtual void Parse_Arguments(PARSE_ARGS& parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Add_OpenGL_Initialization();

private:
    void Read_Grid();
    virtual void Pre_Frame_Extra();
    virtual void Set_Frame_Extra();

    void Command_Prompt();
    void Command_Prompt_Response();
    void Toggle_2d_Mode();
    void Change_Level(int new_level);

    // Components
    OPENGL_COMPONENT_SCALAR_FIELD_2D<T> *sph_cell_weight_component;
    OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* pressure_component,* coarse_pressure_component;
    OPENGL_COMPONENT_PARTICLES_2D<T>* positive_particles_component;
    OPENGL_COMPONENT_PARTICLES_2D<T>* negative_particles_component;
    OPENGL_COMPONENT_PARTICLES_2D<T>* removed_positive_particles_component;
    OPENGL_COMPONENT_PARTICLES_2D<T>* removed_negative_particles_component;
    OPENGL_COMPONENT_BASIC<T,OPENGL_GRID_2D<T> >* grid_component;
    // Options
    std::string basedir;

    GRID<TV> grid,mac_grid,regular_grid;
    bool has_valid_grid;
    bool has_valid_sub_grids;
    bool node_based;

    ARRAY<int> rigid_bodies_no_draw_list;
};
//#####################################################################
// OPENGL_2D_VISUALIZATION
//#####################################################################
template<class T> OPENGL_2D_VISUALIZATION<T>::
OPENGL_2D_VISUALIZATION(STREAM_TYPE stream_type)
    :ANIMATED_VISUALIZATION<T>(stream_type),
    sph_cell_weight_component(0),
    pressure_component(0),coarse_pressure_component(0),positive_particles_component(0),negative_particles_component(0),
    removed_positive_particles_component(0),removed_negative_particles_component(0),
    grid_component(0)
{
    add_axes=false;
}
//#####################################################################
// OPENGL_2D_VISUALIZATION
//#####################################################################
template<class T> OPENGL_2D_VISUALIZATION<T>::
~OPENGL_2D_VISUALIZATION()
{
    delete grid_component;
}
//#####################################################################
// Add_Arguments
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS& parse_args)
{
    basedir=".";
    ANIMATED_VISUALIZATION<T>::Add_Arguments(parse_args);
    parse_args.Add("-rigid_bodies_no_draw",&rigid_bodies_no_draw_list,"id","Do not draw this rigid body (may be repeated)");
    parse_args.Extra_Optional(&basedir,"basedir","base directory");
}
//#####################################################################
// Parse_Arguments
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS& parse_args)
{
#ifdef __linux__
    opengl_window_title="opengl_2d: " + FILE_UTILITIES::Real_Path(basedir);
#endif

    if(FILE_UTILITIES::File_Exists(basedir+"/common/first_frame")) FILE_UTILITIES::Read_From_Text_File(basedir+"/common/first_frame",start_frame);

    ANIMATED_VISUALIZATION<T>::Parse_Arguments(parse_args);

    if(parse_args.unclaimed_arguments)
        parse_args.Print_Usage(true);

    last_frame_filename=basedir+"/common/last_frame";

    // don't override camera script filename if it was already set in base class based on command line argument
    if(camera_script_filename.empty()) camera_script_filename=basedir+"/camera_script";
}
//#####################################################################
// Read_Grid
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Read_Grid()
{
    has_valid_grid=false;
    has_valid_sub_grids=false;
    std::string filename,coarse_filename,sub_filename;

    sub_filename=LOG::sprintf("%s/%d/sub_grids",basedir.c_str(),start_frame);
    filename=LOG::sprintf("%s/%d/levelset",basedir.c_str(),start_frame);
    // For backwards compatibility
    if(!FILE_UTILITIES::File_Exists(filename)) filename=LOG::sprintf("%s/%d/levelset.phi",basedir.c_str(),start_frame);
    if(!FILE_UTILITIES::File_Exists(filename)) filename=LOG::sprintf("%s/%d/levelset_0",basedir.c_str(),start_frame);

    if(FILE_UTILITIES::File_Exists(filename)){
        LOG::cout<<"Reading grid from '"<<filename<<"'..."<<std::endl<<std::flush;
        ARRAY<T,VECTOR<int,2> > phi;
        LEVELSET<TV> levelset(grid,phi);
        FILE_UTILITIES::Read_From_File(stream_type,filename,levelset);
        has_valid_grid=true;}
    else if(FILE_UTILITIES::File_Exists(basedir+"/common/grid")){
        filename=basedir+"/common/grid";
        LOG::cout<<"Reading grid from '"<<filename<<"'..."<<std::flush;
        FILE_UTILITIES::Read_From_File(stream_type,filename,grid);
        has_valid_grid=true;}

    if(has_valid_grid){
        node_based=!grid.Is_MAC_Grid();
        mac_grid=grid.Get_MAC_Grid();regular_grid=grid.Get_Regular_Grid();}
}

//#####################################################################
// Initialize_Components_And_Key_Bindings
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION<T>::Initialize_Components_And_Key_Bindings();
    //opengl_world.Unbind_Keys("abcdDeilLMnoOStTuvVZ 1!2@3#4$5%67&*890 ;'~\t=-`\\^[]");
    Read_Grid();
    std::string filename,filename2;
    opengl_world.Set_Key_Binding_Category_Priority(1);

    // Density
    filename=basedir+"/%d/density";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* density_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_component->opengl_scalar_field.Set_Uniform_Contour_Values(.2568,6.067,(6.067-.2568)/30.);
        density_component->opengl_scalar_field.Update();
        opengl_world.Set_Key_Binding_Category("Density");
        Add_Component(density_component,"Density",'d',BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key('D',density_component->viewer_callbacks.Get("toggle_color_map"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F8),density_component->viewer_callbacks.Get("toggle_draw_mode"));
        opengl_world.Append_Bind_Key('`',density_component->viewer_callbacks.Get("toggle_smooth"));}
    filename=basedir+"/%d/density_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* density_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_gradient_component->opengl_scalar_field.Update();
        density_gradient_component->viewer_callbacks.Get("toggle_draw").func();
        opengl_world.Set_Key_Binding_Category("Density");
        Add_Component(density_gradient_component,"Density Gradient",'8',BASIC_VISUALIZATION<T>::OWNED);}
    // Soot
    filename=basedir+"/%d/soot";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* soot_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.01,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        soot_component->opengl_scalar_field.Set_Uniform_Contour_Values(0,2,.1);
        soot_component->opengl_scalar_field.Update();
        opengl_world.Set_Key_Binding_Category("Soot");
        Add_Component(soot_component,"Soot",'i',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',soot_component->viewer_callbacks.Get("toggle_color_map"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F8),soot_component->viewer_callbacks.Get("toggle_draw_mode"));
        opengl_world.Append_Bind_Key('`',soot_component->viewer_callbacks.Get("toggle_smooth"));}
    filename=basedir+"/%d/internal_energy";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* internal_energy_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Yellow(1,1)));
        internal_energy_component->opengl_scalar_field.Update();
        internal_energy_component->viewer_callbacks.Get("toggle_draw").func();
        opengl_world.Set_Key_Binding_Category("Density");
        Add_Component(internal_energy_component,"Internal_Energy",'E',BASIC_VISUALIZATION<T>::OWNED);}

    // Temperature
    filename=basedir+"/%d/temperature";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        ARRAY<T,VECTOR<int,2> > temp0;FILE_UTILITIES::Read_From_File(stream_type,FILE_UTILITIES::Get_Frame_Filename(filename,start_frame),temp0);
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* temperature_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Red(1,1)));
        temperature_component->opengl_scalar_field.Set_Scale_Range(temp0.Min(),temp0.Max());
        opengl_world.Set_Key_Binding_Category("Temperature");
        Add_Component(temperature_component,"Temperature",'t',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('T',temperature_component->viewer_callbacks.Get("toggle_color_map"));
        opengl_world.Append_Bind_Key('`',temperature_component->viewer_callbacks.Get("toggle_smooth"));}

    // SPH
    filename=basedir+"/%d/sph_cell_weights";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        sph_cell_weight_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,100,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1)));
        sph_cell_weight_component->opengl_scalar_field.Update();
        Add_Component(sph_cell_weight_component,"SPH Cell Weights",'\0',BASIC_VISUALIZATION<T>::OWNED);}

    // Grid vorticity and rasterized vortex particles
    const T vorticity_visualization_min=-10,vorticity_visualization_max=10;
    filename=basedir+"/%d/grid_vorticity";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        opengl_world.Set_Key_Binding_Category("Vorticity");
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(vorticity_visualization_min,vorticity_visualization_max))
            ,"Grid Vorticity",'u',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/grid_vorticity_raw_particles";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        opengl_world.Set_Key_Binding_Category("Vorticity");
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(vorticity_visualization_min,vorticity_visualization_max)),
            "Grid Vorticity Raw Particles",'i',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/grid_vorticity_particles";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        opengl_world.Set_Key_Binding_Category("Vorticity");
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(vorticity_visualization_min,vorticity_visualization_max)),
            "Grid Vorticity Particles",'o',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}

    // Level sets
    filename=basedir+"/%d/levelset";
    OPENGL_COMPONENT_LEVELSET_2D<T>* levelset_component=0;
    if(!FILE_UTILITIES::Frame_File_Exists(filename,start_frame) && !FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_0",start_frame)) filename=basedir+"/%d/levelset.phi"; // for backwards compatiblity
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame) || FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_0",start_frame)){
        levelset_component=new OPENGL_COMPONENT_LEVELSET_2D<T>(stream_type,filename,basedir+"/%d/levelset_%d");
        if(levelset_component->opengl_levelsets.m>1) for(int j=0;j<levelset_component->opengl_levelsets.m;j++){
            levelset_component->opengl_levelsets(j)->draw_cells=false;
            levelset_component->opengl_levelsets(j)->draw_area=false;
            levelset_component->opengl_levelsets(j)->draw_curve=true;
            levelset_component->opengl_levelsets(j)->Update();}
        else{
            levelset_component->opengl_levelsets(0)->Set_Inside_And_Outside_Colors(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red(.5));
            levelset_component->opengl_levelsets(0)->draw_cells=false;
            levelset_component->opengl_levelsets(0)->draw_area=false;
            levelset_component->opengl_levelsets(0)->draw_curve=true;
            levelset_component->opengl_levelsets(0)->Update();}
        opengl_world.Set_Key_Binding_Category("Level Set");
        Add_Component(levelset_component,"Level Set",'l',BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key(';',levelset_component->viewer_callbacks.Get("toggle_draw_mode"));
        opengl_world.Append_Bind_Key('\'',levelset_component->viewer_callbacks.Get("toggle_draw_sign"));
        opengl_world.Append_Bind_Key('L',levelset_component->viewer_callbacks.Get("toggle_color_mode"));
        opengl_world.Append_Bind_Key('n',levelset_component->viewer_callbacks.Get("toggle_normals"));
        opengl_world.Append_Bind_Key('`',levelset_component->viewer_callbacks.Get("toggle_smooth"));
        opengl_world.Append_Bind_Key('>',levelset_component->viewer_callbacks.Get("next_set"));
        opengl_world.Append_Bind_Key('<',levelset_component->viewer_callbacks.Get("previous_set"));
        opengl_world.Append_Bind_Key('M',levelset_component->viewer_callbacks.Get("toggle_draw_multiple_levelsets"));
        opengl_world.Append_Bind_Key('^',levelset_component->viewer_callbacks.Get("toggle_draw_ghost_values"));}
    filename=basedir+"/%d/object_levelset";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_LEVELSET_2D<T>* object_levelset_component=new OPENGL_COMPONENT_LEVELSET_2D<T>(stream_type,filename);
        object_levelset_component->opengl_levelset->outside_color=OPENGL_COLOR::Red(.6);
        object_levelset_component->opengl_levelset->inside_color=OPENGL_COLOR::Blue(.6);
        object_levelset_component->opengl_levelset->draw_cells=true;
        object_levelset_component->opengl_levelset->draw_area=false;
        object_levelset_component->opengl_levelset->draw_curve=false;
        object_levelset_component->opengl_levelset->Update();
        opengl_world.Set_Key_Binding_Category("Object Level Set");
        Add_Component(object_levelset_component,"Object Level Set",'9',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('(',object_levelset_component->viewer_callbacks.Get("toggle_draw_mode"));
        opengl_world.Append_Bind_Key(')',object_levelset_component->viewer_callbacks.Get("toggle_color_mode"));}

    // Particles
    opengl_world.Set_Key_Binding_Category("Particles");
    bool particles_stored_per_cell_uniform=false;
    if(has_valid_grid) particles_stored_per_cell_uniform=true;
    filename=basedir+"/%d/positive_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame) || FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/positive_particles_0",start_frame)){
        positive_particles_component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,basedir+"/%d/positive_particles_%d",true,particles_stored_per_cell_uniform);
        positive_particles_component->particles->template Add_Array<int>(ATTRIBUTE_ID_ID);
        if(!positive_particles_component->Uses_Sets()) positive_particles_component->opengl_points->color=OPENGL_COLOR(1,.5,0);
        Add_Component(positive_particles_component,"Positive particles",'1',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('!',positive_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        opengl_world.Append_Bind_Key('%',positive_particles_component->viewer_callbacks.Get("toggle_draw_radii"));
        opengl_world.Append_Bind_Key('>',positive_particles_component->viewer_callbacks.Get("next_set"));
        opengl_world.Append_Bind_Key('<',positive_particles_component->viewer_callbacks.Get("previous_set"));
        opengl_world.Append_Bind_Key('M',positive_particles_component->viewer_callbacks.Get("toggle_draw_multiple_particle_sets"));}
    filename=basedir+"/%d/negative_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame) || FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/negative_particles_0",start_frame)){
        negative_particles_component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,
            basedir+"/%d/negative_particles_%d",true,particles_stored_per_cell_uniform);
        negative_particles_component->particles->template Add_Array<int>(ATTRIBUTE_ID_ID);
        if(!negative_particles_component->Uses_Sets()) negative_particles_component->opengl_points->color=OPENGL_COLOR(0,.5,1);
        Add_Component(negative_particles_component,"Negative particles",'2',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('@',negative_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        opengl_world.Append_Bind_Key('%',negative_particles_component->viewer_callbacks.Get("toggle_draw_radii"));
        opengl_world.Append_Bind_Key('>',negative_particles_component->viewer_callbacks.Get("next_set"));
        opengl_world.Append_Bind_Key('<',negative_particles_component->viewer_callbacks.Get("previous_set"));
        opengl_world.Append_Bind_Key('M',negative_particles_component->viewer_callbacks.Get("toggle_draw_multiple_particle_sets"));}
    filename=basedir+"/%d/removed_positive_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame) || FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/removed_positive_particles_0",start_frame)){
        removed_positive_particles_component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,basedir+"/%d/removed_positive_particles_%d",
            true,particles_stored_per_cell_uniform);
        removed_positive_particles_component->opengl_points->color=OPENGL_COLOR::Green();
        Add_Component(removed_positive_particles_component,"Removed positive particles",'3',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('#',removed_positive_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        opengl_world.Append_Bind_Key('%',removed_positive_particles_component->viewer_callbacks.Get("toggle_draw_radii"));
        opengl_world.Append_Bind_Key('>',removed_positive_particles_component->viewer_callbacks.Get("next_set"));
        opengl_world.Append_Bind_Key('<',removed_positive_particles_component->viewer_callbacks.Get("previous_set"));
        opengl_world.Append_Bind_Key('M',removed_positive_particles_component->viewer_callbacks.Get("toggle_draw_multiple_particle_sets"));}
    filename=basedir+"/%d/removed_negative_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame) || FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/removed_negative_particles_0",start_frame)){
        removed_negative_particles_component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,
            basedir+"/%d/removed_negative_particles_%d",true,particles_stored_per_cell_uniform);
        removed_negative_particles_component->opengl_points->color=OPENGL_COLOR::Cyan();
        Add_Component(removed_negative_particles_component,"Removed negative particles",'4',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('$',removed_negative_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        opengl_world.Append_Bind_Key('%',removed_negative_particles_component->viewer_callbacks.Get("toggle_draw_radii"));
        opengl_world.Append_Bind_Key('>',removed_negative_particles_component->viewer_callbacks.Get("next_set"));
        opengl_world.Append_Bind_Key('<',removed_negative_particles_component->viewer_callbacks.Get("previous_set"));
        opengl_world.Append_Bind_Key('M',removed_negative_particles_component->viewer_callbacks.Get("toggle_draw_multiple_particle_sets"));}
    filename=basedir+"/%d/particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_2D<T>* particles_component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,"",true,false);
        particles_component->opengl_points->color=OPENGL_COLOR(1,1,1);
        Add_Component(particles_component,"Particles",'P',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);}
    filename=basedir+"/%d/vorticity_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_2D<T>* vorticity_particles_component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,"",false);
        vorticity_particles_component->opengl_points->color=OPENGL_COLOR::Yellow();
        Add_Component(vorticity_particles_component,"Vorticity particles",'O',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/sph_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_2D<T>* sph_particles_component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,"",false);
        sph_particles_component->opengl_points->color=OPENGL_COLOR::Blue();
        Add_Component(sph_particles_component,"SPH particles",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);}

    // Uniform velocities
    OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_component=0;
    OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* ke_component=0;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* vector_velocity_component=0;
    opengl_world.Set_Key_Binding_Category("Velocity");
    filename=basedir+"/%d/velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            vector_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(stream_type,regular_grid,filename);
            vector_velocity_component->opengl_grid_based_vector_field.size=.01;
            vector_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(vector_velocity_component,"Node velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
        else{
            mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
            mac_velocity_component->opengl_mac_velocity_field->size=.01;
            mac_velocity_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Green();
            Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}}
    filename=basedir+"/%d/kinetic_energy";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        ke_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
        ke_component->opengl_mac_velocity_field->size=.01;
        ke_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Green();
        ke_component->opengl_mac_velocity_field->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);
        ke_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(ke_component,"Kinetic Energy",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/mac_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
        mac_velocity_component->opengl_mac_velocity_field->size=.01;
        mac_velocity_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Magenta();
        mac_velocity_component->opengl_mac_velocity_field->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);
        mac_velocity_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/compressible_mac_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
        mac_velocity_component->opengl_mac_velocity_field->size=.01;
        mac_velocity_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Magenta();
        mac_velocity_component->opengl_mac_velocity_field->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);
        mac_velocity_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(mac_velocity_component,"Compressible MAC velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/centered_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        vector_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(stream_type,mac_grid,filename);
        vector_velocity_component->opengl_grid_based_vector_field.size=.01;
        vector_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(vector_velocity_component,"Centered velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/mass_fluxes";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        std::cout<<"Using mass fluxes (from "<<filename<<")"<<std::endl;
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* face_flux_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
        face_flux_component->opengl_mac_velocity_field->size=.01;
        face_flux_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Magenta();
        face_flux_component->opengl_mac_velocity_field->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);
        face_flux_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(face_flux_component,"Mass fluxes",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('M',face_flux_component->viewer_callbacks.Get("toggle_velocity_mode_and_draw"));
        opengl_world.Append_Bind_Key('=',face_flux_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',face_flux_component->viewer_callbacks.Get("decrease_vector_size"));}
    // TODO: integrate these into creation stage...
    if(ke_component){
        opengl_world.Append_Bind_Key('B',ke_component->viewer_callbacks.Get("toggle_velocity_mode_and_draw"));
        opengl_world.Append_Bind_Key('=',ke_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',ke_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',ke_component->viewer_callbacks.Get("toggle_arrowhead"));}
    if(mac_velocity_component){
        opengl_world.Append_Bind_Key('C',mac_velocity_component->viewer_callbacks.Get("toggle_draw"));
        opengl_world.Append_Bind_Key('V',mac_velocity_component->viewer_callbacks.Get("toggle_velocity_mode_and_draw"));
        opengl_world.Append_Bind_Key("<F9>",mac_velocity_component->viewer_callbacks.Get("toggle_draw_divergence"));
        opengl_world.Append_Bind_Key('=',mac_velocity_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',mac_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',mac_velocity_component->viewer_callbacks.Get("toggle_arrowhead"));
        opengl_world.Append_Bind_Key('S',mac_velocity_component->viewer_callbacks.Get("toggle_draw_streamlines"));
        opengl_world.Append_Bind_Key('/',mac_velocity_component->viewer_callbacks.Get("toggle_use_streamline_seed"));
        opengl_world.Append_Bind_Key(']',mac_velocity_component->viewer_callbacks.Get("lengthen_streamlines"));
        opengl_world.Append_Bind_Key('[',mac_velocity_component->viewer_callbacks.Get("shorten_streamlines"));
        opengl_world.Append_Bind_Key('v',mac_velocity_component->viewer_callbacks.Get("toggle_draw_vorticity"));
        opengl_world.Append_Bind_Key('N',mac_velocity_component->viewer_callbacks.Get("normalize_vorticity_color_map"));}
    if(vector_velocity_component){
        opengl_world.Append_Bind_Key('v',vector_velocity_component->viewer_callbacks.Get("toggle_draw"));
        opengl_world.Append_Bind_Key('=',vector_velocity_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',vector_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',vector_velocity_component->viewer_callbacks.Get("toggle_arrowhead"));}

    filename=basedir+"/%d/beta_face";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T>* beta_face_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T>(stream_type,mac_grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.002,OPENGL_COLOR::Gray(1),OPENGL_COLOR::Gray(0)));
        Add_Component(beta_face_component,"Beta Face",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F6),beta_face_component->viewer_callbacks.Get("toggle_draw"));}

    filename=basedir+"/%d/mac_velocities_fuel";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_fuel_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
        mac_velocity_fuel_component->opengl_mac_velocity_field->size=.01;
        mac_velocity_fuel_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Cyan();
        mac_velocity_fuel_component->opengl_mac_velocity_field->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);
        Add_Component(mac_velocity_fuel_component,"MAC velocities fuel",'B',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',mac_velocity_fuel_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',mac_velocity_fuel_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',mac_velocity_fuel_component->viewer_callbacks.Get("toggle_arrowhead"));}

    filename=basedir+"/%d/velocities_ghost_fuel";
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* vector_velocity_ghost_minus_component=0;
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            vector_velocity_ghost_minus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(stream_type,regular_grid,filename);
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.size=.01;
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(vector_velocity_ghost_minus_component,"Node ghost fuel velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key('c',vector_velocity_ghost_minus_component->viewer_callbacks.Get("toggle_draw"));
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_minus_component->viewer_callbacks.Get("toggle_arrowhead"));
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_minus_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_minus_component->viewer_callbacks.Get("decrease_vector_size"));}
        else{
            OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_ghost_minus_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
            mac_velocity_ghost_minus_component->opengl_mac_velocity_field->size=.01;
            mac_velocity_ghost_minus_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Green();
            mac_velocity_ghost_minus_component->Set_Draw(false);
            Add_Component(mac_velocity_ghost_minus_component,"MAC ghost fuel velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key('c',mac_velocity_ghost_minus_component->viewer_callbacks.Get("toggle_velocity_mode_and_draw"));
            opengl_world.Append_Bind_Key('=',mac_velocity_ghost_minus_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',mac_velocity_ghost_minus_component->viewer_callbacks.Get("decrease_vector_size"));
            opengl_world.Append_Bind_Key('h',mac_velocity_ghost_minus_component->viewer_callbacks.Get("toggle_arrowhead"));}}

    filename=basedir+"/%d/velocities_ghost";
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* vector_velocity_ghost_plus_component=0;
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            vector_velocity_ghost_plus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(stream_type,regular_grid,filename);
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.size=.01;
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            vector_velocity_ghost_plus_component->Set_Draw(false);
            Add_Component(vector_velocity_ghost_plus_component,"Node ghost velocities",'b',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_plus_component->viewer_callbacks.Get("toggle_arrowhead"));
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_plus_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_plus_component->viewer_callbacks.Get("decrease_vector_size"));}
        else{
            OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_ghost_plus_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(stream_type,mac_grid,filename);
            mac_velocity_ghost_plus_component->opengl_mac_velocity_field->size=.01;
            mac_velocity_ghost_plus_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Green();
            mac_velocity_ghost_plus_component->Set_Draw(false);
            Add_Component(mac_velocity_ghost_plus_component,"MAC ghost velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key('b',mac_velocity_ghost_plus_component->viewer_callbacks.Get("toggle_velocity_mode_and_draw"));
            opengl_world.Append_Bind_Key('=',mac_velocity_ghost_plus_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',mac_velocity_ghost_plus_component->viewer_callbacks.Get("decrease_vector_size"));
            opengl_world.Append_Bind_Key('h',mac_velocity_ghost_plus_component->viewer_callbacks.Get("toggle_arrowhead"));}}

    if(has_valid_grid && vector_velocity_ghost_minus_component && vector_velocity_ghost_plus_component && levelset_component){
        OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>* two_phase_velocity_magnitude_component=new OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>
            (stream_type,*vector_velocity_ghost_minus_component,*vector_velocity_ghost_plus_component,*levelset_component);
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.plus.size=.01;
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.minus.size=.01;
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.plus.vector_color=OPENGL_COLOR::Magenta();
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.minus.vector_color=OPENGL_COLOR::Cyan();
        Add_Component(two_phase_velocity_magnitude_component,"Two Phase Magnitude",'2',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('\\',two_phase_velocity_magnitude_component->viewer_callbacks.Get("toggle_3d_mode"));
        //opengl_world.Append_Bind_Key('h',two_phase_velocity_magnitude_component->viewer_callbacks.Get("toggle_arrowhead"));
        opengl_world.Append_Bind_Key('=',two_phase_velocity_magnitude_component->viewer_callbacks.Get("increase_point_size"));
        opengl_world.Append_Bind_Key('-',two_phase_velocity_magnitude_component->viewer_callbacks.Get("decrease_point_size"));}

    filename=basedir+"/%d/center_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* center_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(stream_type,mac_grid,filename);
        center_velocity_component->opengl_grid_based_vector_field.size=.01;
        center_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(center_velocity_component,"Centered velocities",'V',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',center_velocity_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',center_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',center_velocity_component->viewer_callbacks.Get("toggle_arrowhead"));}

    opengl_world.Set_Key_Binding_Category("Pressure");
    // TODO: these ramps are leaking memory
    OPENGL_COLOR_MAP<T>* pressure_for_pressure_coupling_color_map=OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(-10000,10000,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1));
    filename=basedir+"/%d/pressure";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        pressure_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,mac_grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1)));
        pressure_component->opengl_scalar_field.Set_Scale_Range(0,30);
        pressure_component->opengl_scalar_field.Set_Uniform_Contour_Values(2,28,(T)26/(T)53);
        Add_Component(pressure_component,"Pressure",'7',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F10),pressure_component->viewer_callbacks.Get("toggle_draw_mode"));}

    opengl_world.Set_Key_Binding_Category("Compressible_Implicit_Pressure");
    filename=basedir+"/%d/compressible_implicit_pressure";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* compressible_implicit_pressure_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,mac_grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1));
        compressible_implicit_pressure_component->opengl_scalar_field.Set_Uniform_Contour_Values(0,20,1);
        Add_Component(compressible_implicit_pressure_component,"Compressible_Implicit_Pressure",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F9),compressible_implicit_pressure_component->viewer_callbacks.Get("toggle_draw"));}
    filename=basedir+"/%d/pressure_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* pressure_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Magenta(1,1)));
        pressure_gradient_component->opengl_scalar_field.Update();
        Add_Component(pressure_gradient_component,"Pressure Gradient",'9',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename=basedir+"/%d/pressure_for_pressure_coupling";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* pressure2_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(stream_type,mac_grid,filename,pressure_for_pressure_coupling_color_map);
        pressure2_component->opengl_scalar_field.Set_Uniform_Contour_Values(-10000,10000,100);
        Add_Component(pressure2_component,"Pressure2",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F7),pressure2_component->viewer_callbacks.Get("toggle_draw"));}

    // Draw grid here so it'll be above particles and pressure
    opengl_world.Set_Key_Binding_Category("Grid");
    if(has_valid_grid){
        OPENGL_GRID_2D<T>* opengl_grid=new OPENGL_GRID_2D<T>(stream_type,*(new GRID<TV>(grid)),OPENGL_COLOR::Gray(.5));
        grid_component=new OPENGL_COMPONENT_BASIC<T,OPENGL_GRID_2D<T> >(stream_type,*opengl_grid);
        Add_Component(grid_component,"Grid",'6',BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('^',grid_component->object.viewer_callbacks.Get("toggle_draw_ghost_values"));}
    opengl_world.Set_Key_Binding_Category("Sub Grids");

    // deformable and rigid bodies
    OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>* rigid_bodies_component=0;
    OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>* deformable_objects_component=0;
    std::string deformable_object_filename=LOG::sprintf("%s/%d/deformable_object_particles",basedir.c_str(),start_frame);
    if(FILE_UTILITIES::File_Exists(deformable_object_filename)){
        deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>(stream_type,basedir+"/",start_frame);
        deformable_objects_component->Set_Vector_Size(.01);
        deformable_objects_component->selectable=true;
        // TODO: what the hell?
        if(!FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/rigid_body_particles",start_frame)){
            }}
    if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/rigid_body_particles",start_frame)){
        rigid_bodies_component=new OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T>(stream_type,basedir);
        rigid_bodies_component->Set_Vector_Size(.01);
        rigid_bodies_component->selectable=true;
        if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/colors",start_frame)) FILE_UTILITIES::template Read_From_File(stream_type,LOG::sprintf("%s/%d/colors",basedir.c_str(),start_frame),rigid_bodies_component->colors);
        for(int i=0;i<rigid_bodies_no_draw_list.m;i++){
            LOG::cout<<"Rigid bodies: not drawing object "<<rigid_bodies_no_draw_list(i)<<std::endl;
            rigid_bodies_component->Set_Draw_Object(rigid_bodies_no_draw_list(i),false);}
        opengl_world.Set_Key_Binding_Category("Rigid Bodies");
        Add_Component(rigid_bodies_component,"Rigid Bodies",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('5',rigid_bodies_component->viewer_callbacks.Get("toggle_draw_mode"));
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->viewer_callbacks.Get("toggle_velocity_vectors"));
        opengl_world.Append_Bind_Key('a',rigid_bodies_component->viewer_callbacks.Get("toggle_articulation_points"));
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->viewer_callbacks.Get("toggle_show_object_names"));
        opengl_world.Append_Bind_Key('=',rigid_bodies_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',rigid_bodies_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('m',rigid_bodies_component->viewer_callbacks.Get("toggle_linear_muscles"));}
    if(deformable_objects_component){
        opengl_world.Set_Key_Binding_Category("Deformable Objects");
        Add_Component(deformable_objects_component,"Deformable Objects",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('8',deformable_objects_component->viewer_callbacks.Get("toggle_draw"));
        opengl_world.Append_Bind_Key('e',deformable_objects_component->viewer_callbacks.Get("cycle_display_mode"));
        opengl_world.Append_Bind_Key('*',deformable_objects_component->viewer_callbacks.Get("toggle_draw_velocities"));
        opengl_world.Append_Bind_Key('=',deformable_objects_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',deformable_objects_component->viewer_callbacks.Get("decrease_vector_size"));}

    std::string soft_constraints_deformable_object_filename=basedir+LOG::sprintf("/%d/soft_constraints_deformable_object_particles",start_frame);
    if(FILE_UTILITIES::File_Exists(soft_constraints_deformable_object_filename)){
        OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>* soft_constraints_deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T>(stream_type,basedir+"/soft_constraints_",start_frame);
        soft_constraints_deformable_objects_component->Set_Vector_Size(.01);
        Add_Component(soft_constraints_deformable_objects_component,"Soft Constraints Deformable Objects",'e',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);}

    opengl_world.Set_Key_Binding_Category("Fluid Boundaries");
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N",start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,bool>* psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,bool>(stream_type,grid,basedir+"/%d/psi_N",
            new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()));
        Add_Component(psi_N_component,"Psi_N points",'\0',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_N_component->viewer_callbacks.Get("toggle_draw"));}
    filename=basedir+"/%d/debug_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>* component=new OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>(stream_type,filename);
        Add_Component(component,"Debug particles",'w',BASIC_VISUALIZATION<T>::SELECTABLE|BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key('W',component->viewer_callbacks.Get("toggle_draw_velocities"));
        opengl_world.Append_Bind_Key('h',component->viewer_callbacks.Get("toggle_arrowhead"));
        opengl_world.Append_Bind_Key('=',component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',component->viewer_callbacks.Get("decrease_vector_size"));}
    filename=basedir+"/%d/residual_energy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_2D<T>* component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,"",false,false);
        Add_Component(component,"Residual Energy",'k',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED);}
    filename=basedir+"/%d/collision_iterators";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_2D<T>* component=new OPENGL_COMPONENT_PARTICLES_2D<T>(stream_type,filename,"",false,false);
        Add_Component(component,"Collision Iterators",'I',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED);}
    filename=basedir+"/%d/psi_D";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>* psi_D_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>(stream_type,mac_grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS);
        Add_Component(psi_D_component,"Psi_D points",'\0',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_D_component->viewer_callbacks.Get("toggle_draw"));}

    filename=basedir+"/%d/euler_psi";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>* psi_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>(stream_type,mac_grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Red()),OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS);
        Add_Component(psi_component,"Psi points",'\0',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F3),psi_component->viewer_callbacks.Get("toggle_draw"));}

    filename=basedir+"/%d/colors";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_INDEXED_COLOR_MAP* colors_color_map=OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map();colors_color_map->Set_Index_Mode(OPENGL_INDEXED_COLOR_MAP::PERIODIC);
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T,int>* psi_colors_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T,int>(stream_type,mac_grid,filename,colors_color_map,OPENGL_SCALAR_FIELD_2D<T,int>::DRAW_POINTS);
        Add_Component(psi_colors_component,"Psi colors",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_colors_component->viewer_callbacks.Get("toggle_draw"));}

    filename=basedir+"/%d/pressure_jumps";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* pressure_jump_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(stream_type,grid,filename);
        pressure_jump_component->opengl_grid_based_vector_field.size=.01;
        pressure_jump_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Magenta();
        Add_Component(pressure_jump_component,"Pressure jumps",'&',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',pressure_jump_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',pressure_jump_component->viewer_callbacks.Get("decrease_vector_size"));}

    filename=basedir+"/%d/thin_shells_grid_visibility";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>* thin_shells_debugging_component=new OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>(stream_type,grid,basedir);
        opengl_world.Set_Key_Binding_Category("Thin Shells");
        Add_Component(thin_shells_debugging_component,"thin shells debugging",'\0',BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F4),thin_shells_debugging_component->viewer_callbacks.Get("toggle_draw_grid_visibility"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),thin_shells_debugging_component->viewer_callbacks.Get("toggle_draw_density_valid_mask"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),thin_shells_debugging_component->viewer_callbacks.Get("toggle_draw_phi_valid_mask"));}

    filename=basedir+"/%d/strain";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>* strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>(stream_type,grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->viewer_callbacks.Get("increase_size"));
        opengl_world.Append_Bind_Key('_',strain_component->viewer_callbacks.Get("decrease_size"));}

    filename=basedir+"/%d/strain_0";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>* strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>(stream_type,grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->viewer_callbacks.Get("increase_size"));
        opengl_world.Append_Bind_Key('_',strain_component->viewer_callbacks.Get("decrease_size"));}

    // Commands and misc
    opengl_world.Set_Key_Binding_Category("Misc.");
    opengl_world.Append_Bind_Key('~',{[this](){Command_Prompt();},"command_prompt"});
    opengl_world.Append_Bind_Key('\\',{[this](){Toggle_2d_Mode();},"toggle_2d_mode"});

    // Initialize priority ordering for selections
    Selection_Priority(OPENGL_SELECTION<T>::DEBUG_PARTICLES_2D)=110;
    Selection_Priority(OPENGL_SELECTION<T>::POINTS_2D)=100;
    Selection_Priority(OPENGL_SELECTION<T>::COMPONENT_PARTICLES_2D)=100;
    Selection_Priority(OPENGL_SELECTION<T>::ARTICULATED_RIGID_BODIES_JOINT_2D)=95;
    Selection_Priority(OPENGL_SELECTION<T>::ARTICULATED_RIGID_BODIES_MUSCLE_2D)=95;
    Selection_Priority(OPENGL_SELECTION<T>::COMPONENT_RIGID_BODIES_2D)=80;
    Selection_Priority(OPENGL_SELECTION<T>::COMPONENT_DEFORMABLE_OBJECT_2D)=80;
    Selection_Priority(OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_2D)=79;
    Selection_Priority(OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_2D)=78;
    Selection_Priority(OPENGL_SELECTION<T>::TRIANGULATED_AREA_VERTEX)=77;
    Selection_Priority(OPENGL_SELECTION<T>::TRIANGULATED_AREA_SEGMENT)=76;
    Selection_Priority(OPENGL_SELECTION<T>::TRIANGULATED_AREA_TRIANGLE)=75;
    Selection_Priority(OPENGL_SELECTION<T>::BEZIER_SPLINE_VERTEX_2D)=74;
    Selection_Priority(OPENGL_SELECTION<T>::BEZIER_SPLINE_SEGMENT_2D)=73;
    Selection_Priority(OPENGL_SELECTION<T>::B_SPLINE_VERTEX_2D)=72;
    Selection_Priority(OPENGL_SELECTION<T>::B_SPLINE_SEGMENT_2D)=71;
    Selection_Priority(OPENGL_SELECTION<T>::GRID_NODE_2D)=70;
    Selection_Priority(OPENGL_SELECTION<T>::GRID_CELL_2D)=60;
}
//#####################################################################
// Add_OpenGL_Initialization
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION<T>::Add_OpenGL_Initialization();
    opengl_world.Set_2D_Mode(true);
}
//#####################################################################
// Pre_Frame_Extra
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Pre_Frame_Extra()
{
    if(grid_component) grid_component->object.Set_Frame(frame);
}
//#####################################################################
// Set_Frame_Extra
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Set_Frame_Extra()
{
    std::string filename=LOG::sprintf("%s/%d/frame_title",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){std::ifstream input(filename.c_str());getline(input,frame_title);}
    else frame_title="";
    filename=LOG::sprintf("%s/%d/time",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){T time;FILE_UTILITIES::Read_From_File(stream_type,filename,time);frame_title=LOG::sprintf("(%.05f) ",time)+frame_title;}
}
//#####################################################################
// Command_Prompt_Response
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Command_Prompt_Response()
{
    if(!opengl_world.prompt_response.empty()){
        std::string command;
        std::istringstream sstream(opengl_world.prompt_response);
        sstream>>command;
        if(command=="s"){
            int id;
            if(sstream>>id){
                OPENGL_SELECTION<T>* selection=0;
                // TODO : handle multiple particle sets
                if(!selection && positive_particles_component) selection=positive_particles_component->Get_Selection_By_Id(id,1);
                if(!selection && negative_particles_component) selection=negative_particles_component->Get_Selection_By_Id(id,1);
                if(!selection && removed_positive_particles_component) selection=removed_positive_particles_component->Get_Selection_By_Id(id,1);
                if(!selection && removed_negative_particles_component) selection=removed_negative_particles_component->Get_Selection_By_Id(id,1);
                if(selection){Set_Current_Selection(selection);Update_OpenGL_Strings();}}}
        else if(command=="p"){
            std::string next_word;
            if(sstream>>next_word){
                if(next_word=="auto"){
                    if(pressure_component){
                        T min_pressure,max_pressure;
                        min_pressure=pressure_component->opengl_scalar_field.values.Min();
                        max_pressure=pressure_component->opengl_scalar_field.values.Max();
                        LOG::cout<<"min pressure="<<min_pressure<<", max pressure="<<max_pressure<<std::endl;
                        pressure_component->opengl_scalar_field.Set_Scale_Range(min_pressure,max_pressure);
                        pressure_component->opengl_scalar_field.Update();}}
                else{
                    T min_pressure=atof(next_word.c_str()),max_pressure;
                    if(sstream>>max_pressure){
                        if(pressure_component){
                            LOG::cout<<"min pressure="<<min_pressure<<", max pressure="<<max_pressure<<std::endl;
                            pressure_component->opengl_scalar_field.Set_Scale_Range(min_pressure,max_pressure);
                            pressure_component->opengl_scalar_field.Update();}}}}}
        else if(command=="skip"){
            int skip;
            sstream>>skip;
            frame_increment=skip;}
        else if(command=="sphcw"){
            if(sph_cell_weight_component){
                std::ostringstream output_stream;int x,y;
                if(sstream>>x && sstream>>y){
                    VECTOR<int,2> index(x,y);
                    if(sph_cell_weight_component->Is_Up_To_Date(frame))
                        if(sph_cell_weight_component->opengl_scalar_field.values.Valid_Index(index))
                            output_stream<<"SPH cell weight at "<<x<<" "<<y<<": "<<sph_cell_weight_component->opengl_scalar_field.values(index)<<std::endl;
                        else output_stream<<x<<" "<<y<<" is not a valid cell index"<<std::endl;
                    else output_stream<<"frame "<<frame<<" is not up to date"<<std::endl;}
                else output_stream<<"invalid input"<<std::endl;
                opengl_world.Add_String(output_stream.str());}}}
}
//#####################################################################
// Toggle_2d_Mode
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Toggle_2d_Mode()
{
    opengl_world.Set_2D_Mode(!opengl_world.mode_2d);
}
//#####################################################################
// Command_Prompt
//#####################################################################
template<class T> void OPENGL_2D_VISUALIZATION<T>::
Command_Prompt()
{
    opengl_world.Prompt_User("Command: ",{[this](){Command_Prompt_Response();},"command_prompt_response"});
}
//#####################################################################
// main
//#####################################################################
int main(int argc,char* argv[])
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    bool type_double=false; // float by default
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Add_Not("-float",&type_double,"Use floats");
    parse_args.Add("-double",&type_double,"Use doubles");
    parse_args.Parse(true);

    if(!type_double)
    {
        ANIMATED_VISUALIZATION<float> *visualization=new OPENGL_2D_VISUALIZATION<float>(STREAM_TYPE(1.f));
        visualization->Initialize_And_Run(parse_args);
        delete visualization;
    }
    else
    {
        ANIMATED_VISUALIZATION<double> *visualization=new OPENGL_2D_VISUALIZATION<double>(STREAM_TYPE(1.));
        visualization->Initialize_And_Run(parse_args);
        delete visualization;
    }

    return 0;
}
