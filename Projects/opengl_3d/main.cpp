//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Kevin Der, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Cynthia Lau, Michael Lentine, Sergey Levine, Frank Losasso, Nick Rasmussen, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Jerry Talton, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Incompressible/Particles/VORTICITY_PARTICLES.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <Dynamics/Particles/SPH_PARTICLES.h>
#include <OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_LEVELSET_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <OpenGL/OpenGL/OPENGL_SLICE_MANAGER.h>
#include <OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEBUG_PARTICLES_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DIAGNOSTICS.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_VORTICITY_PARTICLES_3D.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace PhysBAM;
using namespace std;

template<class T>
class VISUALIZATION:public ANIMATED_VISUALIZATION<T>
{
    typedef VECTOR<T,3> TV;
public:
    using ANIMATED_VISUALIZATION<T>::last_frame_filename;using ANIMATED_VISUALIZATION<T>::camera_script_filename;
    using ANIMATED_VISUALIZATION<T>::opengl_window_title;using ANIMATED_VISUALIZATION<T>::opengl_world;
    using ANIMATED_VISUALIZATION<T>::Add_Component;using ANIMATED_VISUALIZATION<T>::Selection_Priority;
    using ANIMATED_VISUALIZATION<T>::frame_rate;using ANIMATED_VISUALIZATION<T>::start_frame;
    using ANIMATED_VISUALIZATION<T>::frame;using ANIMATED_VISUALIZATION<T>::draw_all_objects_cb;
    using ANIMATED_VISUALIZATION<T>::frame_title;using ANIMATED_VISUALIZATION<T>::Set_Current_Selection;
    using ANIMATED_VISUALIZATION<T>::stream_type;
    VISUALIZATION(STREAM_TYPE stream_type);
    virtual ~VISUALIZATION();

protected:
    virtual void Add_Arguments(PARSE_ARGS &parse_args);
    virtual void Parse_Arguments(PARSE_ARGS &parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Update_OpenGL_Strings();

private:
    void Read_Grid();
    virtual void Pre_Frame_Extra();
    virtual void Set_Frame_Extra();

    // callbacks
    void Command_Prompt();
    void Command_Prompt_Response();
    void Slice_Has_Changed();

    // Components
    OPENGL_COMPONENT_PARTICLES_3D<T>* positive_particles_component;
    OPENGL_COMPONENT_PARTICLES_3D<T>* negative_particles_component;
    OPENGL_COMPONENT_PARTICLES_3D<T>* removed_positive_particles_component;
    OPENGL_COMPONENT_PARTICLES_3D<T>* removed_negative_particles_component;

    // TODO: need better grid control 
    OPENGL_COMPONENT_BASIC<T,OPENGL_GRID_3D<T> >* grid_component;

    // Options
    std::string basedir;
    bool has_diagnostics;

    GRID<TV> grid,mac_grid,regular_grid;
    bool has_valid_grid;
    bool node_based;
    OPENGL_SLICE_MANAGER<T> slice_manager;
    OPENGL_BOX_3D<T>* opengl_box;
    OPENGL_UNIFORM_SLICE<T>* slice;

    ARRAY<int> rigid_bodies_no_draw_list;
    ARRAY<int> deformable_no_draw_list;

    bool allow_caching;
    bool always_add_mac_velocities;
};

// ------------------------------------------------------------------

template<class T> VISUALIZATION<T>::
VISUALIZATION(STREAM_TYPE stream_type)
    :ANIMATED_VISUALIZATION<T>(stream_type),positive_particles_component(0),negative_particles_component(0),
    removed_positive_particles_component(0),removed_negative_particles_component(0),grid_component(0),opengl_box(0),slice(0),
    allow_caching(true),always_add_mac_velocities(false)
{
}

template<class T> VISUALIZATION<T>::
~VISUALIZATION()
{
    if(opengl_box) delete &opengl_box->box;
    delete slice;
}

template<class T> void VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    basedir=".";

    ANIMATED_VISUALIZATION<T>::Add_Arguments(parse_args);

    parse_args.Add_Not("-no_caching",&allow_caching,"Allow caching");
    parse_args.Add("-macvel",&always_add_mac_velocities,"force adding mac velocities component");
    parse_args.Add("-rigid_bodies_no_draw",&rigid_bodies_no_draw_list,"id","Do not draw this rigid body (may be repeated)");
    parse_args.Add("-deformable_no_draw",&deformable_no_draw_list,"id","Do not draw this deformable body (may be repeated)");
    parse_args.Extra_Optional(&basedir,"basedir","base directory");
}

template<class T> void VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
#ifdef __linux__
    opengl_window_title="opengl_3d: " + FILE_UTILITIES::Real_Path(basedir);
#endif

    if(FILE_UTILITIES::File_Exists(basedir+"/common/first_frame")){
       int start_frame_from_file;
       FILE_UTILITIES::Read_From_Text_File(basedir+"/common/first_frame",start_frame_from_file);
       start_frame=std::max(start_frame,start_frame_from_file);}

    ANIMATED_VISUALIZATION<T>::Parse_Arguments(parse_args);

    if(parse_args.unclaimed_arguments)
        parse_args.Print_Usage(true);

    last_frame_filename=basedir+"/common/last_frame";

    // don't override camera script filename if it was already set in base class based on command line argument
    if(camera_script_filename.empty()) camera_script_filename=basedir+"/camera_script";
}

template<class T> void VISUALIZATION<T>::
Read_Grid()
{
    has_valid_grid=false;
    std::string filename,coarse_filename;

    filename=LOG::sprintf("%s/%d/levelset",basedir.c_str(),start_frame);
    coarse_filename=LOG::sprintf("%s/%d/coarse_levelset",basedir.c_str(),start_frame);
    // For backwards compatibility
    if(!FILE_UTILITIES::File_Exists(filename)) filename=LOG::sprintf("%s/%d/levelset.phi",basedir.c_str(),start_frame);

    if(FILE_UTILITIES::File_Exists(filename)){
        std::cout<<"Reading grid from '"<<filename<<"'..."<<std::endl;
        ARRAY<T,VECTOR<int,3> > phi;
        LEVELSET<TV> levelset(grid,phi);
        FILE_UTILITIES::Read_From_File(stream_type,filename,levelset);
        has_valid_grid=true;}
    else if(FILE_UTILITIES::File_Exists(basedir+"/common/grid")){
        filename=basedir+"/common/grid";
        std::cout<<"Reading grid from '"<<filename<<"'..."<<std::endl;
        FILE_UTILITIES::Read_From_File(stream_type,filename,grid);
        has_valid_grid=true;}

    if(has_valid_grid){
        node_based=!grid.Is_MAC_Grid();
        mac_grid=grid.Get_MAC_Grid();regular_grid=grid.Get_Regular_Grid();}
}

//#####################################################################
// Initialize_Components_And_Key_Bindings
//#####################################################################
template<class T> void VISUALIZATION<T>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION<T>::Initialize_Components_And_Key_Bindings();
    opengl_world.Set_Key_Binding_Category_Priority(1);
    opengl_world.Unbind_Keys("abBCdDEeFjJjkKlLMotTvV 1!2@3#4$5%67&89 ^=-`{}\b\\[]~\t");

    std::string filename,filename2,coarse_filename;
    filename=basedir+"/common/sim.param";
    PARAMETER_LIST parameter_list;
    parameter_list.Read(filename);

    frame_rate=parameter_list.template Get_Parameter<int>("frame_rate",frame_rate);
    std::cout<<"Using frame rate "<<frame_rate<<std::endl;


    filename=basedir+"/%d/diagnostics";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_DIAGNOSTICS<T>(stream_type,filename),"Diagnostics",'\0',BASIC_VISUALIZATION<T>::OWNED);

    Read_Grid();
    if(has_valid_grid){
        slice=new OPENGL_UNIFORM_SLICE<T>(opengl_world);
        slice->Initialize(grid);
        slice_manager.slice=slice;
        std::cout<<"Using uniform grid slice"<<std::endl;}

    if(slice_manager.slice)
        slice_manager.Set_Slice_Has_Changed_Callback({[this](){Slice_Has_Changed();},"slice_has_changed"});

    if(has_valid_grid){
        opengl_box=new OPENGL_BOX_3D<T>(stream_type,*(new RANGE<TV>(grid.Domain())),OPENGL_COLOR::Gray(0.5));
        OPENGL_COMPONENT_BASIC<T,OPENGL_BOX_3D<T> >* domain_box_component=new OPENGL_COMPONENT_BASIC<T,OPENGL_BOX_3D<T> >(stream_type,*opengl_box);
        Add_Component(domain_box_component,"Domain box",'6',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);}

    {OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>* deformable_objects_component=0;
    std::string deformable_object_filename=basedir+LOG::sprintf("/%d/deformable_object_structures",start_frame);
    if(FILE_UTILITIES::File_Exists(basedir+"/common/deformable_object_structures") || FILE_UTILITIES::File_Exists(deformable_object_filename)){
        OPENGL_MATERIAL front_material=OPENGL_MATERIAL::Matte(OPENGL_COLOR::Yellow());
        OPENGL_MATERIAL back_material=OPENGL_MATERIAL::Matte(OPENGL_COLOR::Yellow(0.5));
        deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>(stream_type,basedir+"/",start_frame);
        deformable_objects_component->Set_All_Materials(front_material,front_material,back_material);
        if(deformable_no_draw_list.m){
            deformable_objects_component->active_list.Resize(PhysBAM::max(deformable_objects_component->active_list.m,deformable_no_draw_list.Last()),true,true,true);
            deformable_objects_component->active_list.Subset(deformable_no_draw_list).Fill(false);}
        opengl_world.Set_Key_Binding_Category("Deformable Objects");
        Add_Component(deformable_objects_component,"Deformable Objects",'8',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('9',deformable_objects_component->viewer_callbacks.Get("toggle_active_value"));
        opengl_world.Append_Bind_Key('(',deformable_objects_component->viewer_callbacks.Get("show_only_first"));
        opengl_world.Append_Bind_Key('W',deformable_objects_component->viewer_callbacks.Get("toggle_selection_mode"));
        opengl_world.Append_Bind_Key('h',deformable_objects_component->viewer_callbacks.Get("toggle_hide_unselected"));
        opengl_world.Append_Bind_Key('b',deformable_objects_component->viewer_callbacks.Get("toggle_draw_interior"));
        opengl_world.Append_Bind_Key('I',deformable_objects_component->viewer_callbacks.Get("toggle_differentiate_inverted"));
        opengl_world.Append_Bind_Key('t',deformable_objects_component->viewer_callbacks.Get("toggle_draw_subsets"));
        opengl_world.Append_Bind_Key('e',deformable_objects_component->viewer_callbacks.Get("cycle_display_mode"));
        opengl_world.Append_Bind_Key('i',deformable_objects_component->viewer_callbacks.Get("cycle_interaction_pair_display_mode"));
        opengl_world.Append_Bind_Key('$',deformable_objects_component->viewer_callbacks.Get("toggle_velocity_vectors"));
        opengl_world.Append_Bind_Key('=',deformable_objects_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',deformable_objects_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F8),deformable_objects_component->viewer_callbacks.Get("cycle_forces_mode"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F7),deformable_objects_component->viewer_callbacks.Get("cycle_relative_velocity_mode"));
        if(deformable_objects_component->has_embedded_objects || deformable_objects_component->has_soft_bindings)
            opengl_world.Append_Bind_Key('h',deformable_objects_component->viewer_callbacks.Get("cycle_hard_bound_surface_display_mode"));
        if(deformable_objects_component->has_tetrahedralized_volumes){
            opengl_world.Append_Bind_Key('C',deformable_objects_component->viewer_callbacks.Get("cycle_cutaway_mode"));
            opengl_world.Append_Bind_Key('{',deformable_objects_component->viewer_callbacks.Get("decrease_cutaway_fraction"));
            opengl_world.Append_Bind_Key('}',deformable_objects_component->viewer_callbacks.Get("increase_cutaway_fraction"));}
        opengl_world.Append_Bind_Key('Z',deformable_objects_component->viewer_callbacks.Get("highlight_particle"));
        if(slice_manager.slice) slice_manager.Add_Object(deformable_objects_component);}

    if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/rigid_body_particles",start_frame)){
        OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>* rigid_bodies_component=0;
        rigid_bodies_component=new OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>(stream_type,basedir,true);
        rigid_bodies_component->Set_Vector_Size(0.01);
        rigid_bodies_component->Set_Frame(start_frame); // needed before reinitialize so that no draw list will work
        rigid_bodies_component->Reinitialize();
        for(int i=0;i<rigid_bodies_no_draw_list.m;i++){
            std::cout<<"Rigid bodies: not drawing object "<<rigid_bodies_no_draw_list(i)<<std::endl;
            rigid_bodies_component->Set_Draw_Object(rigid_bodies_no_draw_list(i),false);}
        rigid_bodies_component->Reinitialize();
        opengl_world.Set_Key_Binding_Category("Rigid Bodies");
        Add_Component(rigid_bodies_component,"Rigid Bodies",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('5',rigid_bodies_component->viewer_callbacks.Get("toggle_draw_mode"));
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->viewer_callbacks.Get("toggle_velocity_vectors"));
        opengl_world.Append_Bind_Key('a',rigid_bodies_component->viewer_callbacks.Get("toggle_individual_axes"));
        opengl_world.Append_Bind_Key('A',rigid_bodies_component->viewer_callbacks.Get("toggle_articulation_points"));
        opengl_world.Append_Bind_Key('F',rigid_bodies_component->viewer_callbacks.Get("toggle_joint_frames"));
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->viewer_callbacks.Get("toggle_show_object_names"));
        opengl_world.Append_Bind_Key('=',rigid_bodies_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',rigid_bodies_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('M',rigid_bodies_component->viewer_callbacks.Get("toggle_draw_particles"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),rigid_bodies_component->viewer_callbacks.Get("toggle_forces_and_torques"));
        opengl_world.Append_Bind_Key('o',rigid_bodies_component->viewer_callbacks.Get("toggle_one_sided"));
        if(slice_manager.slice) slice_manager.Add_Object(rigid_bodies_component);}

    std::string soft_constraints_deformable_object_filename=basedir+"/soft_constraints_deformable_object_particles";
    if(FILE_UTILITIES::File_Exists(soft_constraints_deformable_object_filename) // TODO(jontg): Not sure what to do here...
        || FILE_UTILITIES::File_Exists(soft_constraints_deformable_object_filename+LOG::sprintf(".%d",start_frame))){
        OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>* soft_constraints_deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T>(stream_type,basedir+"/soft_constraints_",start_frame);
        soft_constraints_deformable_objects_component->selectable=true;
        Add_Component(soft_constraints_deformable_objects_component,"Soft Constraints Deformable Objects",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        if(slice_manager.slice) slice_manager.Add_Object(soft_constraints_deformable_objects_component);}

    filename=basedir+"/%d/levelset";
    coarse_filename=basedir+"/%d/coarse_levelset";
    filename2=allow_caching?(basedir+"/%d/levelset_tesselated"):"";
    // for backwards compatiblity
    if(!FILE_UTILITIES::Frame_File_Exists(filename,start_frame) && !FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_0",start_frame)){
        filename=basedir+"/levelset_%d.phi";
        filename2=allow_caching?(basedir+"/levelset_%d.tri"):"";}
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)||FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_0",start_frame)){
        OPENGL_COMPONENT_LEVELSET_3D<T>* levelset_component=new OPENGL_COMPONENT_LEVELSET_3D<T>(stream_type,filename,filename2,basedir+"/levelset_%d.%d",basedir+"/levelset_tesselated_%d.%d",true);
        opengl_world.Set_Key_Binding_Category("Levelset");
        Add_Component(levelset_component,"Levelset",'l',BASIC_VISUALIZATION<T>::OWNED);
//        levelset_component->Set_Surface_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T).65,1)),OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T).65,1)));
        if(slice_manager.slice) slice_manager.Add_Object(levelset_component);
        opengl_world.Append_Bind_Key('L',levelset_component->viewer_callbacks.Get("toggle_slice_color_mode"));
        opengl_world.Append_Bind_Key("^l",levelset_component->viewer_callbacks.Get("toggle_display_overlay"));
        opengl_world.Append_Bind_Key('`',levelset_component->viewer_callbacks.Get("toggle_smooth_slice"));
        if(levelset_component->Use_Sets()){
            opengl_world.Append_Bind_Key('M',levelset_component->viewer_callbacks.Get("toggle_draw_multiple_levelsets"));
            opengl_world.Append_Bind_Key('>',levelset_component->viewer_callbacks.Get("next_set"));
            opengl_world.Append_Bind_Key('<',levelset_component->viewer_callbacks.Get("previous_set"));}}
    if(FILE_UTILITIES::Frame_File_Exists(coarse_filename,start_frame)){
        OPENGL_COMPONENT_LEVELSET_3D<T>* levelset_component=new OPENGL_COMPONENT_LEVELSET_3D<T>(stream_type,coarse_filename,filename2,basedir+"/coarse_levelset_%d.%d",basedir+"/coarse_levelset_tesselated_%d.%d",true);
        levelset_component->ghost_cells=0;
        opengl_world.Set_Key_Binding_Category("Levelset");
        Add_Component(levelset_component,"Levelset coarse",'.',BASIC_VISUALIZATION<T>::OWNED);
        levelset_component->Set_Surface_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T)1,(T).65)),OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T)1,(T).65)));
        if(slice_manager.slice) slice_manager.Add_Object(levelset_component);
        opengl_world.Append_Bind_Key('>',levelset_component->viewer_callbacks.Get("toggle_slice_color_mode"));
        opengl_world.Append_Bind_Key("^.",levelset_component->viewer_callbacks.Get("toggle_display_overlay"));}

    filename=basedir+"/%d/object_levelset";
    filename2=basedir+"/%d/object_levelset_tesselated";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_LEVELSET_3D<T>* object_levelset_component=new OPENGL_COMPONENT_LEVELSET_3D<T>(stream_type,filename,filename2,"","",true);
        Add_Component(object_levelset_component,"Object Levelset",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        if(slice_manager.slice) slice_manager.Add_Object(object_levelset_component);
        opengl_world.Set_Key_Binding_Category("Object Levelset");
        opengl_world.Append_Bind_Key("^l",object_levelset_component->viewer_callbacks.Get("toggle_draw"));
        opengl_world.Append_Bind_Key('L',object_levelset_component->viewer_callbacks.Get("toggle_slice_color_mode"));
        opengl_world.Append_Bind_Key('`',object_levelset_component->viewer_callbacks.Get("toggle_smooth_slice"));}

    opengl_world.Set_Key_Binding_Category("Density");
    filename=basedir+"/%d/density";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* density_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_component->opengl_scalar_field.Update();
        Add_Component(density_component,"Density",'d',BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key('D',density_component->viewer_callbacks.Get("toggle_color_map"));
        opengl_world.Append_Bind_Key('`',density_component->viewer_callbacks.Get("toggle_smooth_slice"));
        slice_manager.Add_Object(density_component);}

    filename=basedir+"/%d/density_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* density_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_gradient_component->opengl_scalar_field.Update();
        Add_Component(density_gradient_component,"Density Gradient",'8',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        slice_manager.Add_Object(density_gradient_component);}
        
    opengl_world.Set_Key_Binding_Category("Soot");
    filename=basedir+"/%d/soot";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* soot_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.01,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        soot_component->opengl_scalar_field.Update();
        Add_Component(soot_component,"Soot",'i',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',soot_component->viewer_callbacks.Get("toggle_color_map"));
        opengl_world.Append_Bind_Key('`',soot_component->viewer_callbacks.Get("toggle_smooth_slice"));
        slice_manager.Add_Object(soot_component);}

    filename=basedir+"/%d/soot_fuel";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* soot_fuel_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.01,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        soot_fuel_component->opengl_scalar_field.Update();
        Add_Component(soot_fuel_component,"soot_fuel",'f',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',soot_fuel_component->viewer_callbacks.Get("toggle_color_map"));
        opengl_world.Append_Bind_Key('`',soot_fuel_component->viewer_callbacks.Get("toggle_smooth_slice"));
        slice_manager.Add_Object(soot_fuel_component);}

    opengl_world.Set_Key_Binding_Category("Internal Energy");
    filename=basedir+"/%d/internal_energy";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* internal_energy_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        internal_energy_component->opengl_scalar_field.Update();
        Add_Component(internal_energy_component,"Internal Energy",'E',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        slice_manager.Add_Object(internal_energy_component);}

    filename=basedir+"/%d/debug_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>* component=new OPENGL_COMPONENT_DEBUG_PARTICLES_3D<T>(stream_type,filename);
        Add_Component(component,"Debug particles",'w',BASIC_VISUALIZATION<T>::SELECTABLE|BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key('W',component->viewer_callbacks.Get("toggle_draw_velocities"));
        opengl_world.Append_Bind_Key('q',component->viewer_callbacks.Get("show_colored_wireframe"));
        opengl_world.Append_Bind_Key('=',component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',component->viewer_callbacks.Get("decrease_vector_size"));
        if(slice_manager.slice) slice_manager.Add_Object(component);}

    opengl_world.Set_Key_Binding_Category("Temperature");
    filename=basedir+"/%d/temperature";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        ARRAY<T,VECTOR<int,3> > temp0;FILE_UTILITIES::Read_From_File(stream_type,FILE_UTILITIES::Get_Frame_Filename(filename,start_frame),temp0);
        T min_temp=temp0.Min(),max_temp=temp0.Max();min_temp=283.15;
        std::cout<<"Using temperature: [min="<<min_temp<<" max="<<max_temp<<"]"<<std::endl;
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* temperature_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Red(1,1)));
        temperature_component->opengl_scalar_field.Set_Scale_Range(temp0.Min(),temp0.Max());
        Add_Component(temperature_component,"Temperature",'t',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('T',temperature_component->viewer_callbacks.Get("toggle_color_map"));
        slice_manager.Add_Object(temperature_component);}

    opengl_world.Set_Key_Binding_Category("SPH");
    filename=basedir+"/%d/sph_cell_weights";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* sph_cell_weight_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,100,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1)));
        sph_cell_weight_component->opengl_scalar_field.Update();
        Add_Component(sph_cell_weight_component,"SPH Cell Weights",'H',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        slice_manager.Add_Object(sph_cell_weight_component);}

    opengl_world.Set_Key_Binding_Category("Pressure");
    filename=basedir+"/%d/pressure";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COLOR_MAP<T>* pressure_color_map=OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(101635,1e7,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1));
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* pressure_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,mac_grid,filename,pressure_color_map);
        Add_Component(pressure_component,"Pressure",'7',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',pressure_component->viewer_callbacks.Get("toggle_color_map"));
        opengl_world.Append_Bind_Key('`',pressure_component->viewer_callbacks.Get("toggle_smooth_slice"));
        slice_manager.Add_Object(pressure_component);}

    filename=basedir+"/%d/pressure_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* pressure_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(stream_type,grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1000,OPENGL_COLOR::Magenta(0,0),OPENGL_COLOR::Magenta(1)));
        pressure_gradient_component->opengl_scalar_field.Update();
        Add_Component(pressure_gradient_component,"Pressure Gradient",'1',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        slice_manager.Add_Object(pressure_gradient_component);}

    opengl_world.Set_Key_Binding_Category("Velocity");
    { // regular and mac velocity TODO: kill node based stuff?
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* vector_velocity_component=0;
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>* mac_velocity_component=0;
        filename=basedir+"/%d/velocities";
        if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
            if(node_based){
                vector_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(stream_type,regular_grid,filename);
                vector_velocity_component->opengl_grid_based_vector_field.size=0.01;
                vector_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
                Add_Component(vector_velocity_component,"Node velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
                slice_manager.Add_Object(vector_velocity_component);}
            else{
                mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>(stream_type,mac_grid,filename);
                mac_velocity_component->opengl_mac_velocity_field.size=0.1;
                mac_velocity_component->opengl_mac_velocity_field.vector_color=OPENGL_COLOR::Green();
                Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
                slice_manager.Add_Object(mac_velocity_component);}}
        
        filename=basedir+"/%d/mac_velocities";
        if(has_valid_grid && (always_add_mac_velocities || FILE_UTILITIES::Frame_File_Exists(filename,start_frame))){
            mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>(stream_type,mac_grid,filename);
            mac_velocity_component->opengl_mac_velocity_field.size=0.01;
            mac_velocity_component->opengl_mac_velocity_field.vector_color=OPENGL_COLOR::Magenta();
            mac_velocity_component->opengl_mac_velocity_field.Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_3D<T>::FACE_CENTERED);
            Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            slice_manager.Add_Object(mac_velocity_component);}
        
        if(mac_velocity_component){
            opengl_world.Append_Bind_Key('C',mac_velocity_component->viewer_callbacks.Get("toggle_draw"));
            opengl_world.Append_Bind_Key('V',mac_velocity_component->viewer_callbacks.Get("toggle_velocity_mode_and_draw"));
            opengl_world.Append_Bind_Key('v',mac_velocity_component->viewer_callbacks.Get("toggle_draw_vorticity"));
            opengl_world.Append_Bind_Key('N',mac_velocity_component->viewer_callbacks.Get("normalize_vorticity_color_map"));
            opengl_world.Append_Bind_Key('=',mac_velocity_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',mac_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
            opengl_world.Append_Bind_Key('h',mac_velocity_component->viewer_callbacks.Get("toggle_arrowhead"));}
        if(vector_velocity_component){
            opengl_world.Append_Bind_Key('v',vector_velocity_component->viewer_callbacks.Get("toggle_draw"));
            opengl_world.Append_Bind_Key('=',vector_velocity_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',vector_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
            opengl_world.Append_Bind_Key('h',vector_velocity_component->viewer_callbacks.Get("toggle_arrowhead"));}

        filename=basedir+"/%d/centered_velocities";
        if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
            OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* center_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(stream_type,mac_grid,filename);
            center_velocity_component->opengl_grid_based_vector_field.size=0.1;
            center_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(center_velocity_component,"Centered velocities",'B',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key('=',center_velocity_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',center_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
            slice_manager.Add_Object(center_velocity_component);}
    }

    filename=basedir+"/%d/object_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* object_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(stream_type,regular_grid,filename);
        object_velocity_component->opengl_grid_based_vector_field.size=0.01;
        object_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Yellow();
        Add_Component(object_velocity_component,"Object velocities",'%',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',object_velocity_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',object_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',object_velocity_component->viewer_callbacks.Get("toggle_arrowhead"));
        slice_manager.Add_Object(object_velocity_component);}

    filename=basedir+"/%d/forces";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>* force_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>(stream_type,mac_grid,filename);
        force_component->opengl_mac_velocity_field.size=0.01;
        force_component->opengl_mac_velocity_field.vector_color=OPENGL_COLOR::Yellow();
        force_component->opengl_mac_velocity_field.Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_3D<T>::FACE_CENTERED);
        Add_Component(force_component,"fluid control force",'G',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',force_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',force_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',force_component->viewer_callbacks.Get("toggle_arrowhead"));
        slice_manager.Add_Object(force_component);}

    filename=basedir+"/%d/velocities_ghost_fuel";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* vector_velocity_ghost_plus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(stream_type,regular_grid,filename);
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.size=0.01;
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(vector_velocity_ghost_plus_component,"Node ghost plus velocities",'b',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_plus_component->viewer_callbacks.Get("toggle_arrowhead"));
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_plus_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_plus_component->viewer_callbacks.Get("decrease_vector_size"));
            slice_manager.Add_Object(vector_velocity_ghost_plus_component);}}

    filename=basedir+"/%d/velocities_ghost";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* vector_velocity_ghost_minus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(stream_type,regular_grid,filename);
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.size=0.01;
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            vector_velocity_ghost_minus_component->Set_Draw(false);
            Add_Component(vector_velocity_ghost_minus_component,"Node ghost minus velocities",'c',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_minus_component->viewer_callbacks.Get("toggle_arrowhead"));
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_minus_component->viewer_callbacks.Get("increase_vector_size"));
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_minus_component->viewer_callbacks.Get("decrease_vector_size"));
            slice_manager.Add_Object(vector_velocity_ghost_minus_component);}}

    filename=basedir+"/%d/center_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame) && parameter_list.template Get_Parameter<bool>("use_cell_centered_velocities",false)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* center_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(stream_type,mac_grid,filename);
        center_velocity_component->opengl_grid_based_vector_field.size=0.1;
        center_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(center_velocity_component,"Centered velocities",'V',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',center_velocity_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',center_velocity_component->viewer_callbacks.Get("decrease_vector_size"));
        slice_manager.Add_Object(center_velocity_component);}

    opengl_world.Set_Key_Binding_Category("Pressure Jump");
    filename=basedir+"/%d/pressure_jumps";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* pressure_jump_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(stream_type,grid,filename);
        pressure_jump_component->opengl_grid_based_vector_field.size=0.001;
        pressure_jump_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Magenta();
        Add_Component(pressure_jump_component,"Pressure jumps",'&',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',pressure_jump_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',pressure_jump_component->viewer_callbacks.Get("decrease_vector_size"));
        opengl_world.Append_Bind_Key('h',pressure_jump_component->viewer_callbacks.Get("toggle_arrowhead"));
        slice_manager.Add_Object(pressure_jump_component);}

    filename=basedir+"/%d/beta_face";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/beta_face",start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T>* beta_face_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T>(stream_type,mac_grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,(T).002,OPENGL_COLOR::Gray(1),OPENGL_COLOR::Gray(0)));
        Add_Component(beta_face_component,"Beta Face",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F6),beta_face_component->viewer_callbacks.Get("toggle_draw"));
        slice_manager.Add_Object(beta_face_component);}

    opengl_world.Set_Key_Binding_Category("Pressure Jump");

    filename=basedir+"/%d/pseudo_dirichlet";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>* pseudo_dirichlet_component=new OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>(stream_type,grid,filename);
        pseudo_dirichlet_component->Set_Vector_Size(0.1);
        Add_Component(pseudo_dirichlet_component,"pseudo dirichlet",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F3),pseudo_dirichlet_component->viewer_callbacks.Get("toggle_draw"));
        opengl_world.Append_Bind_Key('=',pseudo_dirichlet_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',pseudo_dirichlet_component->viewer_callbacks.Get("decrease_vector_size"));
        slice_manager.Add_Object(pseudo_dirichlet_component);}

    filename=basedir+"/%d/thin_shells_grid_visibility";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>* thin_shells_debugging_component=new OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>(stream_type,grid,basedir);
        Add_Component(thin_shells_debugging_component,"thin shells debugging",'\0',BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F4),thin_shells_debugging_component->viewer_callbacks.Get("toggle_draw_grid_visibility_mode"));
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),thin_shells_debugging_component->viewer_callbacks.Get("toggle_draw_density_valid_mask"));
        slice_manager.Add_Object(thin_shells_debugging_component);}

    if(has_valid_grid){
        OPENGL_GRID_3D<T>* opengl_grid=new OPENGL_GRID_3D<T>(stream_type,*(new GRID<TV>(grid)),OPENGL_COLOR::Gray(0.5));
        opengl_grid->owns_grid=true;
        grid_component=new OPENGL_COMPONENT_BASIC<T,OPENGL_GRID_3D<T> >(stream_type,*opengl_grid);
        opengl_world.Set_Key_Binding_Category("Grid");
        Add_Component(grid_component,"Grid",'6',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('^',grid_component->object.viewer_callbacks.Get("toggle_draw_ghost_values"));
        slice_manager.Add_Object(grid_component);}

    opengl_world.Set_Key_Binding_Category("Particles");

    bool particles_stored_per_cell_uniform=false;
        if(has_valid_grid) particles_stored_per_cell_uniform=true;
    filename=basedir+"/%d/positive_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)||FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/positive_particles_0",start_frame)){
        positive_particles_component=new OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename,basedir+"/positive_particles_%d.%d",true,particles_stored_per_cell_uniform);
        positive_particles_component->opengl_points->color=OPENGL_COLOR(1,0.5,0);
        positive_particles_component->opengl_points->point_size=2;
        Add_Component(positive_particles_component,"Positive particles",'1',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('!',positive_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        if(slice_manager.slice) slice_manager.Add_Object(positive_particles_component);}

    filename=basedir+"/%d/negative_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)||FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/negative_particles_0",start_frame)){
        negative_particles_component=new OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename,basedir+"/negative_particles_%d.%d",true,particles_stored_per_cell_uniform);
        negative_particles_component->opengl_points->color=OPENGL_COLOR(0,0.5,1);
        negative_particles_component->opengl_points->point_size=2;
        Add_Component(negative_particles_component,"Negative particles",'2',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('@',negative_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        if(slice_manager.slice) slice_manager.Add_Object(negative_particles_component);}

    filename=basedir+"/%d/removed_positive_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)||FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/removed_positive_particles_0",start_frame)){
        removed_positive_particles_component=new OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename,basedir+"/removed_positive_particles_%d.%d",true,particles_stored_per_cell_uniform);
        removed_positive_particles_component->opengl_points->color=OPENGL_COLOR::Green();
        removed_positive_particles_component->opengl_points->point_size=2;
        Add_Component(removed_positive_particles_component,"Removed positive particles",'3',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('#',removed_positive_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        if(slice_manager.slice) slice_manager.Add_Object(removed_positive_particles_component);}

    filename=basedir+"/%d/removed_negative_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)||FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/removed_negative_particles_0",start_frame)){
        removed_negative_particles_component=new OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename,basedir+"/removed_negative_particles_%d.%d",true,particles_stored_per_cell_uniform);
        removed_negative_particles_component->opengl_points->color=OPENGL_COLOR::Cyan();
        removed_negative_particles_component->opengl_points->point_size=2;
        Add_Component(removed_negative_particles_component,"Removed negative particles",'4',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('$',removed_negative_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        if(slice_manager.slice) slice_manager.Add_Object(removed_negative_particles_component);}

    filename=basedir+"/%d/spray_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_3D<T>* spray_particles_component=new OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename);
        spray_particles_component->opengl_points->color=OPENGL_COLOR((T).7,(T).8,(T).9);
        spray_particles_component->opengl_points->point_size=1;
        Add_Component(spray_particles_component,"Spray particles",'5',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('%',spray_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        if(slice_manager.slice) slice_manager.Add_Object(spray_particles_component);}

    filename=basedir+"/%d/foam_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_3D<T>* foam_particles_component=new OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename);
        foam_particles_component->opengl_points->color=OPENGL_COLOR(1,0,0);
        foam_particles_component->opengl_points->point_size=1;
        Add_Component(foam_particles_component,"foam particles",'$',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('%',foam_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        if(slice_manager.slice) slice_manager.Add_Object(foam_particles_component);}

    filename=basedir+"/%d/vorticity_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<T>* vorticity_particles_component=new OPENGL_COMPONENT_VORTICITY_PARTICLES_3D<T>(stream_type,filename,true);
        vorticity_particles_component->opengl_points->color=OPENGL_COLOR(155,155,200);
        vorticity_particles_component->Set_Vector_Size((T).01);
        Add_Component(vorticity_particles_component,"Vorticity particles",'j',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::SELECTABLE);
        opengl_world.Append_Bind_Key('J',vorticity_particles_component->viewer_callbacks.Get("toggle_draw_point_numbers"));
        opengl_world.Append_Bind_Key('h',vorticity_particles_component->viewer_callbacks.Get("toggle_arrowhead"));
        opengl_world.Append_Bind_Key("^j",vorticity_particles_component->viewer_callbacks.Get("toggle_draw_velocities"));
        opengl_world.Append_Bind_Key('=',vorticity_particles_component->viewer_callbacks.Get("increase_vector_size"));
        opengl_world.Append_Bind_Key('-',vorticity_particles_component->viewer_callbacks.Get("decrease_vector_size"));
        if(slice_manager.slice) slice_manager.Add_Object(vorticity_particles_component);}

    filename=basedir+"/%d/sph_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_3D<T>* sph_particles_component=new OPENGL_COMPONENT_PARTICLES_3D<T>(stream_type,filename,"",false);
        sph_particles_component->opengl_points->color=OPENGL_COLOR::Blue();
        sph_particles_component->opengl_points->point_size=1;
        Add_Component(sph_particles_component,"SPH particles",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::SELECTABLE);
        if(slice_manager.slice) slice_manager.Add_Object(sph_particles_component);}

    if(has_valid_grid){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>* psi_N_component=0;
        if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N",start_frame))
            psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>(stream_type,grid,basedir+"/%d/psi_N",new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()));
        if(psi_N_component){
            Add_Component(psi_N_component,"Psi_N points",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
            opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_N_component->viewer_callbacks.Get("toggle_draw"));
            slice_manager.Add_Object(psi_N_component);}}

    filename=basedir+"/%d/psi_D";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>* psi_D_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>(stream_type,mac_grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS);
        Add_Component(psi_D_component,"Psi_D points",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_D_component->viewer_callbacks.Get("toggle_draw"));
        slice_manager.Add_Object(psi_D_component);}

    filename=basedir+"/%d/maccormack_cell_mask";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>* maccormack_cell_mask_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>(stream_type,mac_grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS);
        Add_Component(maccormack_cell_mask_component,"Maccormack cell mask points",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F9),maccormack_cell_mask_component->viewer_callbacks.Get("toggle_draw"));
        slice_manager.Add_Object(maccormack_cell_mask_component);}

    filename=basedir+"/%d/maccormack_face_mask";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>* maccormack_face_mask_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>(stream_type,grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()));
        Add_Component(maccormack_face_mask_component,"Maccormack face mask points",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F10),maccormack_face_mask_component->viewer_callbacks.Get("toggle_draw"));
        slice_manager.Add_Object(maccormack_face_mask_component);}

    filename=basedir+"/%d/colors";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_INDEXED_COLOR_MAP* colors_color_map=OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map();colors_color_map->Set_Index_Mode(OPENGL_INDEXED_COLOR_MAP::PERIODIC);
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T,int>* psi_colors_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,int>(stream_type,mac_grid,filename,colors_color_map,OPENGL_SCALAR_FIELD_3D<T,int>::DRAW_POINTS);
        Add_Component(psi_colors_component,"Psi colors",'\0',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_colors_component->viewer_callbacks.Get("toggle_draw"));
        slice_manager.Add_Object(psi_colors_component);}

    filename=basedir+"/%d/strain";
    OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T>* strain_component=0; // TODO: make this not a hack for multiphase
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T>(stream_type,grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->viewer_callbacks.Get("increase_size"));
        opengl_world.Append_Bind_Key('_',strain_component->viewer_callbacks.Get("decrease_size"));
        slice_manager.Add_Object(strain_component);}

    filename=basedir+"/%d/strain_0";
    if(!strain_component && has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T>(stream_type,grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION<T>::OWNED|BASIC_VISUALIZATION<T>::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->viewer_callbacks.Get("increase_size"));
        opengl_world.Append_Bind_Key('_',strain_component->viewer_callbacks.Get("decrease_size"));
        slice_manager.Add_Object(strain_component);}}

    filename=basedir+"/%d/surface.tri";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>* triangulated_surface_component=new OPENGL_COMPONENT_TRIANGULATED_SURFACE<T>(stream_type,filename);
        Add_Component(triangulated_surface_component,"Triangulated Surface",',',BASIC_VISUALIZATION<T>::OWNED);}

    if(slice_manager.slice){
        opengl_world.Set_Key_Binding_Category("Slice Control");
        opengl_world.Append_Bind_Key("^h",{[this](){slice_manager.Toggle_Slice_Mode();},"Toggle 3D/Slice mode"});
        opengl_world.Append_Bind_Key('\\',{[this](){slice_manager.Toggle_Slice_Axis();},"Toggle slice axis"});
        opengl_world.Append_Bind_Key(']',{[this](){slice_manager.Increment_Slice();},"Increment slice"});
        opengl_world.Append_Bind_Key('[',{[this](){slice_manager.Decrement_Slice();},"Decrement slice"});}

    opengl_world.Set_Key_Binding_Category("Misc.");
    opengl_world.Append_Bind_Key('~',{[this](){Command_Prompt();},"command_prompt"});

    opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),draw_all_objects_cb);

    // initialize selection priority (highest on top)
    Selection_Priority(OPENGL_SELECTION<T>::DEBUG_PARTICLES_3D)=110;
    Selection_Priority(OPENGL_SELECTION<T>::POINTS_3D)=100;
    Selection_Priority(OPENGL_SELECTION<T>::COMPONENT_PARTICLES_3D)=100;
    Selection_Priority(OPENGL_SELECTION<T>::ARTICULATED_RIGID_BODIES_JOINT_3D)=95;
    Selection_Priority(OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_VERTEX)=90;
    Selection_Priority(OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_SEGMENT)=89;
    Selection_Priority(OPENGL_SELECTION<T>::TRIANGULATED_SURFACE_TRIANGLE)=88;
    Selection_Priority(OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_VERTEX)=85;
    Selection_Priority(OPENGL_SELECTION<T>::TETRAHEDRALIZED_VOLUME_TETRAHEDRON)=84;
    Selection_Priority(OPENGL_SELECTION<T>::COMPONENT_RIGID_BODIES_3D)=80;
    Selection_Priority(OPENGL_SELECTION<T>::COMPONENT_DEFORMABLE_COLLECTION_3D)=80;
    Selection_Priority(OPENGL_SELECTION<T>::SEGMENTED_CURVE_VERTEX_3D)=79;
    Selection_Priority(OPENGL_SELECTION<T>::SEGMENTED_CURVE_SEGMENT_3D)=78;
    Selection_Priority(OPENGL_SELECTION<T>::GRID_NODE_3D)=70;
    Selection_Priority(OPENGL_SELECTION<T>::GRID_CELL_3D)=60;
    }

template<class T> void VISUALIZATION<T>::
Update_OpenGL_Strings()
{
    ANIMATED_VISUALIZATION<T>::Update_OpenGL_Strings();

    // TODO: slice manager should be a component
    std::ostringstream output_stream;
    if(slice_manager.slice) slice_manager.slice->Print_Slice_Info(output_stream);
    opengl_world.Add_String(output_stream.str());
}

template<class T> void VISUALIZATION<T>::
Pre_Frame_Extra()
{
    if(grid_component) grid_component->Set_Frame(frame);
}

template<class T> void VISUALIZATION<T>::
Set_Frame_Extra()
{
    std::string filename=LOG::sprintf("%s/%d/frame_title",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){std::ifstream input(filename.c_str());getline(input,frame_title);}
    else frame_title="";
    filename=LOG::sprintf("%s/%d/time",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){T time;FILE_UTILITIES::Read_From_File(stream_type,filename,time);frame_title=LOG::sprintf("(%.05f) ",time)+frame_title;}
}

template<class T> void VISUALIZATION<T>::
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
                if(!selection && positive_particles_component) selection=positive_particles_component->Get_Selection_By_Id(id);
                if(!selection && negative_particles_component) selection=negative_particles_component->Get_Selection_By_Id(id);
                if(!selection && removed_positive_particles_component) selection=removed_positive_particles_component->Get_Selection_By_Id(id);
                if(!selection && removed_negative_particles_component) selection=removed_negative_particles_component->Get_Selection_By_Id(id);
                if(selection){Set_Current_Selection(selection); Update_OpenGL_Strings();}}}
        else if(command=="n"){
            int id;
            if(sstream>>id){
                OPENGL_SELECTION<T>* selection=0;
                if(selection){Set_Current_Selection(selection); Update_OpenGL_Strings();}}}
        else if(command=="c"){
            int id;
            if(sstream>>id){
                OPENGL_SELECTION<T>* selection=0;
                if(selection){Set_Current_Selection(selection); Update_OpenGL_Strings();}}}
    }
}

template<class T> void VISUALIZATION<T>::
Command_Prompt()
{
    opengl_world.Prompt_User("Command: ",{[this](){Command_Prompt_Response();},"command_prompt_response"});
}

template<class T> void VISUALIZATION<T>::
Slice_Has_Changed()
{
    Update_OpenGL_Strings();
}

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
        ANIMATED_VISUALIZATION<float> *visualization=new VISUALIZATION<float>(STREAM_TYPE(1.f));
        visualization->Initialize_And_Run(parse_args);
        delete visualization;
    }
    else
    {
        ANIMATED_VISUALIZATION<double> *visualization=new VISUALIZATION<double>(STREAM_TYPE(1.));
        visualization->Initialize_And_Run(parse_args);
        delete visualization;
    }
    return 0;
}
