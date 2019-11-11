//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <OpenGL/OpenGL/OPENGL_AXES.h>
#include <OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <OpenGL/OpenGL/OPENGL_GRID_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_1D.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace PhysBAM;

template<class T>
class OPENGL_1D_VISUALIZATION:public ANIMATED_VISUALIZATION<T>
{
    typedef VECTOR<T,1> TV;
public:
    using ANIMATED_VISUALIZATION<T>::opengl_window_title;
    using ANIMATED_VISUALIZATION<T>::start_frame;
    using ANIMATED_VISUALIZATION<T>::opengl_world;using ANIMATED_VISUALIZATION<T>::camera_script_filename;
    using ANIMATED_VISUALIZATION<T>::component_list;
    using ANIMATED_VISUALIZATION<T>::frame_title;using ANIMATED_VISUALIZATION<T>::Add_Component;
    OPENGL_1D_VISUALIZATION();
    ~OPENGL_1D_VISUALIZATION();

private:
    void Read_Grid();
protected:
    virtual void Add_Arguments(PARSE_ARGS& parse_args);
    virtual void Parse_Arguments(PARSE_ARGS& parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Add_OpenGL_Initialization();
    virtual void Pre_Frame_Extra();
    virtual void Set_Frame_Extra();

    // Options
    VIEWER_DIR viewer_dir{"."};

    bool has_valid_grid;

    bool node_based;

    GRID<TV> grid,mac_grid,regular_grid;
    OPENGL_COMPONENT_BASIC<T,OPENGL_GRID_1D<T> >* grid_component;
};

//#####################################################################
// Function OPENGL_1D_VISUALIZATION
//#####################################################################
template<class T> OPENGL_1D_VISUALIZATION<T>::
OPENGL_1D_VISUALIZATION()
    :ANIMATED_VISUALIZATION<T>(viewer_dir),grid_component(0)
{
    this->opengl_axes->visible=false;
}
//#####################################################################
// Function OPENGL_1D_VISUALIZATION
//#####################################################################
template<class T> OPENGL_1D_VISUALIZATION<T>::
~OPENGL_1D_VISUALIZATION()
{
    delete grid_component;
}
//#####################################################################
// Function Add_Arguments
//#####################################################################
template<class T> void OPENGL_1D_VISUALIZATION<T>::
Add_Arguments(PARSE_ARGS& parse_args)
{
    ANIMATED_VISUALIZATION<T>::Add_Arguments(parse_args);

    parse_args.Extra_Optional(&viewer_dir.output_directory,"basedir","base directory");
}
//#####################################################################
// Function Parse_Arguments
//#####################################################################
template<class T> void OPENGL_1D_VISUALIZATION<T>::
Parse_Arguments(PARSE_ARGS& parse_args)
{
#ifdef __linux__
    opengl_window_title="opengl_1d: " + Real_Path(viewer_dir.output_directory);
#endif

    ANIMATED_VISUALIZATION<T>::Parse_Arguments(parse_args);

    if(parse_args.unclaimed_arguments)
        parse_args.Print_Usage(true);

    // don't override camera script filename if it was already set in base class based on command line argument
    if(camera_script_filename.empty()) camera_script_filename=viewer_dir.output_directory+"/camera_script";
}
//#####################################################################
// Function Read_Grid
//#####################################################################
template<class T> void OPENGL_1D_VISUALIZATION<T>::
Read_Grid()
{
    has_valid_grid=false;
    if(File_Exists(viewer_dir.output_directory+"/common/grid")){
        std::string filename=viewer_dir.output_directory+"/common/grid";
        LOG::cout<<"Reading grid from '"<<filename<<"'..."<<std::endl;
        Read_From_File(filename,grid);
        has_valid_grid=true;}
    if(has_valid_grid){
        if(!grid.MAC_offset){
            regular_grid=grid;
            mac_grid.Initialize(grid.counts-1,grid.Domain(),true);
            node_based=true;}
        else{
            mac_grid=grid;
            regular_grid.Initialize(grid.counts+1,grid.Domain(),false);
            node_based=false;}
        LOG::cout<<"regular grid "<<regular_grid<<" mac grid "<<mac_grid<<" node_based "<<node_based<<std::endl;}
}
//#####################################################################
// Function Initialize_Components
//#####################################################################
template<class T> void OPENGL_1D_VISUALIZATION<T>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION<T>::Initialize_Components_And_Key_Bindings();
    std::string filename;

    Read_Grid();

    opengl_world.Unbind_Keys("delmo7vDEMO&V+-");


    opengl_world.Set_Key_Binding_Category("Compressible");

    if(has_valid_grid){
        OPENGL_GRID_1D<T>* opengl_grid=new OPENGL_GRID_1D<T>(viewer_dir,grid,OPENGL_COLOR::Gray(.5));
        grid_component=new OPENGL_COMPONENT_BASIC<T,OPENGL_GRID_1D<T> >(viewer_dir,*opengl_grid);
        Add_Component(grid_component,"Grid",'6',BASIC_VISUALIZATION<T>::SELECTABLE);}

    filename="u";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "u",'\0',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE);
    filename="rigid_body_particles";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>* rigid_bodies_component=
            new OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T>(viewer_dir);
        Add_Component(rigid_bodies_component,"Rigid Bodies",'5',BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key('n',rigid_bodies_component->viewer_callbacks.Get("toggle_show_object_names"));}
    filename="deformable_object_particles";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>* deformable_geometry_component=new OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_1D<T>(viewer_dir);
        Add_Component(deformable_geometry_component,"Deformable Bodies",'8',BASIC_VISUALIZATION<T>::OWNED);}
    filename="u_exact";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Cyan()),
            "u_exact",'\0',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE);
    filename="levelset";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_LEVELSET_1D<T>(viewer_dir,grid,filename,OPENGL_COLOR::Yellow(),OPENGL_COLOR::Yellow()),
            "levelset",'l',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE);
    // Compressible
    filename="density";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Yellow(),OPENGL_COLOR::Yellow()),
            "density",'1',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE);
    filename="centered_velocities";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Magenta(),OPENGL_COLOR::Magenta()),
            "velocity",'v',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE);
    filename="mac_velocities";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "mac_velocities",'a',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="momentum";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Magenta(),OPENGL_COLOR::Magenta()),
            "momentum",'2',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="energy";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan()),
            "energy",'3',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="internal_energy";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan()),
            "internal energy",'4',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="pressure";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan()),
            "pressure",'7',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE);
    filename="compressible_implicit_pressure";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>* compressible_implicit_pressure_component=
            new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red());
        Add_Component(compressible_implicit_pressure_component,
            "compressible_implicit_pressure",'9',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename="entropy";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Green()),
            "entropy",'e',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="speedofsound";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "speedofsound",'o',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="machnumber";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Magenta()),
            "machnumber",'m',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="velocity_plus_c";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Green()),
            "velocityplusc",'+',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="velocity_minus_c";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Blue()),
            "velocityminusc",'-',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
#if 0
    filename=basedir+"_exact/%d/density";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Cyan()),
            "density_exact",'!',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename=basedir+"_exact/%d/centered_velocities";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Yellow()),
            "velocity_exact",'V',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename=basedir+"_exact/%d/energy";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "energy_exact",'#',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename=basedir+"_exact/%d/internal_energy";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "internal energy_exact",'$',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename=basedir+"_exact/%d/pressure";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "pressure_exact",'&',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename=basedir+"_exact/%d/entropy";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Cyan()),
            "entropy_exact",'E',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename=basedir+"_exact/%d/speedofsound";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Yellow()),
            "speedofsound_exact",'O',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename=basedir+"_exact/%d/machnumber";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "machnumber_exact",'M',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);}
#endif
    filename="p_cavitation";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "p_cavitation",'[',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename="p_internal_energy";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "p_internal_energy",']',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename="density_flux";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        Add_Component(new OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "density_flux",'W',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename="momentum_flux";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        Add_Component(new OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Cyan()),
            "momentum_flux",'S',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);}
    filename="energy_flux";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        Add_Component(new OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Yellow()),
            "energy_flux",'X',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);}

    filename="psi_N";
    if(File_Exists(viewer_dir.current_directory+"/"+filename)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,bool>* psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,bool>(viewer_dir,grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan());
        Add_Component(psi_N_component,"Psi_N points",'\0',BASIC_VISUALIZATION<T>::START_HIDDEN|BASIC_VISUALIZATION<T>::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_N_component->viewer_callbacks.Get("toggle_draw"));}
    filename="psi_D";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,bool>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "psi_D",'D',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);
    filename="euler_psi";
    if(File_Exists(viewer_dir.current_directory+"/"+filename))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,bool>(viewer_dir,grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "euler_psi",'P',BASIC_VISUALIZATION<T>::OWNED | BASIC_VISUALIZATION<T>::SELECTABLE | BASIC_VISUALIZATION<T>::START_HIDDEN);


    // add scaling controls to scalar fields
    opengl_world.Set_Key_Binding_Category("Scaling");
    for(int c=0;c<component_list.m;c++){
        if(OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>* scalar_field_component=dynamic_cast<OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T>*>(component_list(c))){
            opengl_world.Append_Bind_Key('>',scalar_field_component->viewer_callbacks.Get("increase_scale"));
            opengl_world.Append_Bind_Key('<',scalar_field_component->viewer_callbacks.Get("decrease_scale"));}
        else if(OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T>* face_scalar_field_component=dynamic_cast<OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T>*>(component_list(c))){
            opengl_world.Append_Bind_Key('>',face_scalar_field_component->viewer_callbacks.Get("increase_scale"));
            opengl_world.Append_Bind_Key('<',face_scalar_field_component->viewer_callbacks.Get("decrease_scale"));}}
}
//#####################################################################
// Function Set_Frame_Extra
//#####################################################################
template<class T> void OPENGL_1D_VISUALIZATION<T>::
Set_Frame_Extra()
{
    std::string filename=viewer_dir.current_directory+"/frame_title";
    if(File_Exists(filename)){std::ifstream input(filename.c_str());getline(input,frame_title);}
    else frame_title="";
    filename=viewer_dir.current_directory+"/time";
    if(File_Exists(filename)){T time;Read_From_File(filename,time);frame_title=LOG::sprintf("(%.05f) ",time)+frame_title;}
}
//#####################################################################
// Function Pre_Frame_Extra
//#####################################################################
template<class T> void OPENGL_1D_VISUALIZATION<T>::
Pre_Frame_Extra()
{
}
//#####################################################################
// Function Add_OpenGL_Initialization
//#####################################################################
template<class T> void OPENGL_1D_VISUALIZATION<T>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION<T>::Add_OpenGL_Initialization();
    opengl_world.Set_2D_Mode(true);
}
//#####################################################################

int main(int argc,char *argv[])
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
        ANIMATED_VISUALIZATION<float> *visualization=new OPENGL_1D_VISUALIZATION<float>();
        visualization->Initialize_And_Run(parse_args);
        delete visualization;
    }
    else
    {
        ANIMATED_VISUALIZATION<double> *visualization=new OPENGL_1D_VISUALIZATION<double>();
        visualization->Initialize_And_Run(parse_args);
        delete visualization;
    }


    return 0;
}
