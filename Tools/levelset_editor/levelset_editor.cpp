//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Particles/PARTICLES_SUBSET.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/BASIC_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MOUSE_HANDLER.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE_MANAGER.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_3D.h>

#define SIGMA 0.5

using namespace PhysBAM;

template<class T>
class ATTACHMENT_VISUALIZATION:public BASIC_VISUALIZATION,public OPENGL_MOUSE_HANDLER
{
    typedef VECTOR<T,3> TV;
public:
    DEFORMABLE_PARTICLES<TV> particles;
    TRIANGLE_MESH tri_mesh;
    TRIANGULATED_SURFACE<T> tri_surface;
    OPENGL_GRID_3D<T>* opengl_grid_3d;
    OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> >* opengl_grid_3d_component;
    OPENGL_COMPONENT_LEVELSET_3D<T>* opengl_levelset_component;
    std::string data_directory;
    bool paint_mode;
    int paint_radius;
    ARRAY<PAIR<T,VECTOR<int,3> > >copy;
    ARRAY<PAIR<T,VECTOR<int,3> > >selection;
    int last_selection;
    bool dragging;
    bool is_valid;
    bool expanded;
    ARRAY<TV> paint_surface_points;
    ARRAY<TV> template_surface_points;
    GRID<TV> grid;
    ARRAY<T,VECTOR<int,3> > phi, copy_phi;
    LEVELSET_3D<GRID<TV> > levelset,copy_levelset;
    LEVELSET_ADVECTION_3D<T> levelset_advection;
    OPENGL_SELECTION_GRID_CELL_LIST_3D<T> opengl_selection;
    OPENGL_SLICE_MANAGER slice_manager;
    std::string selection_read_filename,selection_write_filename,levelset_copy_filename,levelset_write_filename;
    int min_x,min_y,min_z,max_x,max_y,max_z;
 
    void Clear_Selection()
    {is_valid=expanded=false;selection.Remove_All();opengl_grid_3d->Clear_Highlight();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Clear_Selection);

    void Update_OpenGL_Selection()
    {
        opengl_selection.indicies.Remove_All();
        if(selection.m==0) {opengl_grid_3d->Clear_Highlight();return;}
        for(int i=0;i<selection.m;i++) opengl_selection.indicies.Append(selection(i).y);
        opengl_grid_3d->Highlight_Selection(&opengl_selection);
    }

    void Update_OpenGL_Strings()
    {
        opengl_world.Clear_Strings();
        std::ostringstream output_stream;
        output_stream<<"Paint the Fence: "<<(paint_mode?"On":"Off")<<std::endl;
        output_stream<<"Paint Radius: "<<paint_radius<<std::endl;
        opengl_world.Add_String(output_stream.str());
    }

    void Slice_Has_Changed()
    {
        Update_OpenGL_Selection();
        Update_OpenGL_Strings();
    }
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Slice_Has_Changed);
    
    void Handle_Click(int button,int state,int x,int y,bool ctrl_pressed,bool shift_pressed)
    {
        if(state==GLUT_DOWN && paint_mode){
            if (dragging==false) last_selection=selection.m+1;
            dragging=true;}
        else if(state==GLUT_UP){dragging=false;}
    }

    void Handle_Drag(int x,int y)
    {
        std::cout<<"Dragging at "<<x<<","<<y<<std::endl;
        OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)slice_manager.slice;
        if(paint_mode && slice->mode==OPENGL_SLICE::NO_SLICE){
            is_valid=expanded=false;
            RAY<TV> ray=OPENGL_WORLD::Singleton()->Ray_Through_Normalized_Image_Coordinate(OPENGL_WORLD::Singleton()->Convert_Mouse_Coordinates(x,y));
            if(INTERSECTION::Intersects(ray,tri_surface)){
                std::cout<<"  intersected at "<<ray.Point(ray.t_max)<<std::endl;
                //store selection in the grid
                VECTOR<int,3> grid_index=grid.Index(ray.Point(ray.t_max));
                selection.Append_Unique(PAIR<T,VECTOR<int,3> >(phi(grid_index),grid_index));
                for(int i=1;i<paint_radius;i++){
                    if(abs(phi(grid_index(1)+i,grid_index(2),grid_index(3)))<grid.min_dX)
                        selection.Append_Unique(PAIR<T,VECTOR<int,3> >(phi(grid_index(1)+i,grid_index(2),grid_index(3)),VECTOR<int,3>(grid_index(1)+i,grid_index(2),grid_index(3))));
                    if(abs(phi(grid_index(1),grid_index(2)+i,grid_index(3)))<grid.min_dX)
                        selection.Append_Unique(PAIR<T,VECTOR<int,3> >(phi(grid_index(1),grid_index(2)+i,grid_index(3)),VECTOR<int,3>(grid_index(1),grid_index(2)+i,grid_index(3))));
                    if(abs(phi(grid_index(1),grid_index(2),grid_index(3)+i))<grid.min_dX)
                        selection.Append_Unique(PAIR<T,VECTOR<int,3> >(phi(grid_index(1),grid_index(2),grid_index(3)+i),VECTOR<int,3>(grid_index(1),grid_index(2),grid_index(3)+i)));
                    if(abs(phi(grid_index(1)-i,grid_index(2),grid_index(3)))<grid.min_dX)
                        selection.Append_Unique(PAIR<T,VECTOR<int,3> >(phi(grid_index(1)-i,grid_index(2),grid_index(3)),VECTOR<int,3>(grid_index(1)-i,grid_index(2),grid_index(3))));
                    if(abs(phi(grid_index(1),grid_index(2)-i,grid_index(3)))<grid.min_dX)
                        selection.Append_Unique(PAIR<T,VECTOR<int,3> >(phi(grid_index(1),grid_index(2)-i,grid_index(3)),VECTOR<int,3>(grid_index(1),grid_index(2)-i,grid_index(3))));
                    if(abs(phi(grid_index(1),grid_index(2),grid_index(3)-i))<grid.min_dX)
                        selection.Append_Unique(PAIR<T,VECTOR<int,3> >(phi(grid_index(1),grid_index(2),grid_index(3)-i),VECTOR<int,3>(grid_index(1),grid_index(2),grid_index(3)-i)));}
                Update_OpenGL_Selection();
            }
        }
    }

    ATTACHMENT_VISUALIZATION(const std::string& data_dir, const std::string& prefix)
        :tri_surface(tri_mesh,particles),paint_mode(false),dragging(false),is_valid(false),expanded(false),levelset(grid,phi),copy_levelset(grid,copy_phi),levelset_advection(&levelset),
        opengl_selection(opengl_grid_3d)
    {
        data_directory=data_dir+"/"+prefix;
        opengl_window_title="Levelset Editor";
        levelset_copy_filename=data_directory+"_copy.phi";
        levelset_write_filename=data_directory+"_out.phi";
        selection_read_filename=data_directory+".sel";
        selection_write_filename=data_directory+"_out.sel";
        paint_radius=0;

    }

    void Initialize_Components_And_Key_Bindings() PHYSBAM_OVERRIDE
    {FILE_UTILITIES::Read_From_File<T>(data_directory+".tri",tri_surface);
    FILE_UTILITIES::Read_From_File<T>(data_directory+".phi",levelset);
    copy_phi.Resize(grid.Domain_Indices());
        
    tri_surface.Initialize_Hierarchy();
    tri_surface.mesh.Initialize_Neighbor_Nodes();

    OPENGL_UNIFORM_SLICE* slice=new OPENGL_UNIFORM_SLICE(opengl_world);
    slice->Initialize(grid);
    slice_manager.slice=slice;
    slice_manager.Set_Slice_Has_Changed_Callback(Slice_Has_Changed_CB());
    
    opengl_levelset_component=new OPENGL_COMPONENT_LEVELSET_3D<T>("","","","");
    opengl_levelset_component->opengl_levelset_multiview->Set_Levelset(levelset);
    opengl_levelset_component->opengl_levelset_multiview->Generate_Triangulated_Surface(false,"");
    opengl_levelset_component->opengl_levelset_multiview->Initialize();
    opengl_world.Set_Key_Binding_Category("Levelset");
    Add_Component(opengl_levelset_component,"levelset",'l',BASIC_VISUALIZATION::OWNED);
    opengl_world.Append_Bind_Key('L',opengl_levelset_component->Toggle_Slice_Color_Mode_CB());
    opengl_world.Append_Bind_Key("^l",opengl_levelset_component->Toggle_Display_Overlay_CB());    
    opengl_world.Append_Bind_Key('`',opengl_levelset_component->Toggle_Smooth_Slice_CB());
    if(opengl_levelset_component->Use_Sets()){
        opengl_world.Append_Bind_Key('M',opengl_levelset_component->Toggle_Draw_Multiple_Levelsets_CB());
        opengl_world.Append_Bind_Key('>',opengl_levelset_component->Next_Set_CB());
        opengl_world.Append_Bind_Key('<',opengl_levelset_component->Previous_Set_CB());}
    slice_manager.Add_Object(opengl_levelset_component);
    
    opengl_grid_3d=new OPENGL_GRID_3D<T>(grid);
    opengl_grid_3d->Clear_Highlight();
    opengl_grid_3d->hide_non_selected_grid=true;
    
    opengl_grid_3d_component=new OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> >(*opengl_grid_3d);
    opengl_world.Set_Key_Binding_Category("Grid");
    Add_Component(opengl_grid_3d_component,"grid",'6',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE|BASIC_VISUALIZATION::START_HIDDEN);
    opengl_world.Append_Bind_Key('^',opengl_grid_3d_component->object.Toggle_Draw_Ghost_Values_CB());
    slice_manager.Add_Object(opengl_grid_3d_component);

    opengl_world.Set_Key_Binding_Category("Levelset Editor");
    opengl_world.Append_Bind_Key('p',Toggle_Paint_Mode_CB("Toggle Paint Mode"));
    opengl_world.Append_Bind_Key('a',Select_All_CB("Select All"));
    opengl_world.Append_Bind_Key('c',Copy_Selection_CB("Copy Selection"));
    opengl_world.Append_Bind_Key('u',Undo_Selection_CB("Copy Selection"));
    opengl_world.Append_Bind_Key('x',Cut_Selection_CB("Cut Selection"));
    opengl_world.Append_Bind_Key('v',Paste_Selection_CB("Paste Selection"));
    opengl_world.Append_Bind_Key('d',Delete_Selection_CB("Delete Selection"));
    opengl_world.Append_Bind_Key('f',Fill_Selection_CB("Fill Selection"));
    opengl_world.Append_Bind_Key('e',Expand_Selection_CB("Expand Selection"));
    opengl_world.Append_Bind_Key('i',Expand_Selection_Inside_CB("Expand Selection Inside Surface"));
    opengl_world.Append_Bind_Key('n',Evolve_Selection_CB("Evolve Selection Using Normals"));
    opengl_world.Append_Bind_Key('m',Smooth_Selection_CB("Smooth Selection"));
    opengl_world.Append_Bind_Key('c',Clear_Selection_CB("Clear Selection"));
    opengl_world.Append_Bind_Key('o',Load_Selection_CB("Load Selection"));
    opengl_world.Append_Bind_Key('s',Save_Selection_CB("Save Selection"));
    opengl_world.Append_Bind_Key('r',Save_Selection_As_CB("Save Selection as Levelset"));
    opengl_world.Append_Bind_Key('w',Write_Levelset_CB("Write Levelset"));
    opengl_world.Append_Bind_Key('+',Increase_Paint_Radius_CB("Increase Paint Radius"));
    opengl_world.Append_Bind_Key('-',Decrease_Paint_Radius_CB("Decrease Paint Radius"));
    opengl_world.Set_Key_Binding_Category("Slice Control");
    opengl_world.Append_Bind_Key("^h",slice_manager.Toggle_Slice_Mode_CB("Toggle 3D/Slice mode"));
    opengl_world.Append_Bind_Key('\\',slice_manager.Toggle_Slice_Axis_CB("Toggle slice axis"));
    opengl_world.Append_Bind_Key(']',slice_manager.Increment_Slice_CB("Increment slice"));
    opengl_world.Append_Bind_Key('[',slice_manager.Decrement_Slice_CB("Decrement slice"));}

    virtual ~ATTACHMENT_VISUALIZATION()
    {}

    void Calculate_Bounds()
    {expanded=true;
    min_x=selection(1).y.x;max_x=selection(1).y.x;min_y=selection(1).y.y;max_y=selection(1).y.y;min_z=selection(1).y.z;max_z=selection(1).y.z;
    for(int i=0;i<selection.m;i++){
        if(selection(i).y.x<min_x) min_x=selection(i).y.x;
        if(selection(i).y.x>max_x) max_x=selection(i).y.x;
        if(selection(i).y.y<min_y) min_y=selection(i).y.y;
        if(selection(i).y.y>max_y) max_y=selection(i).y.y;
        if(selection(i).y.z<min_z) min_z=selection(i).y.z;
        if(selection(i).y.z>max_z) max_z=selection(i).y.z;}}

    void Update_Current_Levelset()
    {opengl_levelset_component->opengl_levelset_multiview->Reset_Surface();
    opengl_levelset_component->opengl_levelset_multiview->Initialize();}

    void Create_Copy_Levelset()
    {if(is_valid) return;
    is_valid=true;
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x;i++) for(int j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y;j++) for(int ij=phi.domain.min_corner.z;ij<=phi.domain.max_corner.z;ij++) copy_phi(i,j,ij)=grid.min_dX;
    for(int i=0;i<selection.m;i++) copy_phi(selection(i).y)=selection(i).x;
    copy_levelset.Fast_Marching_Method(0,grid.dX.x*8);}

    void Toggle_Paint_Mode()
    {if(paint_mode){paint_mode=false;OPENGL_WORLD::Singleton()->Set_External_Mouse_Handler();}
    else{paint_mode=true;OPENGL_WORLD::Singleton()->Set_External_Mouse_Handler(this);}
    Update_OpenGL_Strings();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Toggle_Paint_Mode);

    void Select_All()
    {selection.Remove_All();
    for(int i=phi.domain.min_corner.x;i<=phi.domain.max_corner.x;i++) for(int j=phi.domain.min_corner.y;j<=phi.domain.max_corner.y;j++) for(int ij=phi.domain.min_corner.z;ij<=phi.domain.max_corner.z;ij++) if(phi(i,j,ij)<=0) selection.Append(PAIR<T,VECTOR<int,3> >(phi(i,j,ij),VECTOR<int,3>(i,j,ij)));
    Update_OpenGL_Selection();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Select_All);

    void Copy_Selection()
    {copy=selection;}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Copy_Selection);

    void Undo_Selection()
    {if(last_selection==-1) return;
    for(int i=selection.m;i>=last_selection;i--) selection.Remove_Index(i);
    last_selection=-1;
    is_valid=false;
    Update_OpenGL_Selection();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Undo_Selection);

    void Cut_Selection()
    {Copy_Selection();Delete_Selection();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Cut_Selection);

    void Paste_Selection()
    {for(int i=0;i<copy.m;i++) phi(copy(i).y)=copy(i).x;}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Paste_Selection);

    void Delete_Selection()
    {for(int i=0;i<selection.m;i++) phi(selection(i).y)=grid.min_dX;
    Update_Current_Levelset();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Delete_Selection);

    void Fill_Selection()
    {for(int i=0;i<selection.m;i++) phi(selection(i).y)=-grid.min_dX;
    Update_Current_Levelset();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Fill_Selection);

    void Expand_Selection()
    {Calculate_Bounds();
    selection.Remove_All();is_valid=false;
    for(int i=min_x;i<=max_x;i++) for(int j=min_y;j<=max_y;j++) for(int ij=min_z;ij<=max_z;ij++) if(abs(phi(i,j,ij))<grid.min_dX) selection.Append(PAIR<T,VECTOR<int,3> >(phi(i,j,ij),VECTOR<int,3>(i,j,ij)));
    Update_OpenGL_Selection();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Expand_Selection);

    void Expand_Selection_Inside()
    {Calculate_Bounds();
    selection.Remove_All();is_valid=false;
    for(int i=min_x;i<=max_x;i++) for(int j=min_y;j<=max_y;j++) for(int ij=min_z;ij<=max_z;ij++) if(phi(i,j,ij)<grid.min_dX) selection.Append(PAIR<T,VECTOR<int,3> >(phi(i,j,ij),VECTOR<int,3>(i,j,ij)));
    Update_OpenGL_Selection();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Expand_Selection_Inside);

    void Evolve_Selection()
    {PHYSBAM_NOT_IMPLEMENTED();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Evolve_Selection);

    void Smooth_Selection_Response()
    {levelset.Set_Curvature_Motion(SIGMA);
    ARRAY<T,FACE_INDEX<TV::dimension> > V(grid);
    T step=levelset.CFL(V);
    int iterations=1;
    if(!OPENGL_WORLD::Singleton()->prompt_response.empty()){
        std::istringstream sstream(OPENGL_WORLD::Singleton()->prompt_response);
        sstream>>iterations;}
    for(int i=0;i<iterations;i++){
        if(expanded) levelset_advection.Euler_Step_Subset(V,min_x,max_x,min_y,max_y,min_z,max_z,step,0,3);
        else for(int i=0;i<selection.m;i++) levelset_advection.Euler_Step_Cell(V,selection(i).y.x,selection(i).y.y,selection(i).y.z,step,0,3);}
    Update_Current_Levelset();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Smooth_Selection_Response);

    void Smooth_Selection()
    {OPENGL_WORLD::Singleton()->Prompt_User("Enter Number of Iterations: ",Smooth_Selection_Response_CB());}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Smooth_Selection);

    void Load_Selection()
    {FILE_UTILITIES::Read_From_File<T>(selection_read_filename,copy);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Load_Selection);
    
    void Save_Selection()
    {FILE_UTILITIES::Write_To_File<T>(selection_write_filename,selection);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Save_Selection);

    void Save_Selection_As()
    {Create_Copy_Levelset();FILE_UTILITIES::Write_To_File<T>(levelset_copy_filename,copy_levelset);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Save_Selection_As);

    void Write_Levelset()
    {levelset.Fast_Marching_Method(0,grid.dX.x*8);
    FILE_UTILITIES::Write_To_File<T>(levelset_write_filename,levelset);}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Write_Levelset);

    void Increase_Paint_Radius()
    {paint_radius++;
    Update_OpenGL_Strings();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Increase_Paint_Radius);

    void Decrease_Paint_Radius()
    {paint_radius=max(1,paint_radius-1);
    Update_OpenGL_Strings();}
    DEFINE_CALLBACK_CREATOR(ATTACHMENT_VISUALIZATION,Decrease_Paint_Radius);
};
//#####################################################################
// Function main
//#####################################################################
int main(int argc,char *argv[])
{
    std::string data_dir,file;
    if(argc<3) {
        LOG::cout << std::endl << "Usage: ./levelset_editor data_directory file" << std::endl;
        exit(1);
    }
    data_dir=argv[1];
    file=std::string(argv[2]);
    ATTACHMENT_VISUALIZATION<float> visualization(data_dir,file);
    // This going to mess with the args?
    int dummy_argc=1;const char *dummy_argv[]={"paint"};
    visualization.Initialize_And_Run(dummy_argc,(char**)dummy_argv);

    return 0;
}
