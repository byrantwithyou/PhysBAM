//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Log/LOG.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Core/Utilities/VIEWER_DIR.h>
#include <Grid_Tools/Computations/VORTICITY_UNIFORM.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D(const VIEWER_DIR& viewer_dir,const GRID<TV> &grid,const std::string &velocity_filename_input)
    :OPENGL_COMPONENT<T>(viewer_dir,"MAC Velocity Field 2D"),draw_vorticity(false),
    velocity_filename(velocity_filename_input),valid(false),draw_divergence(false),
    draw_streamlines(false),use_seed_for_streamlines(false),opengl_divergence_field(0),
    opengl_streamlines(streamlines),
    min_vorticity(FLT_MAX),max_vorticity(FLT_MIN)
{
    viewer_callbacks.Set("toggle_velocity_mode",{[this](){Toggle_Velocity_Mode();},"Toggle velocity mode"});
    viewer_callbacks.Set("toggle_velocity_mode_and_draw",{[this](){Toggle_Velocity_Mode_And_Draw();},"Toggle velocity mode and draw"});
    viewer_callbacks.Set("increase_vector_size",{[this](){Increase_Vector_Size();},"Increase vector size"});
    viewer_callbacks.Set("decrease_vector_size",{[this](){Decrease_Vector_Size();},"Decrease vector size"});
    viewer_callbacks.Set("toggle_arrowhead",{[this](){Toggle_Arrowhead();},"Toggle arrowhead"});
    viewer_callbacks.Set("toggle_draw_divergence",{[this](){Toggle_Draw_Divergence();},"Toggle draw divergence"});
    viewer_callbacks.Set("toggle_draw_streamlines",{[this](){Toggle_Draw_Streamlines();},"Toggle draw streamlines"});
    viewer_callbacks.Set("toggle_use_streamline_seed",{[this](){Toggle_Use_Streamline_Seed();},"Toggle draw consistent streamlines"});
    viewer_callbacks.Set("lengthen_streamlines",{[this](){Lengthen_Streamlines();},"Lengthen streamlines"});
    viewer_callbacks.Set("shorten_streamlines",{[this](){Shorten_Streamlines();},"Shorten streamlines"});
    viewer_callbacks.Set("toggle_draw_vorticity",{[this](){Toggle_Draw_Vorticity();},"Toggle draw vorticity"});
    viewer_callbacks.Set("normalize_vorticity_color_map",{[this](){Normalize_Vorticity_Color_Map();},"Normalize vorticity map based on current frame"});

    opengl_mac_velocity_field=new OPENGL_MAC_VELOCITY_FIELD_2D<T>(grid);
    number_of_steps=2*opengl_mac_velocity_field->grid.counts.x;
    opengl_vorticity_magnitude=new OPENGL_SCALAR_FIELD_2D<T>(opengl_mac_velocity_field->grid,*new ARRAY<T,VECTOR<int,2> >,OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1),"vorticity_magnitude");

    OPENGL_COLOR_RAMP<T>* ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color(-1e+2,OPENGL_COLOR::Red());
    ramp->Add_Color(-1,OPENGL_COLOR::Yellow());
    ramp->Add_Color(-1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(0,OPENGL_COLOR::Black());
    ramp->Add_Color(1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(1,OPENGL_COLOR::Yellow());
    ramp->Add_Color(1e+2,OPENGL_COLOR::Red());
    opengl_divergence_field=new OPENGL_SCALAR_FIELD_2D<T>(opengl_mac_velocity_field->grid,divergence,ramp,"divergence");

    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
~OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D()
{
    delete opengl_mac_velocity_field;
    delete opengl_divergence_field;
    delete &opengl_vorticity_magnitude->values;
    delete opengl_vorticity_magnitude;
}
//#####################################################################
// Function Print_Cell_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const
{
    if(!valid) return;
    stream<<component_name<<": "<<std::endl;
    opengl_mac_velocity_field->Print_Cell_Selection_Info(stream,cell);
    if(draw_vorticity)
        stream<<"vorticity magnitude = "<<opengl_vorticity_magnitude->values(cell)<<std::endl;
}
//#####################################################################
// Function Print_Node_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const
{
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Set_Frame()
{
    
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT<T>::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Display() const
{
    if(valid){
        if(draw) opengl_mac_velocity_field->Display();
        if(draw_divergence) opengl_divergence_field->Display();
        if(draw_vorticity) opengl_vorticity_magnitude->Display();
        if(draw_streamlines) opengl_streamlines.Display();}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Bounding_Box() const
{
    if(valid && draw) return opengl_mac_velocity_field->Bounding_Box();
    return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Reinitialize()
{
    if(!draw && !draw_divergence) return;
    valid = false;
    std::string tmp_filename=viewer_dir.current_directory+"/"+velocity_filename;
    if(File_Exists(tmp_filename))
        Read_From_File(tmp_filename,opengl_mac_velocity_field->face_velocities);
    else return;
    opengl_mac_velocity_field->Update();
    valid=true;
    Update_Divergence();
    Update_Streamlines();
    if(draw_vorticity){
        Update_Vorticity();
        opengl_vorticity_magnitude->Update();}
}
//#####################################################################
// Function Update_Divergence
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Update_Divergence()
{
    if(draw_divergence && valid){
        GRID<TV>& grid=opengl_mac_velocity_field->grid;
        ARRAY_VIEW<T,VECTOR<int,2> > &u=opengl_mac_velocity_field->u,&v=opengl_mac_velocity_field->v;
        static ARRAY<bool,FACE_INDEX<TV::m> > psi_N;
        static ARRAY<bool,TV_INT> psi_D;
        bool got_all_psi=true;
        std::string psi_N_filename=viewer_dir.current_directory+"/psi_N";
        std::string psi_D_filename=viewer_dir.current_directory+"/psi_D";
        if(File_Exists(psi_N_filename)) Read_From_File(psi_N_filename,psi_N);
        else got_all_psi=false;
        if(File_Exists(psi_D_filename)) Read_From_File(psi_D_filename,psi_D);
        else got_all_psi=false;
        if(!got_all_psi){psi_N.Clean_Memory();psi_D.Clean_Memory();}
        divergence.Resize(grid.counts);
        for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++){
            if(got_all_psi && (psi_D(i,j) || (psi_N.Component(0)(i,j) && psi_N.Component(0)(i+1,j) && psi_N.Component(1)(i,j) && psi_N.Component(1)(i,j+1)))) divergence(i,j)=0;
            else divergence(i,j)=grid.one_over_dX.x*(u(i+1,j)-u(i,j))+grid.one_over_dX.y*(v(i,j+1)-v(i,j));}
        opengl_divergence_field->Update();}
    else divergence.Clean_Memory();
}
//#####################################################################
// Function Update_Streamlines
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Update_Streamlines()
{
    streamlines.Clean_Memory();streamlines.particles.Clean_Memory();streamlines.mesh.Clean_Memory();
    if(!draw_streamlines || !valid) return;
    
    GRID<TV>& grid=opengl_mac_velocity_field->grid;
    int number_of_streamlines=100;
    T step_length=(T).5*grid.dX.Min();

    RANDOM_NUMBERS<T> random;
    if(use_seed_for_streamlines) random.Set_Seed(streamline_seed);
    LINEAR_INTERPOLATION_UNIFORM<TV,T> linear_interpolation;
    ARRAY<T,FACE_INDEX<TV::m> > mac_velocity_field(grid);
    mac_velocity_field.Component(0)=opengl_mac_velocity_field->u;
    mac_velocity_field.Component(1)=opengl_mac_velocity_field->v;
    FACE_LOOKUP_UNIFORM<TV> V_lookup(mac_velocity_field);

    for(int i=0;i<number_of_streamlines;i++){
        int p=streamlines.particles.Add_Element();
        TV X;
        random.Fill_Uniform(X,grid.domain);
        streamlines.particles.X(p)=X;
        for(int step=0;step<number_of_steps;step++){
            TV velocity=linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,X);
            TV X_new=X+step_length*velocity;
            velocity=(T).5*(velocity+linear_interpolation.Clamped_To_Array_Face(grid,V_lookup,X_new));
            X_new=grid.Clamp(X+step_length*velocity);
            int new_particle=streamlines.particles.Add_Element();
            streamlines.particles.X(new_particle)=X_new;
            streamlines.mesh.elements.Append(VECTOR<int,2>(p,new_particle));
            p=new_particle;
            X=X_new;}}
}
//#####################################################################
// Function Toggle_Velocity_Mode
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Velocity_Mode()
{
    opengl_mac_velocity_field->Toggle_Velocity_Mode();
}
//#####################################################################
// Function Toggle_Velocity_Mode_And_Draw
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Velocity_Mode_And_Draw()
{
    if(draw)
    {
        Toggle_Velocity_Mode();
        if((int)opengl_mac_velocity_field->velocity_mode==0) Set_Draw(!draw);
    }
    else Set_Draw(!draw);
}
//#####################################################################
// Function Increase_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Increase_Vector_Size()
{
    opengl_mac_velocity_field->Scale_Vector_Size(1.1);
}
//#####################################################################
// Function Decrease_Vector_Size
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Decrease_Vector_Size()
{
    opengl_mac_velocity_field->Scale_Vector_Size(1/1.1);
}
//#####################################################################
// Function Toggle_Arrowhead
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Arrowhead()
{
    opengl_mac_velocity_field->draw_arrowhead = !opengl_mac_velocity_field->draw_arrowhead;
}
//#####################################################################
// Function Toggle_Draw_Divergence
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Draw_Divergence()
{
    draw_divergence=!draw_divergence;
    Update_Divergence();
}
//#####################################################################
// Function Toggle_Draw_Streamlines
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Draw_Streamlines()
{
    draw_streamlines=!draw_streamlines;
    Update_Streamlines();
}
//#####################################################################
// Function Toggle_Use_Streamline_Seed
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Use_Streamline_Seed()
{
    use_seed_for_streamlines=!use_seed_for_streamlines;
}
//#####################################################################
// Function Set_Streamline_Seed
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Set_Streamline_Seed(const unsigned int seed)
{
    streamline_seed=seed;
}
//#####################################################################
// Function Lengthen_Streamlines
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Lengthen_Streamlines()
{
    number_of_steps+=10;
    Update_Streamlines();
}
//#####################################################################
// Function Shorten_Streamlines
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Shorten_Streamlines()
{
    number_of_steps=max(number_of_steps-10,0);
    Update_Streamlines();
}
//#####################################################################
// Function Toggle_Draw_Vorticity
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Draw_Vorticity()
{
    draw_vorticity=!draw_vorticity;
    if(draw_vorticity) valid=false;
    Reinitialize();
}
//#####################################################################
// Function Normalize_Vorticity_Color_Map
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Normalize_Vorticity_Color_Map()
{
    if(!draw_vorticity) return;
    opengl_vorticity_magnitude->Set_Scale_Range(min_vorticity,max_vorticity);
    valid=false;
    Reinitialize();
}
//#####################################################################
// Function Update_Vorticity
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Update_Vorticity()
{
    GRID<TV>& grid=opengl_mac_velocity_field->grid;
    RANGE<VECTOR<int,2> > domain_indices(grid.Domain_Indices());domain_indices.Change_Size(-VECTOR<int,2>::All_Ones_Vector());
    FACE_LOOKUP_UNIFORM<TV> lookup(opengl_mac_velocity_field->face_velocities);
    opengl_vorticity_magnitude->values.Resize(grid.Domain_Indices());
    for(CELL_ITERATOR<TV> iterator(grid,domain_indices);iterator.Valid();iterator.Next()){VECTOR<int,2> index=iterator.Cell_Index();
        T vorticity_magnitude=VORTICITY_UNIFORM<TV>::Vorticity(grid,lookup,index).Magnitude();
        opengl_vorticity_magnitude->values(index)=vorticity_magnitude;
        min_vorticity=min(min_vorticity,vorticity_magnitude);
        max_vorticity=max(max_vorticity,vorticity_magnitude);}
}
namespace PhysBAM{
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<float>;
template class OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<double>;
}
