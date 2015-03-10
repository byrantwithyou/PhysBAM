//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Avi Robinson-Mosher, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/CELL_ITERATOR.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Grids_Uniform_Computations/VORTICITY_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
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
OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D(STREAM_TYPE stream_type,const GRID<TV> &grid,const std::string &velocity_filename_input)
    :OPENGL_COMPONENT<T>(stream_type,"MAC Velocity Field 2D"),draw_vorticity(false),
    velocity_filename(velocity_filename_input),valid(false),draw_divergence(false),
    draw_streamlines(false),use_seed_for_streamlines(false),opengl_divergence_field(0),
    opengl_streamlines(stream_type,streamlines),psi_N_psi_D_basedir(""),min_vorticity(FLT_MAX),max_vorticity(FLT_MIN)
{
    is_animation=FILE_UTILITIES::Is_Animated(velocity_filename);
    frame_loaded=-1;

    opengl_mac_velocity_field=new OPENGL_MAC_VELOCITY_FIELD_2D<T>(stream_type,grid);
    number_of_steps=2*opengl_mac_velocity_field->grid.counts.x;
    opengl_vorticity_magnitude=new OPENGL_SCALAR_FIELD_2D<T>(stream_type,opengl_mac_velocity_field->grid,*(new ARRAY<T,VECTOR<int,2> >),OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1));

    OPENGL_COLOR_RAMP<T>* ramp=new OPENGL_COLOR_RAMP<T>;
    ramp->Add_Color(-1e+2,OPENGL_COLOR::Red());
    ramp->Add_Color(-1,OPENGL_COLOR::Yellow());
    ramp->Add_Color(-1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(0,OPENGL_COLOR::Black());
    ramp->Add_Color(1e-2,OPENGL_COLOR::Green());
    ramp->Add_Color(1,OPENGL_COLOR::Yellow());
    ramp->Add_Color(1e+2,OPENGL_COLOR::Red());
    opengl_divergence_field=new OPENGL_SCALAR_FIELD_2D<T>(stream_type,opengl_mac_velocity_field->grid,divergence,ramp);

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
// Function Valid_Frame
//#####################################################################
template<class T> bool OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(velocity_filename, frame_input);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* selection) const
{
    if(Is_Up_To_Date(frame)){
        stream<<component_name<<": "<<std::endl;
        opengl_mac_velocity_field->Print_Selection_Info(stream,selection);
        if(draw_vorticity && ((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)){
            VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)->index;
            stream<<"vorticity magnitude = "<<opengl_vorticity_magnitude->values(index)<<std::endl;}}
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT<T>::Set_Frame(frame_input);
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
    else return RANGE<VECTOR<T,3> >::Centered_Box();
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T> void OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>::
Reinitialize()
{
    if(draw || draw_divergence){
        if((is_animation && frame_loaded!=frame) || (!is_animation && frame_loaded < 0)){
            valid = false;
            std::string tmp_filename=FILE_UTILITIES::Get_Frame_Filename(velocity_filename.c_str(), frame);
            if(FILE_UTILITIES::File_Exists(tmp_filename)) FILE_UTILITIES::Read_From_File(stream_type,tmp_filename,opengl_mac_velocity_field->face_velocities);
            else return;
            opengl_mac_velocity_field->Update();
            frame_loaded=frame;
            valid=true;
            Update_Divergence();
            Update_Streamlines();
            if(draw_vorticity){
                Update_Vorticity();
                opengl_vorticity_magnitude->Update();}
        }
    }
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
        if(!psi_N_psi_D_basedir.empty()){
            std::string psi_N_filename=STRING_UTILITIES::string_sprintf("%s/%d/psi_N",psi_N_psi_D_basedir.c_str(),frame);
            std::string psi_D_filename=STRING_UTILITIES::string_sprintf("%s/%d/psi_D",psi_N_psi_D_basedir.c_str(),frame);
            if(FILE_UTILITIES::File_Exists(psi_N_filename)) FILE_UTILITIES::Read_From_File(stream_type,psi_N_filename,psi_N);
            else got_all_psi=false;
            if(FILE_UTILITIES::File_Exists(psi_D_filename)) FILE_UTILITIES::Read_From_File(stream_type,psi_D_filename,psi_D);
            else got_all_psi=false;}
        else got_all_psi=false;
        if(!got_all_psi){psi_N.Clean_Memory();psi_D.Clean_Memory();}
        divergence.Resize(0,grid.counts.x,0,grid.counts.y);
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
        TV X=streamlines.particles.X(p)=random.Get_Uniform_Vector(grid.domain);
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
        if((int)opengl_mac_velocity_field->velocity_mode==0) Toggle_Draw();
    }
    else Toggle_Draw();
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
