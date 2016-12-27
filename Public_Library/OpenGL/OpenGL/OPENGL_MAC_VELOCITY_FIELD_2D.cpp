//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <OpenGL/OpenGL/OPENGL_MAC_VELOCITY_FIELD_2D.h>
#include <OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
using namespace PhysBAM;
using namespace std;
//#####################################################################
// Constructor
//#####################################################################
template<class T> OPENGL_MAC_VELOCITY_FIELD_2D<T>::
OPENGL_MAC_VELOCITY_FIELD_2D(STREAM_TYPE stream_type,const GRID<TV> &grid_input)
    :OPENGL_VECTOR_FIELD_2D<ARRAY<TV> >(stream_type,vector_field,vector_locations),
    grid(grid_input),u(face_velocities.Component(0)),v(face_velocities.Component(1)),active_cells(0),active_faces(0)
{
    PHYSBAM_ASSERT(grid.Is_MAC_Grid());
    Set_Velocity_Mode(CELL_CENTERED);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> OPENGL_MAC_VELOCITY_FIELD_2D<T>::
~OPENGL_MAC_VELOCITY_FIELD_2D()
{
    delete active_cells;
    delete active_faces;
}
//#####################################################################
// Function Set_Velocity_Mode
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Set_Velocity_Mode(VELOCITY_MODE velocity_mode_input)
{
    velocity_mode = velocity_mode_input;
    Update();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Update()
{
    PHYSBAM_ASSERT(grid.Is_MAC_Grid());

    vector_field.Resize(0);
    vector_locations.Resize(0);
    TV_INT ucounts=u.Size(),vcounts=v.Size();
    if(ucounts.x-1!=vcounts.x || ucounts.y!=vcounts.y-1) return;

    int idx=0;
    if(velocity_mode == FACE_CENTERED){
        vector_field.Resize(ucounts.x*ucounts.y + vcounts.x*vcounts.y);
        vector_locations.Resize(ucounts.x*ucounts.y + vcounts.x*vcounts.y);
        for(int i=u.domain.min_corner.x;i<u.domain.max_corner.x;i++) for(int j=u.domain.min_corner.y;j<u.domain.max_corner.y;j++) if(!active_faces||(active_faces->Component(0))(i,j)){
            vector_field(idx)=VECTOR<T,2>(u(i,j),0);vector_locations(idx)=grid.Face(FACE_INDEX<TV::m>(0,TV_INT(i,j)));idx++;}
        for(int i=v.domain.min_corner.x;i<v.domain.max_corner.x;i++) for(int j=v.domain.min_corner.y;j<v.domain.max_corner.y;j++) if(!active_faces||(active_faces->Component(1))(i,j)){
            vector_field(idx)=VECTOR<T,2>(0,v(i,j));vector_locations(idx)=grid.Face(FACE_INDEX<TV::m>(1,TV_INT(i,j)));idx++;}
        vector_field.Resize(idx);
        vector_locations.Resize(idx);}
    else{
        int number_of_cells=(ucounts.x-1)*(vcounts.y-1);
        vector_field.Resize(number_of_cells);
        vector_locations.Resize(number_of_cells);
        for(int i=u.domain.min_corner.x;i<u.domain.max_corner.x-1;i++) for(int j=v.domain.min_corner.y;j<v.domain.max_corner.y-1;j++) if(!active_cells||(*active_cells)(i,j)){
            vector_field(idx)=VECTOR<T,2>((T).5*(u(i,j)+u(i+1,j)),(T).5*(v(i,j)+v(i,j+1)));vector_locations(idx)=grid.Center(TV_INT(i,j));idx++;}
        vector_field.Resize(idx);
        vector_locations.Resize(idx);}
}
//#####################################################################
// Function Print_Cell_Selection_Info
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Print_Cell_Selection_Info(std::ostream& stream,const TV_INT& cell) const
{
    if(grid.Is_MAC_Grid()){
        T u_left,u_right,v_bottom,v_top;
        u_left=u(cell.x,cell.y);
        u_right=u(cell.x+1,cell.y);
        v_bottom=v(cell.x,cell.y);
        v_top=v(cell.x,cell.y+1);
        TV center_velocity(0.5*(u_left+u_right),0.5*(v_bottom+v_top));
        stream<<"    u left = "<<u_left<<",right = "<<u_right<<" avg="<<center_velocity.x<<std::endl;
        stream<<"    v bottom = "<<v_bottom<<",top = "<<v_top<<" avg="<<center_velocity.y<<std::endl;
        T ux=(u_right-u_left)*grid.one_over_dX.x,vy=(v_top-v_bottom)*grid.one_over_dX.y;
        stream<<"    divergence = "<<ux+vy<<" (ux="<<ux<<", vy="<<vy<<")"<<std::endl;}
    // if(selection && selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_2D){
    //     OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *particle_selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)selection;
    //     LINEAR_INTERPOLATION_UNIFORM<TV,T> interpolation;
    //     VECTOR<T,2> interp(interpolation.Clamped_To_Array(grid.Get_Face_Grid(0),u,particle_selection->location),interpolation.Clamped_To_Array(grid.Get_Face_Grid(1),v,particle_selection->location));
    //     stream<<"    @ particle = "<<interp<<std::endl;}
}
//#####################################################################
// Function Print_Node_Selection_Info
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Print_Node_Selection_Info(std::ostream& stream,const TV_INT& node) const
{
}
//#####################################################################
// Convenience functions
//#####################################################################
template<class T> void OPENGL_MAC_VELOCITY_FIELD_2D<T>::
Toggle_Velocity_Mode()
{
    VELOCITY_MODE new_velocity_mode = (VELOCITY_MODE)(((int)velocity_mode+1)%2);
    Set_Velocity_Mode(new_velocity_mode);
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_MAC_VELOCITY_FIELD_2D<float>;
template class OPENGL_MAC_VELOCITY_FIELD_2D<double>;
}
