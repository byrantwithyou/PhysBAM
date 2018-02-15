//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Arrays_Nd/ARRAYS_ND_VIEW.h>
#include <Core/Math_Tools/RANGE.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_3D.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
namespace PhysBAM{
//#####################################################################
// OPENGL_FACE_SCALAR_FIELD_3D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
OPENGL_FACE_SCALAR_FIELD_3D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<3> > &face_values_input,OPENGL_COLOR_MAP<T2> *color_map_input)
    :grid(grid_input),face_values(face_values_input),
    color_map(color_map_input),opengl_points(*new ARRAY<TV>)
{
    PHYSBAM_ASSERT(color_map);
}
//#####################################################################
// ~OPENGL_FACE_SCALAR_FIELD_3D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
~OPENGL_FACE_SCALAR_FIELD_3D()
{
    delete &opengl_points.points;
    delete color_map;
}
//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Display() const
{
    if(face_values.Component(0).domain.Empty()) return;
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    opengl_points.Display();
    glPopMatrix();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Update
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Update()
{
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;
    opengl_points.points.Remove_All();
    opengl_points.Store_Point_Colors(true);
    for(FACE_ITERATOR<TV> it(slice?slice->Get_Slice_Grid():grid,face_values.Number_Of_Ghost_Cells());it.Valid();it.Next()){
        opengl_points.points.Append(it.Location());
        opengl_points.point_colors->Append(color_map->Lookup(face_values(it.Full_Index())));}
}
//#####################################################################
// Function Bool_Update_Helper
//#####################################################################
template<class T> static void 
Bool_Update_Helper(OPENGL_FACE_SCALAR_FIELD_3D<T,bool>& self)
{
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)self.slice;
    self.opengl_points.points.Remove_All();
    for(FACE_ITERATOR<VECTOR<T,3> > it(slice?slice->Get_Slice_Grid():self.grid,self.face_values.Number_Of_Ghost_Cells());it.Valid();it.Next())
        if(self.face_values(it.Full_Index()))
            self.opengl_points.points.Append(it.Location());
}
//#####################################################################
// Update
//#####################################################################
template<> void OPENGL_FACE_SCALAR_FIELD_3D<float,bool>::
Update()
{
    Bool_Update_Helper(*this);
}
//#####################################################################
// Update
//#####################################################################
template<> void OPENGL_FACE_SCALAR_FIELD_3D<double,bool>::
Update()
{
    Bool_Update_Helper(*this);
}
//#####################################################################
// Slice_Has_Changed
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Slice_Has_Changed()
{
    Update();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& output_stream) const
{
    // TODO: this should also interpolate to particles
    if(selected_index.x>=0 && grid.Is_MAC_Grid()){
        FACE_INDEX<TV::m> ix(0,selected_index),iy(1,selected_index),iz(2,selected_index);
        T2 left=face_values(ix);
        T2 bottom=face_values(iy);
        T2 back=face_values(iz);
        ix.index.x++;
        iy.index.y++;
        iz.index.z++;
        T2 right=face_values(ix);
        T2 top=face_values(iy);
        T2 front=face_values(iz);
        output_stream<<"    left = "<<left<<",right = "<<right<<std::endl;
        output_stream<<"    bottom = "<<bottom<<",top = "<<top<<std::endl;
        output_stream<<"    back = "<<back<<",front = "<<front<<std::endl;}
}
//#####################################################################
template class OPENGL_FACE_SCALAR_FIELD_3D<float,int>;
template class OPENGL_FACE_SCALAR_FIELD_3D<float,bool>;
template class OPENGL_FACE_SCALAR_FIELD_3D<float,float>;
template class OPENGL_FACE_SCALAR_FIELD_3D<double,int>;
template class OPENGL_FACE_SCALAR_FIELD_3D<double,bool>;
template class OPENGL_FACE_SCALAR_FIELD_3D<double,double>;
}
