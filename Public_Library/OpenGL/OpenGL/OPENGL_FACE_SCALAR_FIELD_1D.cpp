//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Math_Tools/RANGE.h>
#include <OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_1D.h>
#include <OpenGL/OpenGL/OPENGL_GRID_1D.h>
#include <OpenGL/OpenGL/OPENGL_WORLD.h>
namespace PhysBAM{
//#####################################################################
// OPENGL_FACE_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
OPENGL_FACE_SCALAR_FIELD_1D(STREAM_TYPE stream_type,const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<1> > &face_values_input,OPENGL_COLOR point_color_input,OPENGL_COLOR line_color_input)
:OPENGL_OBJECT<T>(stream_type),grid(grid_input),face_values(face_values_input),point_color(point_color_input),line_color(line_color_input),scale(1)
{
}
//#####################################################################
// ~OPENGL_FACE_SCALAR_FIELD_1D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
~OPENGL_FACE_SCALAR_FIELD_1D()
{
}
//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
Display() const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    line_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_LINE_STRIP);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        OpenGL_Vertex(VECTOR<T,2>(iterator.Location().x,scale*face_values(iterator.Full_Index())));}
    OpenGL_End();
    glColor3f(0,1,1);
    glPointSize(3.0);
    point_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_POINTS);
    for(FACE_ITERATOR<TV> iterator(grid);iterator.Valid();iterator.Next()){
        OpenGL_Vertex(VECTOR<T,2>(iterator.Location().x,scale*face_values(iterator.Full_Index())));}
    OpenGL_End();
    glPopAttrib();
}
template<class T> void
Display_Bool_Helper(const OPENGL_FACE_SCALAR_FIELD_1D<T,bool>& self)
{
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glPointSize(8.0);
    self.point_color.Send_To_GL_Pipeline();
    OpenGL_Begin(GL_POINTS);
    for(FACE_ITERATOR<VECTOR<T,1> > iterator(self.grid);iterator.Valid();iterator.Next())
        if(self.face_values(iterator.Full_Index())){
            OpenGL_Vertex(VECTOR<T,2>(iterator.Location().x,(T)0));}
    OpenGL_End();
    glPopAttrib();
}
template<> void OPENGL_FACE_SCALAR_FIELD_1D<float,bool>::
Display() const
{
    Display_Bool_Helper(*this);
}
template<> void OPENGL_FACE_SCALAR_FIELD_1D<double,bool>::
Display() const
{
    Display_Bool_Helper(*this);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<T,3> > OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_1D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION<T>* selection) const
{
    // TODO: this should also interpolate to particles
    if(selection && selection->type==OPENGL_SELECTION<T>::GRID_CELL_1D && grid.Is_MAC_Grid()){
        VECTOR<int,1> index=((OPENGL_SELECTION_GRID_CELL_1D<T>*)selection)->index;
        FACE_INDEX<TV::m> ix(0,index);
        T2 left=face_values(ix);
        ix.index.x++;
        T2 right=face_values(ix);
        output_stream<<"    left = "<<left<<",right = "<<right<<std::endl;}
}
//#####################################################################
template class OPENGL_FACE_SCALAR_FIELD_1D<float,int>;
template class OPENGL_FACE_SCALAR_FIELD_1D<float,bool>;
template class OPENGL_FACE_SCALAR_FIELD_1D<float,float>;
template class OPENGL_FACE_SCALAR_FIELD_1D<double,int>;
template class OPENGL_FACE_SCALAR_FIELD_1D<double,bool>;
template class OPENGL_FACE_SCALAR_FIELD_1D<double,double>;
}
