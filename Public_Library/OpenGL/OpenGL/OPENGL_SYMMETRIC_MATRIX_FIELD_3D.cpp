//#####################################################################
// Copyright 2004-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <OpenGL/OpenGL/OPENGL_SYMMETRIC_MATRIX_FIELD_3D.h>
#include <OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T>::
Display() const
{
    glPushAttrib(GL_LIGHTING_BIT|GL_TEXTURE_BIT|GL_LINE_BIT);
    glLineWidth(1);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
    OpenGL_Begin(GL_LINES);
    for(int i=0;i<entries.m;i++){
        VECTOR<T,3> node=entries(i).x;MATRIX<T,3> line=size*entries(i).y;VECTOR<bool,3> p=entries(i).z;
        (p.x?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(0),node+line.Column(0));
        (p.y?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(1),node+line.Column(1));
        (p.z?positive_color:negative_color).Send_To_GL_Pipeline();
        OpenGL_Line(node-line.Column(2),node+line.Column(2));}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T>::
Update()
{
    if(field.domain.Empty()){entries.Clean_Memory();return;}

    int m_start=0,m_end=grid.counts.x,n_start=0,n_end=grid.counts.y,mn_start=0,mn_end=grid.counts.z;
    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;
    if(slice && slice->mode!=OPENGL_SLICE::NO_SLICE){
        VECTOR<int,3> domain_start(m_start,n_start,mn_start),domain_end(m_end,n_end,mn_end);
        if((slice->mode == OPENGL_SLICE::CELL_SLICE && (!grid.Is_MAC_Grid() || slice->index < domain_start[slice->axis] || slice->index >= domain_end[slice->axis])) ||
           (slice->mode == OPENGL_SLICE::NODE_SLICE && (grid.Is_MAC_Grid() || slice->index < domain_start[slice->axis] || slice->index >= domain_end[slice->axis]))) return;
        switch(slice->axis){
            case 0:m_start=slice->index;m_end=m_start+1;break;
            case 1:n_start=slice->index;m_end=m_start+1;break;
            case 2:mn_start=slice->index;m_end=m_start+1;break;}}

    int count=0;
    for(int i=m_start;i<m_end;i++)for(int j=m_start;j<n_end;j++)for(int ij=mn_start;ij<mn_end;ij++) if(field(i,j,ij).Frobenius_Norm_Squared())count++;
    entries.Resize(count);count=0;
    DIAGONAL_MATRIX<T,3> D;MATRIX<T,3> U;
    for(int i=m_start;i<m_end;i++)for(int j=m_start;j<n_end;j++)for(int ij=mn_start;ij<mn_end;ij++) if(field(i,j,ij).Frobenius_Norm_Squared()){
        field(i,j,ij).Solve_Eigenproblem(D,U);entries(count++)=TRIPLE<VECTOR<T,3>,MATRIX<T,3>,VECTOR<bool,3> >(grid.X(VECTOR<int,3>(i,j,ij)),U*D,VECTOR<bool,3>(D.x.x>0,D.x.y>0,D.x.z>0));}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_SYMMETRIC_MATRIX_FIELD_3D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<T,3> >(grid.domain);
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_SYMMETRIC_MATRIX_FIELD_3D<float>;
template class OPENGL_SYMMETRIC_MATRIX_FIELD_3D<double>;
}
