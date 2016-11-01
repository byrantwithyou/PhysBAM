//#####################################################################
// Copyright 2003-2008, Eran Guendelman, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_SELECTION
//##################################################################### 
#ifndef __OPENGL_SELECTION__
#define __OPENGL_SELECTION__

#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Math_Tools/constants.h>
#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <Geometry/Topology/SIMPLEX_MESH.h>
#include <OpenGL/OpenGL/OPENGL_COLOR.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
namespace PhysBAM{
namespace OPENGL_SELECTION
{
template<class TV> void Draw_Highlighted_Vertex(const TV& position,int id=-1,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV> void Draw_Highlighted_Segment(const TV& x0,const TV& x1,int id=-1,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV> void Draw_Highlighted_Curve(const ARRAY<VECTOR<TV,2> >& X,int id=-1,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV> void Draw_Highlighted_Triangle_Boundary(const TV& x0,const TV& x1,const TV& x2,int id=-1,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV> void Draw_Highlighted_Tetrahedron_Boundary(const TV& x0,const TV& x1,const TV& x2,const TV& x3,int id=-1,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV> void Draw_Highlighted_Quad(const TV& x00,const TV& x11,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV> void Draw_Highlighted_Quad(const TV& node1,const TV& node2,const TV& node3,const TV& node4,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV> void Draw_Highlighted_Box(const TV& x000,const TV& x111,const OPENGL_COLOR& color=OPENGL_PREFERENCES::selection_highlight_color);
template<class TV,int d> void Draw_Vertices_For_Selection(const SIMPLEX_MESH<d>& mesh,const GEOMETRY_PARTICLES<TV>& particles);
}
}

#endif
