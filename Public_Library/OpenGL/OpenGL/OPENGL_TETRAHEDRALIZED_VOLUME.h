//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_TETRAHEDRALIZED_VOLUME
//#####################################################################
#ifndef __OPENGL_TETRAHEDRALIZED_VOLUME__
#define __OPENGL_TETRAHEDRALIZED_VOLUME__

#include <Core/Vectors/VECTOR_3D.h>
#include <Geometry/Topology/TETRAHEDRON_MESH.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T>
class OPENGL_TETRAHEDRALIZED_VOLUME:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;using OPENGL_OBJECT<T>::World_Space_Box;
    OPENGL_MATERIAL material;
    OPENGL_MATERIAL inverted_material;
    bool use_inverted_material;
    TETRAHEDRON_MESH* mesh;
    GEOMETRY_PARTICLES<VECTOR<T,3> >* particles;
    int current_tetrahedron;
    int current_node;
    int current_boundary_triangle;
    int minimum_valence;
    ARRAY<int> subset;
    ARRAY<int> subset_triangles;
    ARRAY<int> subset_particles;
    bool boundary_only;
    ARRAY<T> spectrum;
    bool draw_subsets;
    int cutaway_mode;
    T cutaway_fraction;
    TETRAHEDRON_MESH cutaway_mesh;
    ARRAY<OPENGL_COLOR>* color_map;
protected:
    bool smooth_normals;
    ARRAY<VECTOR<T,3> >* vertex_normals;
    int selected_vertex;
    int selected_tet;
public:

    OPENGL_TETRAHEDRALIZED_VOLUME(TETRAHEDRON_MESH* mesh_input,GEOMETRY_PARTICLES<VECTOR<T,3> >* particles_input,const OPENGL_MATERIAL& material_input,
        const OPENGL_MATERIAL& inverted_material_input,bool initialize=true,ARRAY<OPENGL_COLOR>* color_map_input=0);
    ~OPENGL_TETRAHEDRALIZED_VOLUME();

    void Display() const override;
    virtual RANGE<TV> Bounding_Box() const override;

    virtual int Get_Selection_Priority(ARRAY_VIEW<GLuint> indices) override;
    virtual bool Set_Selection(ARRAY_VIEW<GLuint> indices,int modifiers) override;
    virtual RANGE<TV> Selection_Bounding_Box() const override;
    void Clear_Selection() override;

    void Set_Boundary_Only(bool boundary_only_input)
    {boundary_only=boundary_only_input;}

    void Draw_Boundary_Only()
    {boundary_only=true;}

    void Draw_Interior()
    {boundary_only=false;}

    void Set_Color_Map(ARRAY<OPENGL_COLOR>* color_map_input)
    {color_map=color_map_input;}

    bool Toggle_Boundary_Only()
    {boundary_only=!boundary_only;return boundary_only;}

    bool Is_Smooth_Normals()
    {return smooth_normals;}

    void Cycle_Cutaway_Mode()
    {cutaway_mode=(cutaway_mode+1)%8;if(cutaway_mode) Update_Cutaway_Plane();}

    void Decrease_Cutaway_Fraction()
    {cutaway_fraction-=(T).05;if(cutaway_fraction<0)cutaway_fraction+=1;if(cutaway_mode) Update_Cutaway_Plane();}

    void Increase_Cutaway_Fraction()
    {cutaway_fraction+=(T).05;if(cutaway_fraction>1)cutaway_fraction-=1;if(cutaway_mode) Update_Cutaway_Plane();}

    bool Toggle_Differentiate_Inverted()
    {use_inverted_material=!use_inverted_material;return use_inverted_material;}

//#####################################################################
    void Draw_Boundary_Triangles(const TETRAHEDRON_MESH& tetrahedron_mesh) const;
    void Draw_Current_Tetrahedron() const;
    void Draw_Subset() const;
    void Draw_Subset_Triangles() const;
    void Draw_Subset_Particles() const;
    void Draw_In_Color_From_Spectrum() const;
    void Draw_In_Color_From_Color_Map() const;
    void Draw_Wireframe_Mesh(const TETRAHEDRON_MESH& tetrahedron_mesh) const;
    void Highlight_Boundary_Nodes_Of_Current_Tetrahedron() const;
    void Highlight_Current_Node() const;
    void Highlight_Nodes_Of_Minimum_Valence() const;
    void Highlight_Boundary_Normal_Vectors_Of_Current_Tetrahedron() const;
    void Highlight_Current_Boundary_Triangle() const;
    void Turn_Smooth_Shading_On() override;
    void Turn_Smooth_Shading_Off() override;
    void Display_Subset();
    void Update_Cutaway_Plane();
    void Print_Selection_Info(std::ostream &output_stream) const override;
    void Print_Selection_Info(std::ostream &output_stream,MATRIX<T,4>* transform) const;
    void Initialize_Vertex_Normals();
protected:
    void Draw_Vertices_For_Selection() const;
    void Draw_Tetrahedra_For_Selection() const;
    int Find_Shortest_Spring(const VECTOR<int,4>& nodes,T& distance,TV& minimum_normal,TV& weights) const;
//#####################################################################
};
}
#endif
