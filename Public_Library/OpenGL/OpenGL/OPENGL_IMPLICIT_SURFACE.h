//#####################################################################
// Copyright 2002, Robert Bridson, Eilene Hao, Geoffrey Irving, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_IMPLICIT_SURFACE
//##################################################################### 
//
//#####################################################################
// Bridson - November 12, 2002
// Molino - October 27, 2002
// Irving - August 7, 2003
// Hao - August 22, 2003
//#####################################################################
#ifndef __OPENGL_IMPLICIT_SURFACE__
#define __OPENGL_IMPLICIT_SURFACE__

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <OpenGL/OpenGL/OPENGL_MATERIAL.h>
#include <OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <cstring> // for gcc 2.96 for memset
namespace PhysBAM{

template<class T>
class OPENGL_IMPLICIT_SURFACE:public OPENGL_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using OPENGL_OBJECT<T>::Send_Transform_To_GL_Pipeline;
    IMPLICIT_OBJECT<TV>& surface;
    bool two_sided;
    OPENGL_MATERIAL front_material;
    OPENGL_MATERIAL back_material;
    int nx,ny,nz; // resolution of grid in eye space
    mutable float *slice1,*slice2;
    mutable unsigned int *up_to_date1,*up_to_date2,up_to_date_counter;

    OPENGL_IMPLICIT_SURFACE(IMPLICIT_OBJECT<TV>& surface_input,const OPENGL_MATERIAL& material_input=OPENGL_MATERIAL::Plastic(OPENGL_COLOR(.6f,.5f,1.f)));
    OPENGL_IMPLICIT_SURFACE(IMPLICIT_OBJECT<TV>& surface_input,const OPENGL_MATERIAL &front_material_input,const OPENGL_MATERIAL& back_material_input);
    virtual ~OPENGL_IMPLICIT_SURFACE();

    virtual RANGE<VECTOR<T,3> > Bounding_Box() const PHYSBAM_OVERRIDE
    {return surface.box;}

    void Set_Two_Sided(bool two_sided_input=true)
    {two_sided=two_sided_input;}

    void Set_Front_Material(const OPENGL_MATERIAL &front_material_input)
    {front_material=front_material_input;}

    void Set_Back_Material(const OPENGL_MATERIAL &back_material_input)
    {back_material=back_material_input;}

//#####################################################################
    void Set_Resolution(const int nx_input,const int ny_input,const int nz_input);
    void Initialize_Slices();
    float inline Phi(VECTOR<T,3> x) const;
    void Display() const PHYSBAM_OVERRIDE;
    void Display_Tetrahedron(const VECTOR<T,3>& x0,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3, 
        const float phi1,const float phi2,const float phi3,const float phi4,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices,ARRAY<GLfloat>& normals) const;
    void Display_Brick(const VECTOR<T,3>& x0,float p1,const VECTOR<T,3>& x1,float p2,const VECTOR<T,3>& x2,float p3,const VECTOR<T,3>& x3,
        float p4,const VECTOR<T,3>& x4,float p5,const VECTOR<T,3>& x5,float p6,const VECTOR<T,3>& x6,float p7,const VECTOR<T,3>& x7,
        float p8,const int parity,ARRAY<typename OPENGL_POLICY<T>::T_GL>& vertices,ARRAY<GLfloat>& normals) const;
//#####################################################################
};
}
#endif
