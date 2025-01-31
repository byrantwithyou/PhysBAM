//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Sergey Koltakov, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_PLANE
//#####################################################################
#ifndef __RENDERING_PLANE__
#define __RENDERING_PLANE__

#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Intersections/RAY_PLANE_INTERSECTION.h>
#include <Geometry/Topology/TRIANGLE_MESH.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_PLANE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::World_Space_Vector;using RENDERING_OBJECT<T>::Object_Space_Point;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection;using RENDERING_OBJECT<T>::World_Space_Point;
    using RENDERING_OBJECT<T>::Object_Space_Ray;

    PLANE<T> plane;
    TV texture_vector1,texture_vector2;

    RENDERING_PLANE()
    {}

    RENDERING_PLANE(const TV& normal_input,const TV& x1_input)
        :plane(normal_input,x1_input)
    {}
    
    RENDERING_PLANE(const TV& x1_input,const TV& x2_input,const TV& x3_input)
        :plane(x1_input,x2_input,x3_input)
    {}

    virtual ~RENDERING_PLANE()
    {}

    bool Intersection(RAY<TV>& ray) const override
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,plane,small_number)){ray.semi_infinite=false;
        ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const override
    {object_tangent=texture_vector1;object_bitangent=texture_vector2;}

    void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const override
    {world_tangent=World_Space_Vector(texture_vector1);world_bitangent=World_Space_Vector(texture_vector2);}

    TV Normal(const TV& location,const int aggregate=0) const override
    {return World_Space_Vector(plane.Normal());}

    bool Inside(const TV& location) const override
    {return plane.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const override
    {return plane.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const override
    {return plane.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const override
    {return World_Space_Point(plane.Surface(Object_Space_Point(location)));}

    T Signed_Distance(const TV& location) const override
    {return plane.Signed_Distance(Object_Space_Point(location));}

    void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const override
    {TV p=object_space_point-plane.x0;
    s=TV::Dot_Product(p,texture_vector1);t=TV::Dot_Product(p,texture_vector2);}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const override
    {TRIANGULATED_SURFACE<T>* surface;
    TV u_vector=texture_vector1,v_vector=texture_vector2;
    GEOMETRY_PARTICLES<TV>* particles=new GEOMETRY_PARTICLES<TV>();
    int vertex_1=particles->Add_Element(),vertex_2=particles->Add_Element(),vertex_3=particles->Add_Element(),vertex_4=particles->Add_Element();
    particles->X(vertex_1)=(plane.x0);particles->X(vertex_2)=(plane.x0+u_vector);
    particles->X(vertex_3)=(plane.x0+u_vector+v_vector);particles->X(vertex_4)=(plane.x0+v_vector);
    ARRAY<VECTOR<int,3> > triangles(2);triangles(0).Set(0,1,3);triangles(1).Set(3,1,2);
    TRIANGLE_MESH* mesh=new TRIANGLE_MESH(particles->Size(),triangles);
    surface=new TRIANGULATED_SURFACE<T>(*mesh,*particles);
    surface->Update_Triangle_List();surface->Update_Vertex_Normals();return surface;};
    
//#####################################################################
};   
}
#endif

