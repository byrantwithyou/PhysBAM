//#####################################################################
// Copyright 2007-2008, Jon Gretarsson, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_TRANSFORMED
//#####################################################################
#ifndef __IMPLICIT_OBJECT_TRANSFORMED__
#define __IMPLICIT_OBJECT_TRANSFORMED__

#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Utilities/PHYSBAM_ATTRIBUTE.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class TV,class TRANSFORM> class IMPLICIT_OBJECT_TRANSFORMED_HELPER;

template<class TV> class IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,FRAME<TV> >
{
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
public:

    IMPLICIT_OBJECT_TRANSFORMED_HELPER(const FRAME<TV>* transform_input)
        :transform(transform_input),owns_transform(false)
    {}

    ~IMPLICIT_OBJECT_TRANSFORMED_HELPER()
    {if(owns_transform) delete transform;}

    const FRAME<TV>* transform;
    bool owns_transform;

    bool Scale_Transform(const T scale_input) PHYSBAM_ALWAYS_INLINE
    {return false;}

    TV World_Space_Point(const TV& object_space_point) const
    {return transform->Frame()*object_space_point;}

    TV Object_Space_Point(const TV& world_space_point) const
    {return transform->Frame().Inverse_Times(world_space_point);}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return transform->Frame().r.Rotate(object_space_vector);}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return transform->Frame().r.Inverse_Rotate(world_space_vector);}

    TV World_Space_Unitless_Vector(const TV& object_space_vector) const
    {return World_Space_Vector(object_space_vector);}

    TV Object_Space_Unitless_Vector(const TV& world_space_vector) const
    {return Object_Space_Vector(world_space_vector);}

    const T& World_Space_Length(const T& object_space_length) const PHYSBAM_ALWAYS_INLINE
    {return object_space_length;}

    const T& Object_Space_Length(const T& world_space_length) const PHYSBAM_ALWAYS_INLINE
    {return world_space_length;}

    RANGE<TV> World_Space_Box(const RANGE<TV>& object_space_box) const
    {if(object_space_box==RANGE<TV>::Full_Box() || object_space_box==RANGE<TV>::Empty_Box()) return object_space_box;
    return T_ORIENTED_BOX(object_space_box,transform->Frame()).Axis_Aligned_Bounding_Box();}

    RANGE<TV> Object_Space_Box(const RANGE<TV>& world_space_box) const
    {if(world_space_box==RANGE<TV>::Full_Box() || world_space_box==RANGE<TV>::Empty_Box()) return world_space_box;
    return T_ORIENTED_BOX(world_space_box,transform->Frame().Inverse()).Axis_Aligned_Bounding_Box();}

    SYMMETRIC_MATRIX<T,TV::m> World_Space_Length_Hessian(const SYMMETRIC_MATRIX<T,TV::m>& object_space_hessian) const
    {return SYMMETRIC_MATRIX<T,TV::m>::Conjugate(transform->Frame().r.Rotation_Matrix(),object_space_hessian);}

    static std::string Name_Helper()
    {return LOG::sprintf("IMPLICIT_OBJECT_TRANSFORMED<VECTOR<T,%d>,%s>",TV::dimension,FRAME<TV>::Static_Name().c_str());}

    template<class READ_TRANSFORM> void Read_Transform_Helper(TYPED_ISTREAM& input,const READ_TRANSFORM& transform_input)
    {PHYSBAM_FATAL_ERROR();}

    void Read_Transform_Helper(TYPED_ISTREAM& input,const FRAME<TV>& transform_input)
    {transform=new FRAME<TV>();Read_Binary(input,(FRAME<TV>&)*transform);owns_transform=true;}

    void Read_Transform_Helper(TYPED_ISTREAM& input)
    {Read_Transform_Helper(input,*transform);}

    template<class READ_TRANSFORM> void Write_Transform_Helper(TYPED_OSTREAM& output,const READ_TRANSFORM& transform_input) const
    {PHYSBAM_FATAL_ERROR();}

    void Write_Transform_Helper(TYPED_OSTREAM& output,const FRAME<TV>& transform_input) const
    {Write_Binary(output,transform_input);}

    void Write_Transform_Helper(TYPED_OSTREAM& output) const
    {Write_Transform_Helper(output,*transform);}
};

template<class TV> class IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:

    IMPLICIT_OBJECT_TRANSFORMED_HELPER(const TV& center_input,T scale_input)
        :center(center_input)
    {
        Set_Scale(scale_input);
    }

    ~IMPLICIT_OBJECT_TRANSFORMED_HELPER()
    {}

    TV center;
    T scale,one_over_scale;
    TV center_adjustment;

    void Set_Scale(const T scale_input)
    {scale=scale_input;one_over_scale=1/scale;center_adjustment=center*(1-scale);}

    bool Scale_Transform(const T scale_factor)
    {Set_Scale(scale*scale_factor);return true;}

    TV World_Space_Point(const TV& object_space_point) const
    {return object_space_point*scale+center_adjustment;}

    TV Object_Space_Point(const TV& world_space_point) const
    {return (world_space_point-center_adjustment)*one_over_scale;}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return object_space_vector*scale;}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return world_space_vector*one_over_scale;}

    const TV& World_Space_Unitless_Vector(const TV& object_space_vector) const PHYSBAM_ALWAYS_INLINE
    {return object_space_vector;}

    const TV& Object_Space_Unitless_Vector(const TV& world_space_vector) const PHYSBAM_ALWAYS_INLINE
    {return world_space_vector;}

    T World_Space_Length(const T& object_space_length) const
    {return object_space_length*scale;}

    T Object_Space_Length(const T& world_space_length) const
    {return world_space_length*one_over_scale;}

    RANGE<TV> World_Space_Box(const RANGE<TV>& object_space_box) const
    {if(object_space_box==RANGE<TV>::Full_Box() || object_space_box==RANGE<TV>::Empty_Box()) return object_space_box;
    return RANGE<TV>(World_Space_Point(object_space_box.min_corner),World_Space_Point(object_space_box.max_corner));}

    RANGE<TV> Object_Space_Box(const RANGE<TV>& world_space_box) const
    {if(world_space_box==RANGE<TV>::Full_Box() || world_space_box==RANGE<TV>::Empty_Box()) return world_space_box;
    return RANGE<TV>(Object_Space_Point(world_space_box.min_corner),Object_Space_Point(world_space_box.max_corner));}

    SYMMETRIC_MATRIX<T,TV::m> World_Space_Length_Hessian(const SYMMETRIC_MATRIX<T,TV::m>& object_space_hessian) const
    {return scale*object_space_hessian;}

    static std::string Name_Helper()
    {return LOG::sprintf("IMPLICIT_OBJECT_TRANSFORMED<VECTOR<T,%d>,T>",TV::dimension);}

    void Read_Transform_Helper(TYPED_ISTREAM& input)
    {Read_Binary(input,center,scale,one_over_scale,center_adjustment);}

    void Write_Transform_Helper(TYPED_OSTREAM& output) const
    {Write_Binary(output,center,scale,one_over_scale,center_adjustment);}
};

template<class TV,class TRANSFORM>
class IMPLICIT_OBJECT_TRANSFORMED:public IMPLICIT_OBJECT<TV>,public IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
public:
    typedef int HAS_TYPED_READ_WRITE;
    typedef IMPLICIT_OBJECT<TV> BASE;
    typedef IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM> BASE_HELPER;
    using BASE::box;using BASE_HELPER::Name_Helper;using BASE_HELPER::Object_Space_Box;using BASE_HELPER::World_Space_Length;using BASE_HELPER::World_Space_Box;
    using BASE_HELPER::Scale_Transform;using BASE_HELPER::Object_Space_Point;using BASE_HELPER::World_Space_Unitless_Vector;using BASE_HELPER::Object_Space_Length;
    using BASE_HELPER::World_Space_Point;using BASE_HELPER::World_Space_Length_Hessian;using BASE_HELPER::World_Space_Vector;using BASE_HELPER::Object_Space_Unitless_Vector;
    using BASE_HELPER::Read_Transform_Helper;using BASE_HELPER::Write_Transform_Helper;

    bool owns_implicit_object;
    IMPLICIT_OBJECT<TV>* object_space_implicit_object;

    IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<TV>* object_space_implicit_object_input,bool owns_implicit_object_input,const TRANSFORM* transform_input);
    IMPLICIT_OBJECT_TRANSFORMED(IMPLICIT_OBJECT<TV>* object_space_implicit_object_input,bool owns_implicit_object_input,const TV& center_input,T scale_input);
    virtual ~IMPLICIT_OBJECT_TRANSFORMED();

    static IMPLICIT_OBJECT_TRANSFORMED* Create_Helper(T* input)
    {return new IMPLICIT_OBJECT_TRANSFORMED(0,false,TV(),1);}

    static IMPLICIT_OBJECT_TRANSFORMED* Create_Helper(FRAME<TV>* input)
    {return new IMPLICIT_OBJECT_TRANSFORMED(0,false,0);}

    template<class CREATE_TYPE> static IMPLICIT_OBJECT_TRANSFORMED* Create_Helper(CREATE_TYPE* input)
    {PHYSBAM_FATAL_ERROR();return 0;}

    static IMPLICIT_OBJECT_TRANSFORMED* Create()
    {return Create_Helper((TRANSFORM*)0);}

    RAY<TV> Object_Space_Ray(const RAY<TV>& world_space_ray) const
    {RAY<TV> transformed_ray(Object_Space_Point(world_space_ray.endpoint),Object_Space_Unitless_Vector(world_space_ray.direction));
    transformed_ray.semi_infinite=world_space_ray.semi_infinite;transformed_ray.t_max=Object_Space_Length(world_space_ray.t_max);
    transformed_ray.aggregate_id=world_space_ray.aggregate_id;
    return transformed_ray;}

    T Value(const TV& location) const
    {return World_Space_Length(object_space_implicit_object(Object_Space_Point(location)));}

    T Extended_Value(const TV& location) const
    {return World_Space_Length(object_space_implicit_object->Extended_Phi(Object_Space_Point(location)));}

    static std::string Static_Name()
    {return Name_Helper();}

    void Read(TYPED_ISTREAM& input) override
    {Read_Transform_Helper(input);std::string name;Read_Binary(input,name);
    object_space_implicit_object=dynamic_cast<IMPLICIT_OBJECT<TV>*>(STRUCTURE<TV>::Create_From_Name(name));
    object_space_implicit_object->Read(input);}

    void Write(TYPED_OSTREAM& output) const override
    {Write_Transform_Helper(output);Write_Binary(output,object_space_implicit_object->Name(),*object_space_implicit_object);}

    bool Intersection(RAY<TV>& ray,const T thickness) const override;
    TV Closest_Point_On_Boundary(const TV& location,const T tolerance=0,const int max_iterations=1,T* distance=0) const override;
    virtual std::string Name() const override;
    RANGE<TV>& Box() override;
    void Update_Box() override;
    void Update_Minimum_Cell_Size(const int maximum_depth=0) override;
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override;
    T operator()(const TV& location) const override;
    T Signed_Distance(const TV& location) const override;
    T Extended_Phi(const TV& location) const override;
    T Phi_Secondary(const TV& location) const override;
    TV Normal(const TV& location,const int aggregate=-1) const override;
    TV Extended_Normal(const TV& location,const int aggregate=-1) const override;
    void Compute_Normals() override;
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) override;
    void Rescale(const T scaling_factor) override;
    void Inflate(const T inflation_distance) override;
    bool Inside(const TV& location,const T thickness_over_two=0) const override;
    bool Lazy_Inside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override;
    bool Outside(const TV& location,const T thickness_over_two=0) const override;
    bool Lazy_Outside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override;
    T Min_Phi() const override;
    TV Velocity(const TV& location) const override;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override;
    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const override
    {
        VECTOR<T,d-1> curvatures=object_space_implicit_object->Principal_Curvatures(Object_Space_Point(X));
        for(int i=0;i<d-1;i++) curvatures(i)=Object_Space_Length(curvatures(i)); // Note: Curvatures transform "backwards"
        return curvatures;
    }
    T Integration_Step(const T phi) const override;
    T Minimum_Cell_Size() const override;
//#####################################################################
};
}
#endif
