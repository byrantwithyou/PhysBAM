//#####################################################################
// Copyright 2007, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_COMBINED
//#####################################################################
#ifndef __IMPLICIT_OBJECT_COMBINED__
#define __IMPLICIT_OBJECT_COMBINED__

#include <Core/Matrices/MATRIX_3X3.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT_COMBINED:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    enum WORKAROUND {d=TV::m};
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    IMPLICIT_OBJECT<TV> *implicit_object1,*implicit_object2;
    bool owns_implicit_object1,owns_implicit_object2;
    T alpha;
    
    IMPLICIT_OBJECT_COMBINED(IMPLICIT_OBJECT<TV>* implicit_object1_input,bool owns_implicit_object1_input,IMPLICIT_OBJECT<TV>* implicit_object2_input, bool owns_implicit_object2_input)
        :implicit_object1(implicit_object1_input),implicit_object2(implicit_object2_input),owns_implicit_object1(owns_implicit_object1_input),owns_implicit_object2(owns_implicit_object2_input)
    {}

    virtual ~IMPLICIT_OBJECT_COMBINED()
    {if(owns_implicit_object1) delete implicit_object1;if(owns_implicit_object2) delete implicit_object2;}

    void Set_Weights(T alpha_input)
    {alpha=alpha_input;}

    void Update_Box() override
    {implicit_object1->Update_Box();implicit_object2->Update_Box();box=implicit_object1->box;
    box.Enlarge_Nonempty_Box_To_Include_Points(implicit_object2->box.min_corner,implicit_object2->box.max_corner);}

    void Update_Minimum_Cell_Size(const int maximum_depth=0) override
    {implicit_object1->Update_Minimum_Cell_Size(maximum_depth);implicit_object2->Update_Minimum_Cell_Size(maximum_depth);}

    // TODO: box is not used.  Is it still needed?
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override
    {return min(implicit_object1->Minimum_Cell_Size_Within_Box(box),implicit_object2->Minimum_Cell_Size_Within_Box(box));} 

    T operator()(const TV& location) const override
    {return (1-alpha)*(*implicit_object1)(location)+alpha*(*implicit_object2)(location);}

    T Signed_Distance(const TV& location) const override
    {return (*this)(location);} // to make this class compatible with the geometry classes

    T Extended_Phi(const TV& location) const override
    {return Extended_Value(location);}

    T Phi_Secondary(const TV& location) const override
    {return (1-alpha)*implicit_object1->Phi_Secondary(location)+alpha*implicit_object2->Phi_Secondary(location);}

    TV Normal(const TV& location,const int aggregate) const override
    {return (1-alpha)*implicit_object1->Normal(location,aggregate)+alpha*implicit_object2->Normal(location,aggregate);}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const override
    {return (1-alpha)*implicit_object1->Extended_Normal(location,aggregate)+alpha*implicit_object2->Extended_Normal(location,aggregate);}

    void Compute_Normals() override
    {implicit_object1->Compute_Normals();implicit_object2->Compute_Normals();}

    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) override
    {implicit_object1->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);implicit_object2->Compute_Cell_Minimum_And_Maximum(recompute_if_exists);}
    
    void Rescale(const T scaling_factor) override
    {implicit_object1->Rescale(scaling_factor);implicit_object2->Rescale(scaling_factor);}

    void Inflate(const T inflation_distance) override
    {implicit_object1->Inflate(inflation_distance);implicit_object2->Inflate(inflation_distance);}

    bool Lazy_Inside(const TV& location,const T contour_value=0) const override
    {return box.Inside(location,-contour_value) && (*this)(location)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& location,T& phi,const T contour_value=0) const override
    {if(box.Inside(location,-contour_value)){phi=(*this)(location);if(phi<=contour_value) return true;} return false;}

    bool Lazy_Inside_Extended_Levelset(const TV& location,const T contour_value=0) const override
    {return box.Inside(location,-contour_value) && Extended_Value(location)<=contour_value;}

    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override
    {if(box.Inside(location,-contour_value)){phi_value=Extended_Phi(location);if(phi_value<=contour_value) return true;} return false;}

    bool Lazy_Outside(const TV& location,const T contour_value=0) const override
    {return !Lazy_Inside(location,contour_value);}

    bool Lazy_Outside_Extended_Levelset(const TV& location,const T contour_value=0) const override
    {return !Lazy_Inside_Extended_Levelset(location,contour_value);}

    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override
    {return !Lazy_Inside_Extended_Levelset_And_Value(location,phi_value,contour_value);}

    TV Velocity(const TV& location) const override
    {return (1-alpha)*implicit_object1->Velocity(location)+alpha*implicit_object2->Velocity(location);}

    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override
    {return (1-alpha)*implicit_object1->Hessian(X)+alpha*implicit_object2->Hessian(X);}

    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const override
    {return (1-alpha)*implicit_object1->Principal_Curvatures(X)+alpha*implicit_object2->Principal_Curvatures(X);}

    T Integration_Step(const T phi) const override
    {return (1-alpha)*implicit_object1->Integration_Step(phi)+alpha*implicit_object2->Integration_Step(phi);}

    T Minimum_Cell_Size() const override
    {return min(implicit_object1->Minimum_Cell_Size(),implicit_object2->Minimum_Cell_Size());}

    T Value(const TV& location) const
    {return (1-alpha)*(*implicit_object1)(location)+alpha*(*implicit_object2)(location);}

    T Extended_Value(const TV& location) const
    {return (1-alpha)*implicit_object1->Extended_Phi(location)+alpha*implicit_object1->Extended_Phi(location);}

    bool Intersection(RAY<TV>& ray,const T thickness) const override
    {PHYSBAM_NOT_IMPLEMENTED();}

    virtual void Read(TYPED_ISTREAM input) override {PHYSBAM_FATAL_ERROR();}
    virtual void Write(TYPED_OSTREAM output) const override {PHYSBAM_FATAL_ERROR();}

//#####################################################################
};
}
#endif

