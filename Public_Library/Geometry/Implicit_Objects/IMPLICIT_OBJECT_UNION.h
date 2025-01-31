//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_UNION
//#####################################################################
#ifndef __IMPLICIT_OBJECT_UNION__
#define __IMPLICIT_OBJECT_UNION__

#include <Core/Math_Tools/RANGE.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV>
class IMPLICIT_OBJECT_UNION:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<T,d-1> T_CURVATURES;
public:
    typedef int HAS_TYPED_READ_WRITE;
    typedef TV VECTOR_T;
    
    using IMPLICIT_OBJECT<TV>::box;
    using IMPLICIT_OBJECT<TV>::use_secondary_interpolation;

    ARRAY<IMPLICIT_OBJECT<TV>*> io;
    ARRAY<bool> owns_io;

    template<class...Args>
    IMPLICIT_OBJECT_UNION(Args&&... args)
    {
        Add_Implicit_Object(args...);
    }
    virtual ~IMPLICIT_OBJECT_UNION();

    template<class...Args>
    void Add_Implicit_Object(IMPLICIT_OBJECT<TV>* o,Args&&... args)
    {io.Append(o);owns_io.Append(true);Add_Implicit_Object(args...);}

    void Add_Implicit_Object(){}

    void Use_Secondary_Interpolation(const bool use_secondary_interpolation_input=true)
    {use_secondary_interpolation=use_secondary_interpolation_input;}

//#####################################################################
    void Update_Box() override;
    void Update_Minimum_Cell_Size(const int maximum_depth=0) override;
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override;
    T operator()(const TV& location) const override;
    T Extended_Phi(const TV& location) const override;
    T Phi_Secondary(const TV& location) const override;
    TV Normal(const TV& location,const int aggregate=-1) const override;
    TV Extended_Normal(const TV& location,const int aggregate=-1) const override;
    void Compute_Normals() override;
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) override;
    void Rescale(const T scaling_factor) override;
    void Translate(const TV& translation) override;
    void Inflate(const T inflation_distance) override;
    bool Inside(const TV& location,const T thickness_over_two=0) const override;
    bool Outside(const TV& location,const T thickness_over_two=0) const override;
    bool Boundary(const TV& location,const T thickness_over_two) const override;
    bool Lazy_Inside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override;
    bool Lazy_Outside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override;
    bool Intersection(RAY<TV>& ray,const T thickness=0) const override;
    TV Closest_Point_On_Boundary(const TV& location,const T tolerance=0,const int max_iterations=1,T* distance=0) const override;
    TV Velocity(const TV& location) const override;
    // the following only exist in 3d
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override;
    T_CURVATURES Principal_Curvatures(const TV& X) const override;
    T Integration_Step(const T phi) const override;
    T Minimum_Cell_Size() const override;
    int Active_Levelset(const TV& X) const;
    void Read(TYPED_ISTREAM input) override;
    void Write(TYPED_OSTREAM output) const override;
//#####################################################################
};
template<class TV,class ...Args>
inline IMPLICIT_OBJECT_UNION<TV>* Unite(IMPLICIT_OBJECT<TV>* o,Args&&... io_list)
{
    return new IMPLICIT_OBJECT_UNION<TV>(o,io_list...);
}
}
#endif
