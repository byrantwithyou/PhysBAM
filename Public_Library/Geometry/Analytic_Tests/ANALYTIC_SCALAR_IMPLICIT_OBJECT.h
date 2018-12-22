//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_SCALAR_IMPLICIT_OBJECT__
#define __ANALYTIC_SCALAR_IMPLICIT_OBJECT__

#include <Geometry/Analytic_Tests/ANALYTIC_SCALAR.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class TV>
struct ANALYTIC_SCALAR_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    ANALYTIC_SCALAR<TV> * as;
    T cell_size=0;
    T time=0;

    ANALYTIC_SCALAR_IMPLICIT_OBJECT(ANALYTIC_SCALAR<TV> * as): as(as) {}

    static ANALYTIC_SCALAR_IMPLICIT_OBJECT* Create()
    {return new ANALYTIC_SCALAR_IMPLICIT_OBJECT(0);}

    void Update_Box() override
    {}

    void Update_Minimum_Cell_Size(const int maximum_depth) override
    {}

    T operator()(const TV& location) const override
    {return as->f(location,time);}

    T Extended_Phi(const TV& location) const override
    {return as->f(location,time);}

    TV Normal(const TV& location,const int aggregate=-1) const override
    {return as->dX(location,time);}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const override
    {return as->dX(location,time);}

    void Compute_Normals() override
    {}

    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) override
    {}

    bool Lazy_Inside(const TV& location,const T contour_value=0) const override
    {return as->f(location,time)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override
    {phi_value=as->f(location,time);return phi_value<=contour_value;}

    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override
    {return Lazy_Inside(unclamped_X,contour_value);}

    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override
    {return Lazy_Inside_And_Value(unclamped_X,phi_value,contour_value);}

    bool Lazy_Outside(const TV& location,const T contour_value=0) const override
    {return !Lazy_Inside(location,contour_value);}

    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override
    {return !Lazy_Inside_Extended_Levelset(unclamped_X,contour_value);}

    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override
    {return !Lazy_Inside_Extended_Levelset_And_Value(unclamped_X,phi_value,contour_value);}

    virtual std::string Name() const override {return Static_Name();}
    static std::string Static_Name()
    {return "ANALYTIC_SCALAR_IMPLICIT_OBJECT<TV>";}

    std::string Extension() const override {return "";}

    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override
    {return cell_size;}

    T Minimum_Cell_Size() const override
    {return cell_size;}

    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override
    {return as->ddX(X,time);}

    void Read(TYPED_ISTREAM input) override
    {PHYSBAM_FATAL_ERROR();}

    void Write(TYPED_OSTREAM output) const override
    {PHYSBAM_FATAL_ERROR();}
};
template<class TV>
ANALYTIC_SCALAR_IMPLICIT_OBJECT<TV>* Make_IO(ANALYTIC_SCALAR<TV>* as)
{
    return new ANALYTIC_SCALAR_IMPLICIT_OBJECT<TV>{as};
}
}
#endif
