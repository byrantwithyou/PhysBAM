//#####################################################################
// Copyright 2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __ANALYTIC_IMPLICIT_OBJECT__
#define __ANALYTIC_IMPLICIT_OBJECT__

#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class T_ANALYTIC>
class ANALYTIC_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<typename T_ANALYTIC::VECTOR_T>
{
    typedef typename T_ANALYTIC::VECTOR_T TV;
    typedef typename TV::SCALAR T;
public:
    using IMPLICIT_OBJECT<TV>::box;
    T_ANALYTIC analytic;
    T cell_size;

    ANALYTIC_IMPLICIT_OBJECT(const T_ANALYTIC& analytic_input)
        :analytic(analytic_input),cell_size(0)
    {Update_Box();}

    static ANALYTIC_IMPLICIT_OBJECT* Create()
    {return new ANALYTIC_IMPLICIT_OBJECT(T_ANALYTIC());}

    void Update_Box() override
    {box=analytic.Bounding_Box();if(box!=RANGE<TV>::Empty_Box() && box!=RANGE<TV>::Full_Box()) box.Scale_About_Center((T)1.1);}

    void Update_Minimum_Cell_Size(const int maximum_depth) override
    {}

    T operator()(const TV& location) const override
    {return analytic.Signed_Distance(location);}

    T Extended_Phi(const TV& location) const override
    {return analytic.Signed_Distance(location);}

    TV Normal(const TV& location,const int aggregate=-1) const override
    {return analytic.Normal(location);}

    TV Extended_Normal(const TV& location,const int aggregate=-1) const override
    {return analytic.Normal(location);}

    void Compute_Normals() override
    {}

    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) override
    {}

    bool Lazy_Inside(const TV& location,const T contour_value=0) const override
    {return box.Inside(location,-contour_value)&&analytic.Signed_Distance(location)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override
    {if(box.Inside(location,-contour_value)){phi_value=analytic.Signed_Distance(location);if(phi_value<=contour_value) return true;}
    return false;}

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
    {return LOG::sprintf("ANALYTIC_IMPLICIT_OBJECT<%s>",T_ANALYTIC::Name().c_str());}

    std::string Extension() const override {return "";}

    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const override
    {return cell_size;}

    T Minimum_Cell_Size() const override
    {return cell_size;}

    VECTOR<T,TV::dimension-1> Principal_Curvatures(const TV& X) const override
    {return analytic.Principal_Curvatures(X);}

    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override
    {return analytic.Hessian(X);}

    void Read(TYPED_ISTREAM& input) override
    {Read_Binary(input,analytic);Update_Box();}

    void Write(TYPED_OSTREAM& output) const override
    {Write_Binary(output,analytic);}

//#####################################################################
};
}
#endif
