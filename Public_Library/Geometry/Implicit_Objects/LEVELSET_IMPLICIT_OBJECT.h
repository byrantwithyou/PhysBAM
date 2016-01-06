//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __LEVELSET_IMPLICIT_OBJECT__
#define __LEVELSET_IMPLICIT_OBJECT__

#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class LEVELSET_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<T,d-1> T_PRINCIPAL_CURVATURES;
public:
    typedef int HAS_TYPED_READ_WRITE;
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    LEVELSET<TV> levelset;
    ARRAY<TV,TV_INT>* V;
    INTERPOLATION_UNIFORM<TV,TV>* velocity_interpolation;
private:
    static LINEAR_INTERPOLATION_UNIFORM<TV,TV> default_velocity_interpolation;
protected:
    T minimum_cell_size;
    bool need_destroy_data;
public:

    LEVELSET_IMPLICIT_OBJECT(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input);
    virtual ~LEVELSET_IMPLICIT_OBJECT();

    void Read(TYPED_ISTREAM& input) override
    {Read_Binary(input,levelset);Update_Box();Update_Minimum_Cell_Size();}

    void Write(TYPED_OSTREAM& output) const override
    {Write_Binary(output,levelset);}

//###########################################################################
    static LEVELSET_IMPLICIT_OBJECT<TV>* Create();
    void Set_Custom_Secondary_Interpolation(INTERPOLATION_UNIFORM<TV,T>& interpolation);
    void Set_Custom_Normal_Interpolation(INTERPOLATION_UNIFORM<TV,TV>& interpolation);
    void Set_Custom_Velocity_Interpolation(INTERPOLATION_UNIFORM<TV,TV>& interpolation);
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
    void Inflate(const T inflation_distance) override;
    bool Lazy_Inside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override;
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override;
    bool Lazy_Outside(const TV& location,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const override;
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const override;
    TV Velocity(const TV& location) const override;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override;
    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const override;
    void Rescale(const T scaling_factor) override;
    void Translate(const TV& translation) override;
    virtual std::string Name() const override {return Static_Name();}
    virtual std::string Extension() const override {return Static_Extension();}
    static std::string Static_Name();
    static std::string Static_Extension();
    T Integration_Step(const T phi) const override;
    T Minimum_Cell_Size() const override;
//###########################################################################
};
}
#endif
