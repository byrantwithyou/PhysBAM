//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SMOOTH_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __SMOOTH_LEVELSET_IMPLICIT_OBJECT__
#define __SMOOTH_LEVELSET_IMPLICIT_OBJECT__

#include <Tools/Grids_Uniform_Interpolation/CUBIC_SPLINE_INTERPOLATION_UNIFORM.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class SMOOTH_LEVELSET_IMPLICIT_OBJECT:public LEVELSET_IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<T,d-1> T_PRINCIPAL_CURVATURES;
public:
    typedef int HAS_TYPED_READ_WRITE;
    typedef LEVELSET_IMPLICIT_OBJECT<TV> BASE;
public:
    using BASE::levelset;using BASE::box;

    SMOOTH_LEVELSET_IMPLICIT_OBJECT(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input);
    virtual ~SMOOTH_LEVELSET_IMPLICIT_OBJECT();

//###########################################################################
    static SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>* Create();
    TV Normal(const TV& location,const int aggregate=-1) const override;
    TV Extended_Normal(const TV& location,const int aggregate=-1) const override;
    void Compute_Normals() override;
    SYMMETRIC_MATRIX<T,TV::m> Hessian(const TV& X) const override;
    VECTOR<T,d-1> Principal_Curvatures(const TV& X) const override;
    virtual std::string Name() const override {return Static_Name();}
    virtual std::string Extension() const override {return Static_Extension();}
    static std::string Static_Name();
    static std::string Static_Extension();
//###########################################################################
};
}
#endif
