//#####################################################################
// Copyright 2002-2009, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIBODY_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __MULTIBODY_LEVELSET_IMPLICIT_OBJECT__
#define __MULTIBODY_LEVELSET_IMPLICIT_OBJECT__

#include <Tools/Interpolation/INTERPOLATION_FORWARD.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/MATRIX_POLICY.h>
#include <Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class TV>
class MULTIBODY_LEVELSET_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
public:
    typedef IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;

    ARRAY<IMPLICIT_OBJECT<TV>*>* levelsets;
protected:
    T minimum_cell_size;
public:
    bool need_destroy_data;

    MULTIBODY_LEVELSET_IMPLICIT_OBJECT(ARRAY<IMPLICIT_OBJECT<TV>*>* levelsets_input);
    virtual ~MULTIBODY_LEVELSET_IMPLICIT_OBJECT();

    static MULTIBODY_LEVELSET_IMPLICIT_OBJECT* Create();
    void Delete_Pointers_And_Clean_Memory();
    void Update_Box() PHYSBAM_OVERRIDE;
    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE;
    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box) const  PHYSBAM_OVERRIDE; // TODO: make this check overlap with the grids.
    T operator()(const TV& location) const PHYSBAM_OVERRIDE;
    T Phi_With_Index(const TV& location,int& levelset_index) const;
    T Extended_Phi(const TV& location) const PHYSBAM_OVERRIDE;
    T Phi_Secondary(const TV& location) const PHYSBAM_OVERRIDE; // TODO: implement this, but need Phi_Secondary_Extended()!
    TV Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE;
    TV Extended_Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE;
    void Compute_Normals()  PHYSBAM_OVERRIDE;
    void Compute_Cell_Minimum_And_Maximum(const bool recompute_if_exists=true) PHYSBAM_OVERRIDE;
    void Inflate(const T inflation_distance) PHYSBAM_OVERRIDE;
    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE;
    T Minimum_Cell_Size() const  PHYSBAM_OVERRIDE;
    bool Lazy_Inside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Inside_And_Value(const TV& location,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Outside(const TV& location,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const PHYSBAM_OVERRIDE;
    T Min_Phi() const PHYSBAM_OVERRIDE;
    void Rescale(const T scaling_factor) PHYSBAM_OVERRIDE;
    void Translate(const TV& translation) PHYSBAM_OVERRIDE;
    VECTOR<T,TV::dimension-1> Principal_Curvatures(const TV& X) const PHYSBAM_OVERRIDE;
    virtual std::string Name() const PHYSBAM_OVERRIDE;
    static std::string Static_Name();
    virtual std::string Extension() const PHYSBAM_OVERRIDE;
    static std::string Static_Extension();
    void Read(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE; // TODO -- fix to read/write levelsets
    void Write(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE;
//###########################################################################
};
}
#endif
