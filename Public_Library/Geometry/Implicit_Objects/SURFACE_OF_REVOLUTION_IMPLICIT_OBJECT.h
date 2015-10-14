//#####################################################################
// Copyright 2007, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT__
#define __SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT__

#include <Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <Tools/Images/IMAGE.h>
#include <Geometry/Implicit_Objects/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class T_input>
class SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
    using IMPLICIT_OBJECT<TV>::box;
public:
    GRID<VECTOR<T,2> > slice_grid;
    ARRAY<T,VECTOR<int,2> > slice_phi;
    LEVELSET<VECTOR<T,2> > slice_levelset;
    T width,height,tolerance;

    SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT(std::string& filename,const T width_input,const T height_input,
        const T tolerance_input=1e-6);
    void Update_Box() override;
    T Integration_Step(const T phi) const override;
    T operator()(const TV& X) const override;
    TV Normal(const TV& X,const int aggregate=-1) const override;
    virtual void Read(TYPED_ISTREAM& input) override;
    virtual void Write(TYPED_OSTREAM& output) const override;
//#####################################################################
};
}
#endif
