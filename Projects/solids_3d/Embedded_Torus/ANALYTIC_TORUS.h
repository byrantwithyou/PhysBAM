//#####################################################################
// Copyright 2003-2004, Geoffrey Irving, Neil "Rock and Roll" Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_TORUS
//#####################################################################
#ifndef __ANALYTIC_TORUS__
#define __ANALYTIC_TORUS__

#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template <class T>
class ANALYTIC_TORUS
{
public:
    T Radius;
    T radius;
    VECTOR_3D<T> center;
    BOX_3D<T> bounding_box;

    ANALYTIC_TORUS(){Radius=(T).6;radius=(T).3;Update_Bounding_Box();}

    ANALYTIC_TORUS(VECTOR_3D<T> &center_input,T Radius_input,T radius_input):center(center_input),Radius(Radius_input),radius(radius_input)
    {}

    ~ANALYTIC_TORUS(){}

    void Update_Bounding_Box()
    {T expansion_ratio=(T)1.1;bounding_box.xmin=center.x-(radius+Radius);bounding_box.xmax=center.x+(radius+Radius);   
    bounding_box.ymin=center.y-(radius+Radius);bounding_box.ymax=center.y+(radius+Radius);bounding_box.zmin=center.z-radius;bounding_box.zmax=center.z+radius;
    bounding_box.Scale_About_Center(expansion_ratio);}

    //d = r - sqrt(x*x + y*y + z*z + R*R - 2*R*sqrt(x*x + y*y))
    T Phi(const VECTOR_3D<T> &x) const
    {VECTOR_3D<T> y=x-center;return (sqrt(y.Magnitude_Squared() + Radius*Radius - (T)2*Radius*sqrt(sqr(y.x)+sqr(y.y))) - radius);}

//#####################################################################
};
}
#endif
