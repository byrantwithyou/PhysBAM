//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BSPLINE
//#####################################################################
#ifndef BSPLINE__
#define BSPLINE__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/clamp.h>
#include <Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T,class T2>
class BSPLINE
{
public:
    ARRAY<T> control_points_times;
    ARRAY<T2> control_points;
    int k;
    bool closed;

    BSPLINE(const ARRAY<T>& control_points_times,const ARRAY<T2>& control_points,const int order=1);
    virtual ~BSPLINE();

    T Range()
    {return End_Time()-Start_Time();}

    T Start_Time()
    {return control_points_times(k);}

    T End_Time()
    {return control_points_times(control_points_times.m-(k-1));}

//#####################################################################
    virtual T2 Evaluate(const T t);
    virtual T2 Clamped_Evaluate(const T t);
    T Basis_Function(const int i,const int k,const T t);
    void Clamp_End_Points();
    virtual void Create_Closed_Points();
    void Normalize_Control_Points();
    void Print_Control_Points_And_Times();
//#####################################################################
};
}
#endif
