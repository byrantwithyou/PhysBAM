//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_LAGRANGE_2D
//#####################################################################
#ifndef __GRID_LAGRANGE_2D__
#define __GRID_LAGRANGE_2D__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Compressible/Lagrange_Equations/GRID_LAGRANGE.h>
namespace PhysBAM{

template<class T>
class GRID_LAGRANGE_2D:public GRID_LAGRANGE<T>
{
    typedef VECTOR<int,2> TV_INT;
    typedef VECTOR<T,2> TV;
public:
    int m;       // # of points: dimension 1
    int n;       // # of points: dimension 2
    ARRAY<T,TV_INT>& x; // x spatial location
    ARRAY<T,TV_INT>& y; // y spatial location

    GRID_LAGRANGE_2D(const int m_input,const int n_input,ARRAY<T,TV_INT>& x_input,ARRAY<T,TV_INT>& y_input)
        :m(m_input),n(n_input),x(x_input),y(y_input)
    {}

//#####################################################################
    void Euler_Step(const ARRAY<T,TV_INT>& u,const ARRAY<T,TV_INT>& v,const T dt);
    void Get_Lengths(ARRAY<T,TV_INT>& L0,ARRAY<T,TV_INT>& L1);
    void Get_Areas(ARRAY<T,TV_INT>& A);
    void Get_Normals(ARRAY<T,TV_INT>& N1_x,ARRAY<T,TV_INT>& N1_y,ARRAY<T,TV_INT>& N2_x,ARRAY<T,TV_INT>& N2_y);
    void Get_Centers(ARRAY<T,TV_INT>& C_x,ARRAY<T,TV_INT>& C_y);
    void Get_Midpoints(ARRAY<T,TV_INT>& M1_x,ARRAY<T,TV_INT>& M1_y,ARRAY<T,TV_INT>& M2_x,ARRAY<T,TV_INT>& M2_y);
    void Get_Sub_Zone_Lengths(ARRAY<T,TV_INT>& LL0,ARRAY<T,TV_INT>& LL1,ARRAY<T,TV_INT>& LL2,ARRAY<T,TV_INT>& LL3);
    void Get_Sub_Zone_Areas(ARRAY<T,TV_INT>& AA0,ARRAY<T,TV_INT>& AA1,ARRAY<T,TV_INT>& AA2,ARRAY<T,TV_INT>& AA3);
    void Get_Sub_Zone_Normals(ARRAY<T,TV_INT>& NN1_x,ARRAY<T,TV_INT>& NN1_y,ARRAY<T,TV_INT>& NN2_x,ARRAY<T,TV_INT>& NN2_y,ARRAY<T,TV_INT>& NN3_x,ARRAY<T,TV_INT>& NN3_y,ARRAY<T,TV_INT>& NN4_x,ARRAY<T,TV_INT>& NN4_y);
//#####################################################################
};
}
#endif
