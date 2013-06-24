//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEPARABLE_UNIFORM
//#####################################################################
#ifndef __ADVECTION_SEPARABLE_UNIFORM__
#define __ADVECTION_SEPARABLE_UNIFORM__

#include <Tools/Advection/ADVECTION.h>
#include <Tools/Arrays/ARRAYS_FORWARD.h>
#include <Tools/Grids_Uniform_Advection/ADVECTION_UNIFORM_FORWARD.h>
#include <Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <Tools/Math_Tools/constants.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_AVERAGING> // T_AVERAGING=AVERAGING_UNIFORM<T_GRID>
class ADVECTION_SEPARABLE_UNIFORM:public ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;typedef typename T_AVERAGING::FACE_LOOKUP T_FACE_LOOKUP;
public:
    using ADVECTION<T_GRID,T2,typename T_AVERAGING::FACE_LOOKUP>::Update_Advection_Equation_Cell;

    ADVECTION_SEPARABLE_UNIFORM()
    {}

    virtual ~ADVECTION_SEPARABLE_UNIFORM()
    {}

    static T ENO(const T D1) // 1st order
    {return D1;}

    static T ENO(const T dx,const T D1,const T D2_left,const T D2_right) // 2rd order
    {return D1+minmag(D2_left,D2_right)*dx;}

    static T ENO(const T dx,const T D1,const T D2_left,const T D2_right,const T D3_left,const T D3_center,const T D3_right) // 3rd order
    {if(abs(D2_left) <= abs(D2_right)) return D1+dx*(D2_left+dx*2*minmag(D3_left,D3_center));else return D1+dx*(D2_right-dx*minmag(D3_right,D3_center));}

    static T WENO(const T v1,const T v2,const T v3,const T v4,const T v5,const T epsilon)
    {T s1=(T)thirteen_over_twelve*sqr(v1-2*v2+v3)+(T).25*sqr(v1-4*v2+3*v3),s2=(T)thirteen_over_twelve*sqr(v2-2*v3+v4)+(T).25*sqr(v2-v4),
        s3=(T)thirteen_over_twelve*sqr(v3-2*v4+v5)+(T).25*sqr(3*v3-4*v4+v5);
    T a1=(T).1/sqr(epsilon+s1),a2=(T).6/sqr(epsilon+s2),a3=(T).3/sqr(epsilon+s3),one_over_sum_a=1/(a1+a2+a3);
    T w1=a1*one_over_sum_a,w2=a2*one_over_sum_a,w3=1-w1-w2;
    return (T)one_sixth*(w1*(2*v1-7*v2+11*v3)+w2*(-v2+5*v3+2*v4)+w3*(2*v3+5*v4-v5));}

//#####################################################################
    virtual void Set_Order(const int order_input=3){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Set_Epsilon(const T epsilon_input=1e-6){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Set_Axis(const int axis){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Compute_Epsilon(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Update_Advection_Equation_Node(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const ARRAY<TV,TV_INT>& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& V,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0);
    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T>& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max);
    virtual void Advection_Solver(const int m,const T dx,const ARRAY<T2,VECTOR<int,1> >& Z,const ARRAY<T,VECTOR<int,1> >& u,ARRAY<T2,VECTOR<int,1> >& u_Zx){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
