//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Data_Structures/PAIR.h>
#include <Tools/Grids_Uniform_Interpolation/CUBIC_SPLINE_INTERPOLATION_UNIFORM.h>
#include <Tools/Matrices/DIAGONAL_MATRIX.h>
using namespace PhysBAM;
//#####################################################################
// Function Interpolate
//#####################################################################
template<class T,class T2> static T2
Cubic_Interpolation(const T2 x[4],T a)
{
    T2 d20=x[2]-x[0],d30=x[3]-x[0],d21=x[2]-x[1];
    return (T2).5*(((d30-3*d21)*a-(d30-5*d21+d20))*a+d20)*a+x[1];
}
//#####################################################################
// Function Cubic_Interpolation_Diff
//#####################################################################
template<class T,class T2> static T2
Cubic_Interpolation_Diff(const T2 x[4],T a)
{
    T2 d20=x[2]-x[0],d30=x[3]-x[0],d21=x[2]-x[1];
    return (T2).5*((3*(d30-3*d21)*a-2*(d30-5*d21+d20))*a+d20);
}
//#####################################################################
// Function Cubic_Interpolation_Taylor
//#####################################################################
template<class T,class T2> static void
Cubic_Interpolation_Taylor(const T2 x[4],T a,T2& fa,T2& dfa)
{
    T2 d20=x[2]-x[0],d30=x[3]-x[0],d21=x[2]-x[1];
    T2 r=d30-3*d21;
    T2 s=r-2*d21+d20;
    fa=(T2).5*((r*a-s)*a+d20)*a+x[1];
    dfa=(T2).5*((3*r*a-2*s)*a+d20);
}
//#####################################################################
// Function Cubic_Interpolation_Taylor
//#####################################################################
template<class T,class T2> static void
Cubic_Interpolation_Taylor(const T2 x[4],T a,T2& fa,T2& dfa,T2& ddfa)
{
    T2 d20=x[2]-x[0],d30=x[3]-x[0],d21=x[2]-x[1];
    T2 r=d30-3*d21;
    T2 s=r-2*d21+d20;
    fa=(T2).5*((r*a-s)*a+d20)*a+x[1];
    ddfa=3*r*a-s;
    dfa=(T2).5*((ddfa-s)*a+d20);
}
//#####################################################################
// Function Cubic_Interpolation_Diff2
//#####################################################################
template<class T,class T2> static T2
Cubic_Interpolation_Diff2(const T2 x[4],T a)
{
    T2 d20=x[2]-x[0],d30=x[3]-x[0],d21=x[2]-x[1];
    return 3*(d30-3*d21)*a-(d30-5*d21+d20);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
CUBIC_SPLINE_INTERPOLATION_UNIFORM()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
~CUBIC_SPLINE_INTERPOLATION_UNIFORM()
{
}
namespace{
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class TV,class T2> static T2
From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const TV& X,const VECTOR<int,1>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    return Cubic_Interpolation(&u(index-1),w.x);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class TV,class T2> static T2
From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const TV& X,const VECTOR<int,2>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    const T2* b=&u(index-1);
    T2 x[4];
    for(int i=0;i<4;i++,b+=u.stride.x)
        x[i]=Cubic_Interpolation(b,w.y);
    return Cubic_Interpolation(x,w.x);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class TV,class T2> static T2
From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const TV& X,const VECTOR<int,3>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    int dx=u.stride.x;
    const T2* b=&u(index-1),*c=b;
    T2 x[4],y[4];
    for(int i=0;i<4;i++,b+=dx,c=b){
        for(int j=0;j<4;j++,c+=u.stride.y)
            y[j]=Cubic_Interpolation(c,w.z);
        x[i]=Cubic_Interpolation(y,w.y);}
    return Cubic_Interpolation(x,w.x);
}
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> T2 CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const
{
    return From_Base_Node_Helper(grid,u,X,index);
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> T2 CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    return From_Base_Node_Helper(grid,u,X,INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X));
}
namespace{
//#####################################################################
// Function From_Base_Node_Gradient_Helper
//#####################################################################
template<class TV,class T2> static VECTOR<T2,1>
From_Base_Node_Gradient_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const TV& X,const VECTOR<int,1>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    return VECTOR<T2,1>(Cubic_Interpolation_Diff(&u(index-1),w.x))*grid.one_over_dX;
}
//#####################################################################
// Function From_Base_Node_Gradient
//#####################################################################
template<class TV,class T2> static VECTOR<T2,2>
From_Base_Node_Gradient_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const TV& X,const VECTOR<int,2>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    const T2* b=&u(index-1);
    T2 x[4],x_y[4];
    for(int i=0;i<4;i++,b+=u.stride.x)
        Cubic_Interpolation_Taylor(b,w.y,x[i],x_y[i]);
    return VECTOR<T2,2>(Cubic_Interpolation_Diff(x,w.x),Cubic_Interpolation(x_y,w.x))*grid.one_over_dX;
}
//#####################################################################
// Function From_Base_Node_Gradient
//#####################################################################
template<class TV,class T2> static VECTOR<T2,3>
From_Base_Node_Gradient_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const TV& X,const VECTOR<int,3>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    const T2* b=&u(index-1),*c=b;
    T2 x[4],y[4],x_z[4],y_z[4],x_y[4];
    for(int i=0;i<4;i++,b+=u.stride.x,c=b){
        for(int j=0;j<4;j++,c+=u.stride.y)
            Cubic_Interpolation_Taylor(c,w.z,y[j],y_z[j]);
        Cubic_Interpolation_Taylor(y,w.y,x[i],x_y[i]);
        x_z[i]=Cubic_Interpolation(y_z,w.y);}
    return VECTOR<T2,3>(Cubic_Interpolation_Diff(x,w.x),Cubic_Interpolation(x_y,w.x),Cubic_Interpolation(x_z,w.x))*grid.one_over_dX;
}
}
//#####################################################################
// Function Clamped_To_Array_Gradient
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> VECTOR<T2,TV::m> CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Gradient(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    return From_Base_Node_Gradient_Helper(grid,u,X,INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X));
}
namespace{
//#####################################################################
// Function From_Base_Node_Hessian_Helper
//#####################################################################
template<class TV,class T2> static SYMMETRIC_MATRIX<T2,1>
From_Base_Node_Hessian_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const TV& X,const VECTOR<int,1>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    SYMMETRIC_MATRIX<T2,1> S(Cubic_Interpolation_Diff2(&u(index-1),w.x));
    S.x00*=sqr(grid.one_over_dX.x);
    return S;
}
//#####################################################################
// Function From_Base_Node_Hessian
//#####################################################################
template<class TV,class T2> static SYMMETRIC_MATRIX<T2,2>
From_Base_Node_Hessian_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const TV& X,const VECTOR<int,2>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    const T2* b=&u(index-1);
    T2 x[4],x_y[4],x_yy[4];
    for(int i=0;i<4;i++,b+=u.stride.x)
        Cubic_Interpolation_Taylor(b,w.y,x[i],x_y[i],x_yy[i]);
    SYMMETRIC_MATRIX<T2,2> S(Cubic_Interpolation_Diff2(x,w.x),Cubic_Interpolation_Diff(x_y,w.x),Cubic_Interpolation(x_yy,w.x));
    S.x00*=sqr(grid.one_over_dX.x);
    S.x10*=grid.one_over_dX.x*grid.one_over_dX.y;
    S.x11*=sqr(grid.one_over_dX.y);
    return S;
}
//#####################################################################
// Function From_Base_Node_Hessian
//#####################################################################
template<class TV,class T2> static SYMMETRIC_MATRIX<T2,3>
From_Base_Node_Hessian_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const TV& X,const VECTOR<int,3>& index)
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    const T2* b=&u(index-1),*c=b;
    T2 x[4],y[4],x_z[4],y_z[4],x_y[4],x_zz[4],y_zz[4],x_yy[4],x_yz[4];
    for(int i=0;i<4;i++,b+=u.stride.x,c=b){
        for(int j=0;j<4;j++,c+=u.stride.y)
            Cubic_Interpolation_Taylor(c,w.z,y[j],y_z[j],y_zz[j]);
        Cubic_Interpolation_Taylor(y,w.y,x[i],x_y[i],x_yy[i]);
        Cubic_Interpolation_Taylor(y_z,w.y,x_z[i],x_yz[i]);
        x_zz[i]=Cubic_Interpolation(y_zz,w.y);}
    SYMMETRIC_MATRIX<T2,3> S(
        Cubic_Interpolation_Diff2(x,w.x),
        Cubic_Interpolation_Diff(x_y,w.x),
        Cubic_Interpolation_Diff(x_z,w.x),
        Cubic_Interpolation(x_yy,w.x),
        Cubic_Interpolation(x_yz,w.x),
        Cubic_Interpolation(x_zz,w.x));
    S.x00*=sqr(grid.one_over_dX.x);
    S.x11*=sqr(grid.one_over_dX.y);
    S.x22*=sqr(grid.one_over_dX.z);
    S.x10*=grid.one_over_dX.x*grid.one_over_dX.y;
    S.x20*=grid.one_over_dX.x*grid.one_over_dX.z;
    S.x21*=grid.one_over_dX.y*grid.one_over_dX.z;
    return S;
}
}
//#####################################################################
// Function Clamped_To_Array_Hessian
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> SYMMETRIC_MATRIX<T2,TV::m> CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Hessian(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    return From_Base_Node_Hessian_Helper(grid,u,X,INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X));
}
namespace PhysBAM{
template class CUBIC_SPLINE_INTERPOLATION_UNIFORM<VECTOR<float,1>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > >;
template class CUBIC_SPLINE_INTERPOLATION_UNIFORM<VECTOR<float,2>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class CUBIC_SPLINE_INTERPOLATION_UNIFORM<VECTOR<float,3>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
template class CUBIC_SPLINE_INTERPOLATION_UNIFORM<VECTOR<double,1>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > >;
template class CUBIC_SPLINE_INTERPOLATION_UNIFORM<VECTOR<double,2>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class CUBIC_SPLINE_INTERPOLATION_UNIFORM<VECTOR<double,3>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
}
