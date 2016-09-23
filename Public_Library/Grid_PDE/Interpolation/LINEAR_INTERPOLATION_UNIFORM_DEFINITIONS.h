#include <Core/Matrices/DIAGONAL_MATRIX.h>
#include <Core/Matrices/MATRIX_1X1.h>
#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Function Clamped_To_Array_No_Extrema
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> T2 LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array_No_Extrema(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    TV X_clamped=clamp(X,grid.X(index),grid.X(index+1));
    return From_Base_Node(grid,u,X_clamped,index);
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> T2 LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    return From_Base_Node(grid,u,X,index);
}
//#####################################################################
// Function Clamped_To_Array_Weights
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<VECTOR<int,TV::m>,typename TV::SCALAR> > LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    return From_Base_Node_Weights(grid,u,X,index);
}
//#####################################################################
// Function Extrema_Clamped_To_Array
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> VECTOR<T2,2> LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Extrema_Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u_min,X);
    return Extrema_From_Base_Node(grid,u_min,u_max,X,index);
}
namespace{
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> T2
From_Base_Node_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    int i=index.x;
    return LINEAR_INTERPOLATION<T,T2>::Linear(u(index),u(i+1),(X-grid.X(index))*grid.one_over_dX);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> T2
From_Base_Node_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    int i=index.x,j=index.y;
    return LINEAR_INTERPOLATION<T,T2>::Bilinear(u(index),u(i+1,j),u(i,j+1),u(i+1,j+1),(X-grid.X(index))*grid.one_over_dX);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> T2
From_Base_Node_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    const T2 *base=&u(index),*basep1=base+u.stride.x;
    return LINEAR_INTERPOLATION<T,T2>::Trilinear(*base,*basep1,base[u.stride.y],basep1[u.stride.y],base[1],basep1[1],base[u.stride.y+1],basep1[u.stride.y+1],(X-grid.X(index))*grid.one_over_dX);
}
}
namespace{
//#####################################################################
// Function From_Base_Node_Weights_Helper
//#####################################################################
template<class T,class T2> ARRAY<PAIR<VECTOR<int,1>,T> >
From_Base_Node_Weights_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    typedef VECTOR<int,1> TV_INT;
    int i=index.x;
    T x=((X-grid.X(index))*grid.one_over_dX).x;
    ARRAY<PAIR<TV_INT,T> > weights;
    weights.Append(PAIR<TV_INT,T>(index,1-x));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1),x));
    return weights;
}
//#####################################################################
// Function From_Base_Node_Weights_Helper
//#####################################################################
template<class T,class T2> ARRAY<PAIR<VECTOR<int,2>,T> >
From_Base_Node_Weights_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,2> TV;
    int i=index.x;
    int j=index.y;
    TV w=(X-grid.X(index))*grid.one_over_dX;
    ARRAY<PAIR<TV_INT,T> > weights;
    weights.Append(PAIR<TV_INT,T>(index,(1-w.x)*(1-w.y)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j),w.x*(1-w.y)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j+1),(1-w.x)*w.y));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j+1),w.x*w.y));
    return weights;
}
//#####################################################################
// Function From_Base_Node_Weights_Helper
//#####################################################################
template<class T,class T2> ARRAY<PAIR<VECTOR<int,3>,T> >
From_Base_Node_Weights_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    typedef VECTOR<int,3> TV_INT;typedef VECTOR<T,3> TV;
    int i=index.x;
    int j=index.y;
    int k=index.z;
    TV w=(X-grid.X(index))*grid.one_over_dX;
    ARRAY<PAIR<TV_INT,T> > weights;
    weights.Append(PAIR<TV_INT,T>(index,(1-w.x)*(1-w.y)*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j,k),w.x*(1-w.y)*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j+1,k),(1-w.x)*w.y*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j,k+1),(1-w.x)*(1-w.y)*w.z));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j+1,k),w.x*w.y*(1-w.z)));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j,k+1),w.x*(1-w.y)*w.z));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i,j+1,k+1),(1-w.x)*w.y*w.z));
    weights.Append(PAIR<TV_INT,T>(TV_INT(i+1,j+1,k+1),w.x*w.y*w.z));
    return weights;
}
//#####################################################################
// Function Extrema_From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,2>
Extrema_From_Base_Node_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u_min,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u_max,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    int i=index.x;
    return VECTOR<T2,2>(Componentwise_Min(u_min(i),u_min(i+1)),Componentwise_Max(u_max(i),u_max(i+1)));
}
//#####################################################################
// Function Extrema_From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,2>
Extrema_From_Base_Node_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u_min,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u_max,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    int i=index.x,j=index.y;
    return VECTOR<T2,2>(Componentwise_Min(u_min(i,j),u_min(i+1,j),u_min(i,j+1),u_min(i+1,j+1)),Componentwise_Max(u_max(i,j),u_max(i+1,j),u_max(i,j+1),u_max(i+1,j+1)));
}
//#####################################################################
// Function Extrema_From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,2>
Extrema_From_Base_Node_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u_min,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u_max,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    const T2 *min_base=&u_min(index),*max_base=&u_max(index),*min_basep1=min_base+u_min.stride.x,*max_basep1=max_base+u_max.stride.x;
    return VECTOR<T2,2>(
        Componentwise_Min(*min_base,*min_basep1,min_base[u_min.stride.y],min_basep1[u_min.stride.y],min_base[1],min_basep1[1],min_base[u_min.stride.y+1],min_basep1[u_min.stride.y+1]),
        Componentwise_Max(*max_base,*max_basep1,max_base[u_max.stride.y],max_basep1[u_max.stride.y],max_base[1],max_basep1[1],max_base[u_max.stride.y+1],max_basep1[u_max.stride.y+1]));
}
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> T2 LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const
{
    return From_Base_Node_Helper(grid,u,X,index);
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<VECTOR<int,TV::m>,typename TV::SCALAR> > LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X,const TV_INT& index) const
{
    return From_Base_Node_Weights_Helper(grid,u,X,index);
}
//#####################################################################
// Function Extrema_From_Base_Node
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> VECTOR<T2,2> LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Extrema_From_Base_Node(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u_min,const ARRAYS_ND_BASE<T2,TV_INT>& u_max,const TV& X,const TV_INT& index) const
{
    return Extrema_From_Base_Node_Helper(grid,u_min,u_max,X,index);
}
//#####################################################################
// Function From_Block_Face_Component
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> typename TV::SCALAR LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
{
    return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Interpolate_Face_Component(axis,block,u,X);
}
//#####################################################################
// Function From_Block_Face_Component_Weights
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<FACE_INDEX<TV::m>,typename TV::SCALAR> > LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
From_Block_Face_Component_Weights(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
{
    return LINEAR_INTERPOLATION_MAC_HELPER<TV>::Interpolate_Face_Component_Weights(axis,block,u,X);
}
namespace{
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,1>
From_Base_Node_Gradient_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& w,const VECTOR<int,1>& index)
{
    return LINEAR_INTERPOLATION<T,T2>::Linear_Gradient(u(index),u(index.x+1),w);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,2>
From_Base_Node_Gradient_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& w,const VECTOR<int,2>& index)
{
    int i=index.x,j=index.y;
    return LINEAR_INTERPOLATION<T,T2>::Bilinear_Gradient(u(index),u(i+1,j),u(i,j+1),u(i+1,j+1),w);
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> VECTOR<T2,3>
From_Base_Node_Gradient_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& w,const VECTOR<int,3>& index)
{
    const T2 *base=&u(index),*basep1=base+u.stride.x;
    return LINEAR_INTERPOLATION<T,T2>::Trilinear_Gradient(*base,*basep1,base[u.stride.y],basep1[u.stride.y],base[1],basep1[1],base[u.stride.y+1],basep1[u.stride.y+1],w);
}
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> VECTOR<T2,TV::m> LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Gradient(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    VECTOR<T2,TV::m> V=From_Base_Node_Gradient_Helper(grid,u,(X-grid.X(index))*grid.one_over_dX,index);
    for(int i=0;i<TV::m;i++) V(i)*=grid.one_over_dX(i);
    return V;
}
namespace{
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> SYMMETRIC_MATRIX<T2,1>
From_Base_Node_Hessian_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& w,const VECTOR<int,1>& index)
{
    int i=index.x;
    SYMMETRIC_MATRIX<T2,1> S=LINEAR_INTERPOLATION<T,T2>::Linear_Hessian(u(index),u(i+1),w);
    S.x00*=sqr(grid.one_over_dX.x);
    return S;
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> SYMMETRIC_MATRIX<T2,2>
From_Base_Node_Hessian_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& w,const VECTOR<int,2>& index)
{
    int i=index.x,j=index.y;
    SYMMETRIC_MATRIX<T2,2> S=LINEAR_INTERPOLATION<T,T2>::Bilinear_Hessian(u(index),u(i+1,j),u(i,j+1),u(i+1,j+1),w);
    S.x00*=sqr(grid.one_over_dX.x);
    S.x10*=grid.one_over_dX.x*grid.one_over_dX.y;
    S.x11*=sqr(grid.one_over_dX.y);
    return S;
}
//#####################################################################
// Function From_Base_Node_Helper
//#####################################################################
template<class T,class T2> SYMMETRIC_MATRIX<T2,3>
From_Base_Node_Hessian_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& w,const VECTOR<int,3>& index)
{
    const T2 *base=&u(index),*basep1=base+u.stride.x;
    SYMMETRIC_MATRIX<T2,3> S=LINEAR_INTERPOLATION<T,T2>::Trilinear_Hessian(*base,*basep1,base[u.stride.y],basep1[u.stride.y],base[1],basep1[1],base[u.stride.y+1],basep1[u.stride.y+1],w);
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
// Function Clamped_To_Array
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> SYMMETRIC_MATRIX<T2,TV::m> LINEAR_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Hessian(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    TV_INT index=INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_End_Minus_One(grid,u,X);
    return From_Base_Node_Hessian_Helper(grid,u,(X-grid.X(index))*grid.one_over_dX,index);
}

