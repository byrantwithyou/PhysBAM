#include <Core/Data_Structures/PAIR.h>
#include <Grid_PDE/Interpolation/CUBIC_MONOTONIC_INTERPOLATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
CUBIC_MONOTONIC_INTERPOLATION_UNIFORM()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
~CUBIC_MONOTONIC_INTERPOLATION_UNIFORM()
{
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP,class T> T2
From_Base_Node(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(u(index.x),u(index.x+1),u(index.x+2),u(index.x+3),(X.x-grid.X(index+1).x)*grid.one_over_dX.x);
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP,class T> ARRAY<PAIR<VECTOR<int,TV::m>,typename TV::SCALAR> >
From_Base_Node_Weights(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index)
{
    typedef VECTOR<int,1> TV_INT;
    ARRAY<PAIR<TV_INT,T> > weights;
    PHYSBAM_NOT_IMPLEMENTED();
    return weights;
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP,class T> T2
From_Base_Node(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    VECTOR<T,2> alpha=X-grid.X(index+VECTOR<int,2>::All_Ones_Vector());alpha*=grid.one_over_dX;
    T2 interpolated_in_x[4],value[4];
    for(int b=0;b<4;b++){
        for(int a=0;a<4;a++) value[a]=u(index.x+a,index.y+b);
        interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(value[0],value[1],value[2],value[3],alpha.x);}
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP,class T> ARRAY<PAIR<VECTOR<int,TV::m>,typename TV::SCALAR> >
From_Base_Node_Weights(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index)
{
    typedef VECTOR<int,2> TV_INT;
    ARRAY<PAIR<TV_INT,T> > weights;
    PHYSBAM_NOT_IMPLEMENTED();
    return weights;
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP,class T> ARRAY<PAIR<VECTOR<int,TV::m>,typename TV::SCALAR> >
From_Base_Node_Weights(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index)
{
    typedef VECTOR<int,3> TV_INT;
    ARRAY<PAIR<TV_INT,T> > weights;
    PHYSBAM_NOT_IMPLEMENTED();
    return weights;
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP,class T> T2
From_Base_Node(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,
    const VECTOR<int,3>& index)
{
    VECTOR<T,3> alpha=X-grid.X(index+VECTOR<int,3>::All_Ones_Vector());alpha*=grid.one_over_dX;
    T2 interpolated_in_x[4],interpolated_in_y[4],value[4];
    for(int c=0;c<4;c++){
        for(int b=0;b<4;b++){
            for(int a=0;a<4;a++) value[a]=u(index.x+a,index.y+b,index.z+c);
            interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(value[0],value[1],value[2],value[3],alpha.x);}
        interpolated_in_y[c]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);}
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_y[0],interpolated_in_y[1],interpolated_in_y[2],interpolated_in_y[3],alpha.z);
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> T2 CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    return ::From_Base_Node(*this,grid,u,X,INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X)-TV_INT::All_Ones_Vector());
}
//#####################################################################
// Function Clamped_To_Array_Weights
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<VECTOR<int,TV::m>,typename TV::SCALAR> > CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Weights(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    return ::From_Base_Node_Weights(*this,grid,u,X,INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X)-TV_INT::All_Ones_Vector());
}
//#####################################################################
// Function From_Base_Node_Periodic
//#####################################################################
template<class T,class TV,class T2,class T_FACE_LOOKUP> T2
    From_Base_Node_Periodic(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,1> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,
        const VECTOR<int,1>& index)
{
    VECTOR<T,1> alpha=(X-grid.X(index+1))*grid.one_over_dX;
    T2 value[4];
    VECTOR<int,1> I;
    for(int a=0;a<4;a++){
        I.x=index.x+a;
        while(I.x<0) I.x+=grid.counts.x;
        while(I.x>=grid.counts.x) I.x-=grid.counts.x;
        value[a]=u(I);}
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(value[0],value[1],value[2],value[3],alpha.x);
}
//#####################################################################
// Function From_Base_Node_Periodic
//#####################################################################
template<class T,class TV,class T2,class T_FACE_LOOKUP> T2
    From_Base_Node_Periodic(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,2> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,
        const VECTOR<int,2>& index)
{
    VECTOR<T,2> alpha=(X-grid.X(index+1))*grid.one_over_dX;
    T2 interpolated_in_x[4],value[4];
    VECTOR<int,2> I;
    for(int b=0;b<4;b++){
        I.y=index.y+b;
        while(I.y<0) I.y+=grid.counts.y;
        while(I.y>=grid.counts.y) I.y-=grid.counts.y;
        for(int a=0;a<4;a++){
            I.x=index.x+a;
            while(I.x<0) I.x+=grid.counts.x;
            while(I.x>=grid.counts.x) I.x-=grid.counts.x;
            value[a]=u(I);}
        interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(value[0],value[1],value[2],value[3],alpha.x);}
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);
}
//#####################################################################
// Function From_Base_Node_Periodic
//#####################################################################
template<class T,class TV,class T2,class T_FACE_LOOKUP> T2
    From_Base_Node_Periodic(const CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>& cub,const GRID<VECTOR<T,3> >& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,
        const VECTOR<int,3>& index)
{
    VECTOR<T,3> alpha=(X-grid.X(index+1))*grid.one_over_dX;
    T2 interpolated_in_x[4],interpolated_in_y[4],value[4];
    VECTOR<int,3> I;
    for(int c=0;c<4;c++){
        I.z=index.z+c;
        while(I.z<0) I.z+=grid.counts.z;
        while(I.z>=grid.counts.z) I.z-=grid.counts.z;
        for(int b=0;b<4;b++){
            I.y=index.y+b;
            while(I.y<0) I.y+=grid.counts.y;
            while(I.y>=grid.counts.y) I.y-=grid.counts.y;
            for(int a=0;a<4;a++){
                I.x=index.x+a;
                while(I.x<0) I.x+=grid.counts.x;
                while(I.x>=grid.counts.x) I.x-=grid.counts.x;
                value[a]=u(I);}
            interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(value[0],value[1],value[2],value[3],alpha.x);}
        interpolated_in_y[c]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);}
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_y[0],interpolated_in_y[1],interpolated_in_y[2],interpolated_in_y[3],alpha.z);
}
//#####################################################################
// Function Periodic
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> T2 CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Periodic(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    return ::From_Base_Node_Periodic(*this,grid,u,X,grid.Index(X)-1);
}
//#####################################################################
// Function From_Block_Face_Component_Helper
//#####################################################################
template<class TV,class TV_INT,class T,class LOOKUP> static T
From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const LOOKUP& u,const VECTOR<T,1>& X,const TV_INT& index)
{
    CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T,LOOKUP> cub;
    TV_INT indicies[4];
    for(int i=0;i<4;i++) indicies[i]=index+i;
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(u(axis,indicies[0]),u(axis,indicies[1]),u(axis,indicies[2]),u(axis,indicies[3]),(X.x-grid.Face(FACE_INDEX<TV::m>(axis,index+1)).x)*grid.one_over_dX.x);
}
//#####################################################################
// Function From_Block_Face_Component_Helper
//#####################################################################
template<class TV,class TV_INT,class T,class LOOKUP> static T
From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const LOOKUP& u,const VECTOR<T,2>& X,const TV_INT& index)
{
    CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T,LOOKUP> cub;
    VECTOR<T,2> alpha=X-grid.Face(FACE_INDEX<TV::m>(axis,index+VECTOR<int,2>::All_Ones_Vector()));
    alpha*=grid.one_over_dX;
    T interpolated_in_x[4];
    T value[4];
    for(int b=0;b<4;b++){
        for(int a=0;a<4;a++){
            VECTOR<int,2> ind(index.x+a,index.y+b);
            value[a]=u(axis,ind);}
        interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(value[0],value[1],value[2],value[3],alpha.x);}
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);
}
//#####################################################################
// Function From_Block_Face_Component_Helper
//#####################################################################
template<class TV,class TV_INT,class T,class LOOKUP> static T
From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const LOOKUP& u,const VECTOR<T,3>& X,const TV_INT& index)
{
    CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T,LOOKUP> cub;
    VECTOR<T,3> alpha=X-grid.Face(FACE_INDEX<TV::m>(axis,index+VECTOR<int,3>::All_Ones_Vector()));
    alpha*=grid.one_over_dX;
    T interpolated_in_x[4],interpolated_in_y[4];
    T value[4];
    for(int c=0;c<4;c++){
        for(int b=0;b<4;b++){
            for(int a=0;a<4;a++){
                VECTOR<int,3> ind(index.x+a,index.y+b,index.z+c);
                value[a]=u(axis,ind);}
            interpolated_in_x[b]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(value[0],value[1],value[2],value[3],alpha.x);}
        interpolated_in_y[c]=cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_x[0],interpolated_in_x[1],interpolated_in_x[2],interpolated_in_x[3],alpha.y);}
    return cub.cubic_mn_interpolation.Cubic_MN_Monotonic(interpolated_in_y[0],interpolated_in_y[1],interpolated_in_y[2],interpolated_in_y[3],alpha.z);
}
//#####################################################################
// Function Base_Index
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> VECTOR<int,TV::m> CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Base_Index(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& u,const TV& X) const
{
    return INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::Clamped_Index_Interior_End_Minus_One(grid,u,X)-TV_INT::All_Ones_Vector();
}
//#####################################################################
// Function Baes_Index_Face
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> VECTOR<int,TV::m> CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
Base_Index_Face(const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,int axis,const TV& X) const
{
    return grid.Cell(X+(T).5*(TV::Axis_Vector(axis)-3)*grid.dX);
}
//#####################################################################
// Function From_Block_Face_Component
//#####################################################################
template<class TV,class T2,class T_FACE_LOOKUP> typename TV::SCALAR CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP>::
From_Block_Face_Component(const int axis,const GRID<TV>& grid,const BLOCK_UNIFORM<TV>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
{
    return From_Block_Face_Component_Helper(axis,grid,u,X,Base_Index_Face(grid,u,axis,X));
}

namespace PhysBAM{
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<float,1>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<float,2>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<float,3>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<float,2>,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<float,2>,VECTOR<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<float,3>,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<double,1>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<double,2>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<double,3>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<double,2>,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<double,2>,VECTOR<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class CUBIC_MONOTONIC_INTERPOLATION_UNIFORM<VECTOR<double,3>,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
}
