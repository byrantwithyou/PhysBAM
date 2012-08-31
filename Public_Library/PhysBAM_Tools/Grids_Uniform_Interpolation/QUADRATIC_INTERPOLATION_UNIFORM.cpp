//#####################################################################
// Copyright 2010, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/QUADRATIC_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
QUADRATIC_INTERPOLATION_UNIFORM()
    :ghost_cells(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
~QUADRATIC_INTERPOLATION_UNIFORM()
{
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> typename T_GRID::VECTOR_INT QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Base_Index(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    TV_INT index(rint((X-grid.domain.min_corner)*grid.one_over_dX-grid.MAC_offset));
    return clamp(index,u.domain.min_corner+1,u.domain.max_corner-2);
}
//#####################################################################
// function Clamped_To_Array
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    return From_Base_Node(grid,u,X,Base_Index(grid,u,X));
}
//#####################################################################
// Function Clamped_To_Array_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
{
    return From_Base_Node_Weights(grid,u,X,Base_Index(grid,u,X));
}
//#####################################################################
// Function Quadratic_Interpolation
//#####################################################################
template<class T> static T Quadratic_Interpolation(const T* x,T w)
{
    T a=(T).5*(x[0]+x[2])-x[1],b=(T).5*(x[2]-x[0]),c=x[1];
    return (a*w+b)*w+c;
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    return Quadratic_Interpolation(&u(index-1),w.x);
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    const T2* b=&u(index-1);
    T2 x[3];
    for(int i=0;i<3;i++,b+=u.counts.y)
        x[i]=Quadratic_Interpolation(b,w.y);
    return Quadratic_Interpolation(x,w.x);
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const
{
    TV w=(X-grid.X(index))*grid.one_over_dX;
    int dx=u.counts.y*u.counts.z;
    const T2* b=&u(index-1),*c=b;
    T2 x[3],y[3];
    for(int i=0;i<3;i++,b+=dx,c=b){
        for(int j=0;j<3;j++,c+=u.counts.z)
            y[j]=Quadratic_Interpolation(c,w.z);
        x[i]=Quadratic_Interpolation(y,w.y);}
    return Quadratic_Interpolation(x,w.x);
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,1> >& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return ARRAY<PAIR<TV_INT,T> >();
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,2> >& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return ARRAY<PAIR<TV_INT,T> >();
}
//#####################################################################
// Function From_Base_Node_Weights
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> ARRAY<PAIR<typename T_GRID::VECTOR_INT,typename T_GRID::VECTOR_T::SCALAR> > QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Base_Node_Weights_Helper(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,VECTOR<int,3> >& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return ARRAY<PAIR<TV_INT,T> >();
}
//#####################################################################
// Function Clamped_To_Array
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> typename T_GRID::VECTOR_INT QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
Base_Index_Face(const T_GRID& grid,const typename T_FACE_LOOKUP::LOOKUP& u,int axis,const TV& X) const
{
    TV offset;
    offset(axis)=(T).5;
    TV_INT I(floor((X-grid.domain.min_corner)*grid.one_over_dX+offset));
    RANGE<TV_INT> range=grid.Domain_Indices(ghost_cells);
    range.max_corner(axis)++;
    return clamp(I,range.min_corner+1,range.max_corner-2);
}
//#####################################################################
// Function From_Block_Face_Component
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> typename T_GRID::VECTOR_T::SCALAR QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
{
    return From_Block_Face_Component_Helper(axis,grid,u,X,Base_Index_Face(grid,u,axis,X));
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const VECTOR<T,1>& X,const VECTOR<int,1>& index) const
{
    TV w=(X-grid.Axis_X_Face(FACE_INDEX<TV::m>(axis,index)))*grid.one_over_dX;
    T2 x[3]={u(axis,index-1),u(axis,index),u(axis,index+1)};
    return Quadratic_Interpolation(x,w.x);
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const VECTOR<T,2>& X,const VECTOR<int,2>& index) const
{
    TV w=(X-grid.Axis_X_Face(FACE_INDEX<TV::m>(axis,index)))*grid.one_over_dX;
    TV_INT b(index-1);
    T2 x[3];
    for(int i=0;i<3;i++,b.x++){
        T2 y[3]={u(axis,b),u(axis,b+TV_INT(0,1)),u(axis,b+TV_INT(0,2))};
        x[i]=Quadratic_Interpolation(y,w.y);}
    return Quadratic_Interpolation(x,w.x);
}
//#####################################################################
// Function From_Base_Node
//#####################################################################
template<class T_GRID,class T2,class T_FACE_LOOKUP> T2 QUADRATIC_INTERPOLATION_UNIFORM<T_GRID,T2,T_FACE_LOOKUP>::
From_Block_Face_Component_Helper(const int axis,const GRID<TV>& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const VECTOR<T,3>& X,const VECTOR<int,3>& index) const
{
    TV w=(X-grid.Axis_X_Face(FACE_INDEX<TV::m>(axis,index)))*grid.one_over_dX;
    TV_INT b(index-1),c=b;
    T2 x[3],y[3];
    for(int i=0;i<3;i++,b.x++,c=b){
        for(int j=0;j<3;j++,c.y++){
            T2 z[3]={u(axis,c),u(axis,c+TV_INT(0,0,1)),u(axis,c+TV_INT(0,0,2))};
            y[j]=Quadratic_Interpolation(z,w.z);}
        x[i]=Quadratic_Interpolation(y,w.y);}
    return Quadratic_Interpolation(x,w.x);
}
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template float QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >::Clamped_To_Array(GRID<VECTOR<float,1> > const&,ARRAYS_ND_BASE<float,VECTOR<int,1> > const&,VECTOR<float,1> const&) const;
template float QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::Clamped_To_Array(GRID<VECTOR<float,2> > const&,ARRAYS_ND_BASE<float,VECTOR<int,2> > const&,VECTOR<float,2> const&) const;
template float QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::Clamped_To_Array(GRID<VECTOR<float,3> > const&,ARRAYS_ND_BASE<float,VECTOR<int,3> > const&,VECTOR<float,3> const&) const;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::QUADRATIC_INTERPOLATION_UNIFORM();
template QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::~QUADRATIC_INTERPOLATION_UNIFORM();
template double QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::Clamped_To_Array(GRID<VECTOR<double,2> > const&,ARRAYS_ND_BASE<double,VECTOR<int,2> > const&,VECTOR<double,2> const&) const;
template double QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::Clamped_To_Array(GRID<VECTOR<double,3> > const&,ARRAYS_ND_BASE<double,VECTOR<int,3> > const&,VECTOR<double,3> const&) const;
#endif
