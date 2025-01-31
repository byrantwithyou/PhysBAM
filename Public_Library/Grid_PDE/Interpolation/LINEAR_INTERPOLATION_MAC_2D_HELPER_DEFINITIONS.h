//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/PAIR.h>
#include <Core/Math_Tools/Componentwise_Min_Max.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER.h>
using namespace PhysBAM;
template<class T> template<class T_BLOCK_2,class T_FACE_LOOKUP> T LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,2> >::
Interpolate_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX) // between 0 and 1 in the dual cell
{
    if(DX.x<.5) return LINEAR_INTERPOLATION<T,T>::Bilinear(
        block.Face_X_Value(face_velocities,0),block.Face_X_Value(face_velocities,1),
        block.Face_X_Value(face_velocities,3),block.Face_X_Value(face_velocities,4),VECTOR<T,2>(DX.x+(T).5,DX.y));
    else return LINEAR_INTERPOLATION<T,T>::Bilinear(
        block.Face_X_Value(face_velocities,1),block.Face_X_Value(face_velocities,2),
        block.Face_X_Value(face_velocities,4),block.Face_X_Value(face_velocities,5),VECTOR<T,2>(DX.x-(T).5,DX.y));
}
template<class T> template<class T_BLOCK_2,class T_FACE_LOOKUP> T LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,2> >::
Interpolate_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX) // between 0 and 1 in the dual cell
{
    if(DX.y<.5) return LINEAR_INTERPOLATION<T,T>::Bilinear(
        block.Face_Y_Value(face_velocities,0),block.Face_Y_Value(face_velocities,1),
        block.Face_Y_Value(face_velocities,2),block.Face_Y_Value(face_velocities,3),VECTOR<T,2>(DX.x,DX.y+(T).5));
    else return LINEAR_INTERPOLATION<T,T>::Bilinear(
        block.Face_Y_Value(face_velocities,2),block.Face_Y_Value(face_velocities,3),
        block.Face_Y_Value(face_velocities,4),block.Face_Y_Value(face_velocities,5),VECTOR<T,2>(DX.x,DX.y-(T).5));
}
template<class T> template<class T_BLOCK_2,class T_FACE_LOOKUP> ARRAY<PAIR<FACE_INDEX<2>,T> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,2> >::
Interpolate_Face_X_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX) // between 0 and 1 in the dual cell
{
    ARRAY<PAIR<FACE_INDEX<2>,T> > weights;
    if(DX.x<.5){VECTOR<T,2> w(DX.x+(T).5,DX.y);
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(0),(1-w.x)*(1-w.y)));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(1),w.x*(1-w.y)));
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(3),(1-w.x)*w.y));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(4),w.x*w.y));return weights;}
    else{VECTOR<T,2> w(DX.x-(T).5,DX.y);
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(1),(1-w.x)*(1-w.y)));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(2),w.x*(1-w.y)));
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(4),(1-w.x)*w.y));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_X_Index(5),w.x*w.y));return weights;}
}
template<class T> template<class T_BLOCK_2,class T_FACE_LOOKUP> ARRAY<PAIR<FACE_INDEX<2>,T> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,2> >::
Interpolate_Face_Y_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX) // between 0 and 1 in the dual cell
{
    ARRAY<PAIR<FACE_INDEX<2>,T> > weights;
    if(DX.y<.5){VECTOR<T,2> w(DX.x,DX.y+(T).5);
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(0),(1-w.x)*(1-w.y)));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(1),w.x*(1-w.y)));
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(2),(1-w.x)*w.y));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(3),w.x*w.y));return weights;}
    else{VECTOR<T,2> w(DX.x,DX.y-(T).5);
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(2),(1-w.x)*(1-w.y)));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(3),w.x*(1-w.y)));
        weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(4),(1-w.x)*w.y));weights.Append(PAIR<FACE_INDEX<2>,T>(block.Face_Y_Index(5),w.x*w.y));return weights;}
}
template<class T> template<class T_BLOCK_2,class T_FACE_LOOKUP> VECTOR<T,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,2> >::
Extrema_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,2>& DX) // between 0 and 1 in the dual cell
{
    if(DX.x<.5) return VECTOR<T,2>(
        Componentwise_Min(block.Face_X_Value(u_min,0),block.Face_X_Value(u_min,1),block.Face_X_Value(u_min,3),block.Face_X_Value(u_min,4)),
        Componentwise_Max(block.Face_X_Value(u_max,0),block.Face_X_Value(u_max,1),block.Face_X_Value(u_max,3),block.Face_X_Value(u_max,4)));
    else return VECTOR<T,2>(
        Componentwise_Min(block.Face_X_Value(u_min,1),block.Face_X_Value(u_min,2),block.Face_X_Value(u_min,4),block.Face_X_Value(u_min,5)),
        Componentwise_Max(block.Face_X_Value(u_max,1),block.Face_X_Value(u_max,2),block.Face_X_Value(u_max,4),block.Face_X_Value(u_max,5)));
}
template<class T> template<class T_BLOCK_2,class T_FACE_LOOKUP> VECTOR<T,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,2> >::
Extrema_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,2>& DX) // between 0 and 1 in the dual cell
{
    if(DX.y<.5) return VECTOR<T,2>(
        Componentwise_Min(block.Face_Y_Value(u_min,0),block.Face_Y_Value(u_min,1),block.Face_Y_Value(u_min,2),block.Face_Y_Value(u_min,3)),
        Componentwise_Max(block.Face_Y_Value(u_max,0),block.Face_Y_Value(u_max,1),block.Face_Y_Value(u_max,2),block.Face_Y_Value(u_max,3)));
    else return VECTOR<T,2>(
        Componentwise_Min(block.Face_Y_Value(u_min,2),block.Face_Y_Value(u_min,3),block.Face_Y_Value(u_min,4),block.Face_Y_Value(u_min,5)),
        Componentwise_Max(block.Face_Y_Value(u_max,2),block.Face_Y_Value(u_max,3),block.Face_Y_Value(u_max,4),block.Face_Y_Value(u_max,5)));
}
