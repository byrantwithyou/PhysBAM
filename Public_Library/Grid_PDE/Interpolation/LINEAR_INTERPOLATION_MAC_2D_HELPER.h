//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_MAC_2D_HELPER
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_MAC_2D_HELPER__
#define __LINEAR_INTERPOLATION_MAC_2D_HELPER__

#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/PAIR.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Grid_Tools/Arrays/FACE_ARRAYS.h>
#include <Grid_Tools/Grids/BLOCK_UNIFORM.h>
namespace PhysBAM{

template<class TV> class LINEAR_INTERPOLATION_MAC_HELPER;

template<class T>
class LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename GRID<TV>::BLOCK T_BLOCK;
    typedef TV_INT T_INDEX;
public:
    VECTOR<T,2> base,center,one_over_DX;
    T u2,u5,slope_u01,slope_u12,slope_u34,slope_u45; // standard y-x major ordering
    T v2,v5,slope_v01,slope_v12,slope_v34,slope_v45; // x-y major ordering for symmetry with u
    
    LINEAR_INTERPOLATION_MAC_HELPER(const T_BLOCK& block,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities);

    VECTOR<T,2> Interpolate_Face(const VECTOR<T,2>& X) const;

    template<class T_FACE_LOOKUP>
    static VECTOR<T,2> Interpolate_Face(const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& X)
    {return Interpolate_Face_Transformed(block,face_velocities,Transformed(block,X));}

    template<class T_FACE_LOOKUP>
    static T Interpolate_Face_Component(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& X)
    {switch(axis){
        case 0:return Interpolate_Face_X_Transformed(block,face_velocities,Transformed(block,X));
        default:assert(axis==1);return Interpolate_Face_Y_Transformed(block,face_velocities,Transformed(block,X));}}

    template<class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<2>,T> > Interpolate_Face_Component_Weights(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& X)
    {switch(axis){
        case 0:return Interpolate_Face_X_Transformed_Weights(block,face_velocities,Transformed(block,X));
        default:assert(axis==1);return Interpolate_Face_Y_Transformed_Weights(block,face_velocities,Transformed(block,X));}}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,2> Interpolate_Face_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX)
    {return VECTOR<T,2>(Interpolate_Face_X_Transformed(block,face_velocities,DX),Interpolate_Face_Y_Transformed(block,face_velocities,DX));}

    // assumes face_velocities are 0 where not valid
    static VECTOR<T,2> Interpolate_Face_Normalized(const T_BLOCK& block,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid,const VECTOR<T,2>& X,
        const VECTOR<T,2>& default_value=TV());

    template<class T_FACE_LOOKUP>
    static VECTOR<VECTOR<T,2>,2> Extrema_Face(const T_BLOCK& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,2>& X)
    {return Extrema_Face_Transformed(block,u_min,u_max,Transformed(block,X));}

    template<class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_Component(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,2>& X)
    {switch(axis){
          case 0:return Extrema_Face_X_Transformed(block,u_min,u_max,Transformed(block,X));
          default:assert(axis==1);return Extrema_Face_Y_Transformed(block,u_min,u_max,Transformed(block,X));}}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<VECTOR<T,2>,2> Extrema_Face_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,2>& DX)
    {VECTOR<T,2> x_extrema=Extrema_Face_X_Transformed(block,u_min,u_max,DX),y_extrema=Extrema_Face_Y_Transformed(block,u_min,u_max,DX);
    return VECTOR<VECTOR<T,2>,2>(VECTOR<T,2>(x_extrema.x,y_extrema.x),VECTOR<T,2>(x_extrema.y,y_extrema.y));}

private:
    static void Block_Transfer(const T_BLOCK& source_block,const ARRAY<bool,FACE_INDEX<TV::m> >& source_values,const BLOCK_UNIFORM<TV>& destination_block,
        ARRAY<T,FACE_INDEX<2> >& destination_values)
    {for(int i=0;i<GRID<TV>::number_of_faces_per_block/TV::m;i++){
        destination_block.Face_X_Reference(destination_values,i)=(T)source_block.Face_X_Value(source_values,i);
        destination_block.Face_Y_Reference(destination_values,i)=(T)source_block.Face_Y_Value(source_values,i);}}

    static VECTOR<T,2> Transformed(const T_BLOCK& block,const VECTOR<T,2>& X)
    {return (X-block.Minimum_Corner())*block.One_Over_DX();}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static T Interpolate_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static T Interpolate_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<2>,T> > Interpolate_Face_X_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<2>,T> > Interpolate_Face_Y_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,2>& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,2>& DX); // between 0 and 1 in the dual cell
    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_Y_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,2>& DX); // between 0 and 1 in the dual cell
public:

    template<class T_ARRAYS_2>
    static typename T_ARRAYS_2::ELEMENT Interpolate_Cell(const T_BLOCK& block,const T_ARRAYS_2& cell_value,const VECTOR<T,2>& X)
    {return LINEAR_INTERPOLATION<T,typename T_ARRAYS_2::ELEMENT>::Bilinear(cell_value(block.Cell(0)),cell_value(block.Cell(1)),cell_value(block.Cell(2)),cell_value(block.Cell(3)),
        Transformed(block,X));}

//#####################################################################
};
}
#endif
