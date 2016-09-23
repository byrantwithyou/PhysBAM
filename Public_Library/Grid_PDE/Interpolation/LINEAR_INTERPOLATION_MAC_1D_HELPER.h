//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_MAC_1D_HELPER
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_MAC_1D_HELPER__
#define __LINEAR_INTERPOLATION_MAC_1D_HELPER__

#include <Core/Log/DEBUG_UTILITIES.h>
#include <Core/Math_Tools/Componentwise_Min_Max.h>
#include <Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <Grid_Tools/Grids/BLOCK_UNIFORM.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
namespace PhysBAM{

template<class TV> class LINEAR_INTERPOLATION_MAC_HELPER;

template<class T>
class LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,1> >
{   
    typedef VECTOR<T,1> TV;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename GRID<TV>::BLOCK T_BLOCK;
    typedef TV_INT T_INDEX;
public:
    VECTOR<T,1> center;
    T u2,slope_u01,slope_u12; // standard y-x major ordering
    
    LINEAR_INTERPOLATION_MAC_HELPER(const T_BLOCK& block,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
        :center(block.Center())
    {
        FACE_LOOKUP_UNIFORM<TV> face_velocities_lookup(face_velocities);
        T one_over_dx=block.One_Over_DX().x;
        static const int rotated_face_x[3]={0,1,2};
        u2=block.Face_X_Value(face_velocities_lookup,rotated_face_x[1]);
        slope_u01=one_over_dx*(u2-block.Face_X_Value(face_velocities_lookup,rotated_face_x[0]));
        slope_u12=one_over_dx*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[2])-u2);
    }

    VECTOR<T,1> Interpolate_Face(const VECTOR<T,1>& X) const
    {return VECTOR<T,1>(X.x<center.x?LINEAR_INTERPOLATION<T,T>::Linear(center.x,u2,slope_u01,X)
                                     :LINEAR_INTERPOLATION<T,T>::Linear(center.x,u2,slope_u12,X));}

    template<class T_FACE_LOOKUP>
    static VECTOR<T,1> Interpolate_Face(const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,1>& X)
    {return Interpolate_Face_Transformed(block,face_velocities,Transformed(block,X));}

    template<class T_FACE_LOOKUP>
    static T Interpolate_Face_Component(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,1>& X)
    {assert(axis==0);return Interpolate_Face_X_Transformed(block,face_velocities,Transformed(block,X));}

    template<class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<1>,T> > Interpolate_Face_Component_Weights(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,1>& X)
    {assert(axis==0);return Interpolate_Face_X_Transformed_Weights(block,face_velocities,Transformed(block,X));}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,1> Interpolate_Face_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,1>& DX)
    {return VECTOR<T,1>(Interpolate_Face_X_Transformed(block,face_velocities,DX));}

    // assumes face_velocities are 0 where not valid
    static VECTOR<T,1> Interpolate_Face_Normalized(const T_BLOCK& block,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid,const VECTOR<T,1>& X,
         const VECTOR<T,1>& default_value=TV())
    {static const GRID<TV> valid_values_grid=GRID<TV>(TV_INT(2),RANGE<TV>::Unit_Box()).Get_MAC_Grid_At_Regular_Positions();
    static const BLOCK_UNIFORM<TV> valid_values_block(valid_values_grid,VECTOR<int,1>(2));
    ARRAY<T,FACE_INDEX<1> > valid_values(valid_values_grid);Block_Transfer(block,face_velocities_valid,valid_values_block,valid_values);
    VECTOR<T,1> DX=Transformed(block,X),velocity=Interpolate_Face_Transformed(block,face_velocities,DX),weight=Interpolate_Face_Transformed(valid_values_block,valid_values,DX);
    return VECTOR<T,1>(weight.x==0?default_value.x:velocity.x/weight.x);}

    template<class T_FACE_LOOKUP>
    static VECTOR<VECTOR<T,1>,2> Extrema_Face(const T_BLOCK& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,1>& X)
    {return Extrema_Face_Transformed(block,u_min,u_max,Transformed(block,X));}

    template<class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_Component(const int axis,const T_BLOCK& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,1>& X)
    {assert(axis==0);return Extrema_Face_X_Transformed(block,u_min,u_max,Transformed(block,X));}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<VECTOR<T,1>,2> Extrema_Face_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,1>& DX)
    {VECTOR<T,2> x_extrema=Extrema_Face_X_Transformed(block,u_min,u_max,DX);return VECTOR<VECTOR<T,1>,2>(VECTOR<T,1>(x_extrema.x),VECTOR<T,1>(x_extrema.y));}

private:
    static void Block_Transfer(const T_BLOCK& source_block,const ARRAY<bool,FACE_INDEX<TV::m> >& source_values,const BLOCK_UNIFORM<TV>& destination_block,ARRAY<T,FACE_INDEX<1> >& destination_values)
    {for(int i=0;i<GRID<TV>::number_of_faces_per_block/TV::m;i++)
        destination_block.Face_X_Reference(destination_values,i)=(T)source_block.Face_X_Value(source_values,i);}
        
    static VECTOR<T,1> Transformed(const T_BLOCK& block,const VECTOR<T,1>& X)
    {return (X-block.Minimum_Corner())*block.One_Over_DX();}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static T Interpolate_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,1>& DX) // between 0 and 1 in the dual cell
    {if(DX.x<.5) return LINEAR_INTERPOLATION<T,T>::Linear(block.Face_X_Value(face_velocities,0),block.Face_X_Value(face_velocities,1),VECTOR<T,1>(DX.x+(T).5));
    else return LINEAR_INTERPOLATION<T,T>::Linear(block.Face_X_Value(face_velocities,1),block.Face_X_Value(face_velocities,2),VECTOR<T,1>(DX.x-(T).5));}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static ARRAY<PAIR<FACE_INDEX<1>,T> > Interpolate_Face_X_Transformed_Weights(const T_BLOCK_2& block,const T_FACE_LOOKUP& face_velocities,const VECTOR<T,1>& DX) // between 0 and 1 in the dual cell
    {ARRAY<PAIR<FACE_INDEX<1>,T> > weights;
    if(DX.x<.5){T w=DX.x+(T).5;weights.Append(PAIR<FACE_INDEX<1>,T>(block.Face_X_Index(0),1-w));weights.Append(PAIR<FACE_INDEX<1>,T>(block.Face_X_Index(1),w));}
    else{T w=DX.x-(T).5;weights.Append(PAIR<FACE_INDEX<1>,T>(block.Face_X_Index(1),1-w));weights.Append(PAIR<FACE_INDEX<1>,T>(block.Face_X_Index(2),w));}return weights;}

    template<class T_BLOCK_2,class T_FACE_LOOKUP>
    static VECTOR<T,2> Extrema_Face_X_Transformed(const T_BLOCK_2& block,const T_FACE_LOOKUP& u_min,const T_FACE_LOOKUP& u_max,const VECTOR<T,1>& DX) // between 0 and 1 in the dual cell
    {if(DX.x<.5) return VECTOR<T,2>(Componentwise_Min(block.Face_X_Value(u_min,0),block.Face_X_Value(u_min,1)),Componentwise_Max(block.Face_X_Value(u_max,0),block.Face_X_Value(u_max,1)));
    else return VECTOR<T,2>(Componentwise_Min(block.Face_X_Value(u_min,1),block.Face_X_Value(u_min,2)),Componentwise_Max(block.Face_X_Value(u_max,1),block.Face_X_Value(u_max,2)));}
public:

    template<class T_ARRAYS_2>
    static typename T_ARRAYS_2::ELEMENT Interpolate_Cell(const T_BLOCK& block,const T_ARRAYS_2& cell_value,const VECTOR<T,1>& X)
    {return LINEAR_INTERPOLATION<T,typename T_ARRAYS_2::ELEMENT>::Linear(cell_value(block.Cell(0)),cell_value(block.Cell(1)),Transformed(block,X));}

//#####################################################################
};
}
#endif
