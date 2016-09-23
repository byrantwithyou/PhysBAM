//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Vectors/VECTOR_3D.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER.h>
#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER_DEFINITIONS.h>
using namespace PhysBAM;
template<class T> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,3> >::
LINEAR_INTERPOLATION_MAC_HELPER(const T_BLOCK& block,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities)
    :base(block.Minimum_Corner()),center(block.Center()),one_over_DX(block.One_Over_DX())
{
    FACE_LOOKUP_UNIFORM<TV> face_velocities_lookup(face_velocities);
    static const int rotated_face_x[12]={0,1,2,3,4,5,6,7,8,9,10,11};
    u2=block.Face_X_Value(face_velocities_lookup,rotated_face_x[1]);u5=block.Face_X_Value(face_velocities_lookup,rotated_face_x[4]);
    u8=block.Face_X_Value(face_velocities_lookup,rotated_face_x[7]);u00=block.Face_X_Value(face_velocities_lookup,rotated_face_x[10]);
    slope_u01=one_over_DX.x*(u2-block.Face_X_Value(face_velocities_lookup,rotated_face_x[0]));
    slope_u12=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[2])-u2);
    slope_u34=one_over_DX.x*(u5-block.Face_X_Value(face_velocities_lookup,rotated_face_x[3]));
    slope_u45=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[5])-u5);
    slope_u67=one_over_DX.x*(u8-block.Face_X_Value(face_velocities_lookup,rotated_face_x[6]));
    slope_u78=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[8])-u8);
    slope_u10_11=one_over_DX.x*(u00-block.Face_X_Value(face_velocities_lookup,rotated_face_x[9]));
    slope_u11_12=one_over_DX.x*(block.Face_X_Value(face_velocities_lookup,rotated_face_x[11])-u00);
    static const int rotated_face_y[12]={0,2,4,1,3,5,6,8,10,7,9,11};
    v2=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[1]);v5=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[4]);
    v8=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[7]);v00=block.Face_Y_Value(face_velocities_lookup,rotated_face_y[10]);
    slope_v01=one_over_DX.y*(v2-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[0]));
    slope_v12=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[2])-v2);
    slope_v34=one_over_DX.y*(v5-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[3]));
    slope_v45=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[5])-v5);
    slope_v67=one_over_DX.y*(v8-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[6]));
    slope_v78=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[8])-v8);
    slope_v10_11=one_over_DX.y*(v00-block.Face_Y_Value(face_velocities_lookup,rotated_face_y[9]));
    slope_v11_12=one_over_DX.y*(block.Face_Y_Value(face_velocities_lookup,rotated_face_y[11])-v00);
    static const int rotated_face_z[12]={0,4,8,1,5,9,2,6,10,3,7,11};
    w2=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[1]);w5=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[4]);
    w8=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[7]);w00=block.Face_Z_Value(face_velocities_lookup,rotated_face_z[10]);
    slope_w01=one_over_DX.z*(w2-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[0]));
    slope_w12=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[2])-w2);
    slope_w34=one_over_DX.z*(w5-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[3]));
    slope_w45=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[5])-w5);
    slope_w67=one_over_DX.z*(w8-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[6]));
    slope_w78=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[8])-w8);
    slope_w10_11=one_over_DX.z*(w00-block.Face_Z_Value(face_velocities_lookup,rotated_face_z[9]));
    slope_w11_12=one_over_DX.z*(block.Face_Z_Value(face_velocities_lookup,rotated_face_z[11])-w00);
}
template<class T> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,3> >::
~LINEAR_INTERPOLATION_MAC_HELPER()
{
}
template<class T> VECTOR<T,3> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,3> >::
Interpolate_Face(const TV& X) const
{
    TV yxz(X.y,X.x,X.z),zxy(X.z,X.x,X.y);
    return TV(X.x<center.x?LINEAR_INTERPOLATION<T,T>::Trilinear(u2,u5,u8,u00,one_over_DX.y,one_over_DX.z,center.x,base.y,base.z,slope_u01,slope_u34,slope_u67,slope_u10_11,X)
        :LINEAR_INTERPOLATION<T,T>::Trilinear(u2,u5,u8,u00,one_over_DX.y,one_over_DX.z,center.x,base.y,base.z,slope_u12,slope_u45,slope_u78,slope_u11_12,X),
        X.y<center.y?LINEAR_INTERPOLATION<T,T>::Trilinear(v2,v5,v8,v00,one_over_DX.x,one_over_DX.z,center.y,base.x,base.z,slope_v01,slope_v34,slope_v67,slope_v10_11,yxz)
        :LINEAR_INTERPOLATION<T,T>::Trilinear(v2,v5,v8,v00,one_over_DX.x,one_over_DX.z,center.y,base.x,base.z,slope_v12,slope_v45,slope_v78,slope_v11_12,yxz),
        X.z<center.z?LINEAR_INTERPOLATION<T,T>::Trilinear(w2,w5,w8,w00,one_over_DX.x,one_over_DX.y,center.z,base.x,base.y,slope_w01,slope_w34,slope_w67,slope_w10_11,zxy)
        :LINEAR_INTERPOLATION<T,T>::Trilinear(w2,w5,w8,w00,one_over_DX.x,one_over_DX.y,center.z,base.x,base.y,slope_w12,slope_w45,slope_w78,slope_w11_12,zxy));
}
template<class T> void LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,3> >::
Block_Transfer(const T_BLOCK& source_block,const ARRAY<bool,FACE_INDEX<TV::m> >& source_values,const BLOCK_UNIFORM<TV>& destination_block,ARRAY<T,FACE_INDEX<3> >& destination_values)
{
    for(int i=0;i<GRID<TV>::number_of_faces_per_block/TV::m;i++){
        destination_block.Face_X_Reference(destination_values,i)=(T)source_block.Face_X_Value(source_values,i);
        destination_block.Face_Y_Reference(destination_values,i)=(T)source_block.Face_Y_Value(source_values,i);
        destination_block.Face_Z_Reference(destination_values,i)=(T)source_block.Face_Z_Value(source_values,i);}
}
template<class T> VECTOR<T,3> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<T,3> >::
Interpolate_Face_Normalized(const T_BLOCK& block,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const ARRAY<bool,FACE_INDEX<TV::m> >& face_velocities_valid,const TV& X,const TV& default_value)
{
    static const GRID<TV> valid_values_grid=GRID<TV>(TV_INT()+2,RANGE<TV>::Unit_Box()).Get_MAC_Grid_At_Regular_Positions();
    static const BLOCK_UNIFORM<TV> valid_values_block(valid_values_grid,VECTOR<int,3>(2,2,2));
    ARRAY<T,FACE_INDEX<3> > valid_values(valid_values_grid);Block_Transfer(block,face_velocities_valid,valid_values_block,valid_values);
    TV DX=Transformed(block,X),velocity=Interpolate_Face_Transformed(block,face_velocities,DX),weight=Interpolate_Face_Transformed(valid_values_block,valid_values,DX);
    return TV(weight.x==0?default_value.x:velocity.x/weight.x,weight.y==0?default_value.y:velocity.y/weight.y,weight.z==0?default_value.z:velocity.z/weight.z);
}
namespace PhysBAM{
template class LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >;
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<float,3> >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >(BLOCK_UNIFORM<VECTOR<float,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > const&,VECTOR<float,3> const&);
template class LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >;
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_HELPER<VECTOR<double,3> >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >(BLOCK_UNIFORM<VECTOR<double,3> > const&,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > const&,VECTOR<double,3> const&);
}
