//#####################################################################
// Copyright 2009-2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER_DEFINITIONS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER_DEFINITIONS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
using namespace PhysBAM;
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,VECTOR<float,2> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,VECTOR<float,2> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > >::LOOKUP const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,VECTOR<double,2> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,VECTOR<double,2> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,VECTOR<double,2> const&);
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,VECTOR<double,2> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > >::LOOKUP const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > >::LOOKUP const&,VECTOR<double,3> const&);
