#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM_DEFINITIONS.h>
#include <Incompressible/Collisions_And_Interactions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>
#include <Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>

namespace PhysBAM{
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<float,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_COLLIDABLE_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<GRID<VECTOR<double,3> > > >;
}
