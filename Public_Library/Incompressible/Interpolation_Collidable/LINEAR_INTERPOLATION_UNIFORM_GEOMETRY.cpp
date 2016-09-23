#include <Grid_PDE/Interpolation/LINEAR_INTERPOLATION_UNIFORM_DEFINITIONS.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM.h>
#include <Incompressible/Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_UNIFORM.h>

namespace PhysBAM{
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,1>,float,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<float,1>,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,1>,float,FACE_LOOKUP_COLLIDABLE_UNIFORM<VECTOR<float,1>,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,float,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,float,FACE_LOOKUP_COLLIDABLE_UNIFORM<VECTOR<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,3>,float,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,3>,float,FACE_LOOKUP_COLLIDABLE_UNIFORM<VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,1>,double,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<double,1>,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,1>,double,FACE_LOOKUP_COLLIDABLE_UNIFORM<VECTOR<double,1>,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,double,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,double,FACE_LOOKUP_COLLIDABLE_UNIFORM<VECTOR<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,3>,double,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,3>,double,FACE_LOOKUP_COLLIDABLE_UNIFORM<VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > > >;
}
