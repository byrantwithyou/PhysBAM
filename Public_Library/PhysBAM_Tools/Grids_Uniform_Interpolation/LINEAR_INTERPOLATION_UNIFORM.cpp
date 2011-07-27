#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM_DEFINITIONS.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>

template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,SYMMETRIC_MATRIX<float,2>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,1>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,SYMMETRIC_MATRIX<float,2>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,2>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,4>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,SYMMETRIC_MATRIX<float,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,5>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,SYMMETRIC_MATRIX<double,2>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,1>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,SYMMETRIC_MATRIX<double,2>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,2>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,4>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,SYMMETRIC_MATRIX<double,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,5>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >;
template class LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >;
#endif

