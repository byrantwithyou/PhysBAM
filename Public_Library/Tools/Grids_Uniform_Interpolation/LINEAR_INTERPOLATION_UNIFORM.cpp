#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM_DEFINITIONS.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>

namespace PhysBAM{
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,1>,MATRIX<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,1>,SYMMETRIC_MATRIX<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,1>,VECTOR<double,1>,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,1>,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,1>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,MATRIX<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,SYMMETRIC_MATRIX<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,VECTOR<double,2>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,VECTOR<double,4>,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,2>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,3>,MATRIX<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,3>,SYMMETRIC_MATRIX<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,3>,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,3>,VECTOR<double,5>,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<double,3>,double,FACE_LOOKUP_UNIFORM<VECTOR<double,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,1>,MATRIX<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,1>,SYMMETRIC_MATRIX<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,1>,VECTOR<float,1>,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,1>,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,1>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,1> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,MATRIX<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,SYMMETRIC_MATRIX<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,VECTOR<float,2>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,VECTOR<float,4>,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,2>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,2> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,3>,MATRIX<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,3>,SYMMETRIC_MATRIX<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,3>,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,3>,VECTOR<float,5>,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
template class LINEAR_INTERPOLATION_UNIFORM<VECTOR<float,3>,float,FACE_LOOKUP_UNIFORM<VECTOR<float,3> > >;
}
