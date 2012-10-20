#include "MATRIX_TEST_DEFINITIONS.h"
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
typedef float T;
#else
typedef double T;
#endif
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<DIAGONAL_MATRIX<T,2> >(DIAGONAL_MATRIX<T,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<DIAGONAL_MATRIX<T,3> >(DIAGONAL_MATRIX<T,3>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,1,1> >(MATRIX<T,1,1>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,1,2> >(MATRIX<T,1,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,1,3> >(MATRIX<T,1,3>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,1,4> >(MATRIX<T,1,4>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,1,5> >(MATRIX<T,1,5>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,1,6> >(MATRIX<T,1,6>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,2,1> >(MATRIX<T,2,1>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,2,2> >(MATRIX<T,2,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,2,3> >(MATRIX<T,2,3>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,2,4> >(MATRIX<T,2,4>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,2,5> >(MATRIX<T,2,5>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,2,6> >(MATRIX<T,2,6>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,3,1> >(MATRIX<T,3,1>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,3,2> >(MATRIX<T,3,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,3,3> >(MATRIX<T,3,3>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,3,4> >(MATRIX<T,3,4>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,3,5> >(MATRIX<T,3,5>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,3,6> >(MATRIX<T,3,6>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,4,1> >(MATRIX<T,4,1>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,4,2> >(MATRIX<T,4,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,4,3> >(MATRIX<T,4,3>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,4,4> >(MATRIX<T,4,4>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,4,5> >(MATRIX<T,4,5>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,4,6> >(MATRIX<T,4,6>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,5,1> >(MATRIX<T,5,1>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,5,2> >(MATRIX<T,5,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,5,3> >(MATRIX<T,5,3>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,5,4> >(MATRIX<T,5,4>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,5,5> >(MATRIX<T,5,5>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,5,6> >(MATRIX<T,5,6>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,6,1> >(MATRIX<T,6,1>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,6,2> >(MATRIX<T,6,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,6,3> >(MATRIX<T,6,3>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,6,4> >(MATRIX<T,6,4>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,6,5> >(MATRIX<T,6,5>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX<T,6,6> >(MATRIX<T,6,6>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<MATRIX_MXN<T> >(MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<SYMMETRIC_MATRIX<T,2> >(SYMMETRIC_MATRIX<T,2>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_One_Size<SYMMETRIC_MATRIX<T,3> >(SYMMETRIC_MATRIX<T,3>,int) const;
template bool MATRIX_TESTS<T>::Test<MATRIX_MXN<T> >(bool,std::string const&,bool&,MATRIX_MXN<T> const&) const;
template bool MATRIX_TESTS<T>::Test<MATRIX_MXN<T>,MATRIX_MXN<T> >(bool,std::string const&,bool&,MATRIX_MXN<T> const&,MATRIX_MXN<T> const&) const;
template bool MATRIX_TESTS<T>::Test<MATRIX_MXN<T>,MATRIX_MXN<T>,MATRIX_MXN<T> >(bool,std::string const&,bool&,MATRIX_MXN<T> const&,MATRIX_MXN<T> const&,MATRIX_MXN<T> const&) const;
