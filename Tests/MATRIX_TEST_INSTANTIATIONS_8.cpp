#include "MATRIX_TEST_DEFINITIONS.h"
#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
typedef float T;
#else
typedef double T;
#endif
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<DIAGONAL_MATRIX<T,2>,MATRIX_MXN<T> >(DIAGONAL_MATRIX<T,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<DIAGONAL_MATRIX<T,3>,MATRIX_MXN<T> >(DIAGONAL_MATRIX<T,3>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,1,1>,MATRIX_MXN<T> >(MATRIX<T,1,1>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,1,2>,MATRIX_MXN<T> >(MATRIX<T,1,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,1,3>,MATRIX_MXN<T> >(MATRIX<T,1,3>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,1,4>,MATRIX_MXN<T> >(MATRIX<T,1,4>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,1,5>,MATRIX_MXN<T> >(MATRIX<T,1,5>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,1,6>,MATRIX_MXN<T> >(MATRIX<T,1,6>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,2,1>,MATRIX_MXN<T> >(MATRIX<T,2,1>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,2,2>,MATRIX_MXN<T> >(MATRIX<T,2,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,2,3>,MATRIX_MXN<T> >(MATRIX<T,2,3>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,2,4>,MATRIX_MXN<T> >(MATRIX<T,2,4>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,2,5>,MATRIX_MXN<T> >(MATRIX<T,2,5>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,2,6>,MATRIX_MXN<T> >(MATRIX<T,2,6>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,3,1>,MATRIX_MXN<T> >(MATRIX<T,3,1>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,3,2>,MATRIX_MXN<T> >(MATRIX<T,3,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,3,3>,MATRIX_MXN<T> >(MATRIX<T,3,3>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,3,4>,MATRIX_MXN<T> >(MATRIX<T,3,4>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,3,5>,MATRIX_MXN<T> >(MATRIX<T,3,5>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,3,6>,MATRIX_MXN<T> >(MATRIX<T,3,6>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,4,1>,MATRIX_MXN<T> >(MATRIX<T,4,1>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,4,2>,MATRIX_MXN<T> >(MATRIX<T,4,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,4,3>,MATRIX_MXN<T> >(MATRIX<T,4,3>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,4,4>,MATRIX_MXN<T> >(MATRIX<T,4,4>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,4,5>,MATRIX_MXN<T> >(MATRIX<T,4,5>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,4,6>,MATRIX_MXN<T> >(MATRIX<T,4,6>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,5,1>,MATRIX_MXN<T> >(MATRIX<T,5,1>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,5,2>,MATRIX_MXN<T> >(MATRIX<T,5,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,5,3>,MATRIX_MXN<T> >(MATRIX<T,5,3>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,5,4>,MATRIX_MXN<T> >(MATRIX<T,5,4>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,5,5>,MATRIX_MXN<T> >(MATRIX<T,5,5>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,5,6>,MATRIX_MXN<T> >(MATRIX<T,5,6>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,6,1>,MATRIX_MXN<T> >(MATRIX<T,6,1>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,6,2>,MATRIX_MXN<T> >(MATRIX<T,6,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,6,3>,MATRIX_MXN<T> >(MATRIX<T,6,3>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,6,4>,MATRIX_MXN<T> >(MATRIX<T,6,4>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,6,5>,MATRIX_MXN<T> >(MATRIX<T,6,5>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX<T,6,6>,MATRIX_MXN<T> >(MATRIX<T,6,6>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<MATRIX_MXN<T>,MATRIX_MXN<T> >(MATRIX_MXN<T>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<SYMMETRIC_MATRIX<T,2>,MATRIX_MXN<T> >(SYMMETRIC_MATRIX<T,2>,MATRIX_MXN<T>,int) const;
template bool MATRIX_TESTS<T>::Arbitrary_Test_Two_Sizes<SYMMETRIC_MATRIX<T,3>,MATRIX_MXN<T> >(SYMMETRIC_MATRIX<T,3>,MATRIX_MXN<T>,int) const;
