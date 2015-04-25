//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//###################################################################
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/PRIMITIVE_MATRICES.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Tensors/TENSOR.h>
#include <Tools/Utilities/PROCESS_UTILITIES.h>
#include "TENSOR_TESTS.h"
#include "TENSOR_TESTS_DEFINITIONS.h"
using namespace PhysBAM;
template void TENSOR_TESTS<double>::Test_3<3,1,2>();
template void TENSOR_TESTS<double>::Test_3<3,1,3>();
template void TENSOR_TESTS<double>::Test_3<3,2,1>();
template void TENSOR_TESTS<double>::Test_3<3,2,2>();
