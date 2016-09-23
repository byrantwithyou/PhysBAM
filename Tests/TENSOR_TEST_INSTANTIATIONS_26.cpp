//#####################################################################
// Copyright 2015.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//###################################################################
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Tensors/TENSOR.h>
#include "TENSOR_TESTS.h"
#include "TENSOR_TESTS_DEFINITIONS.h"
using namespace PhysBAM;
template void TENSOR_TESTS<double>::Test_4<3,1,3,1>();
template void TENSOR_TESTS<double>::Test_4<3,1,3,2>();
template void TENSOR_TESTS<double>::Test_4<3,1,3,3>();
template void TENSOR_TESTS<double>::Test_4<3,2,1,1>();
