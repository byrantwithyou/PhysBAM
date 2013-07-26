//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
namespace PhysBAM{
template ARRAY<int>::~ARRAY();
template ARRAY<float>::~ARRAY();
template ARRAY<double>::~ARRAY();
template void ARRAY<int>::Resize_Helper(int, bool, bool, const int&);
template void ARRAY<float>::Resize_Helper(int, bool, bool, const float&);
template void ARRAY<double>::Resize_Helper(int, bool, bool, const double&);
}
