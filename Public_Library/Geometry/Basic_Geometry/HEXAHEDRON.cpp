//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Geometry/Basic_Geometry/HEXAHEDRON.h>
using namespace PhysBAM;
//#####################################################################
// Function Volume
//#####################################################################
template<class T> T HEXAHEDRON<T>::
Volume() const
{
    return Volume(x0,x1,x2,x3,x4,x5,x6,x7);
}
//#####################################################################
// Function Signed_Volume
//#####################################################################
template<class T> T HEXAHEDRON<T>::
Signed_Volume() const
{
    return Signed_Volume(x0,x1,x2,x3,x4,x5,x6,x7);
}
//#####################################################################
namespace PhysBAM{
template class HEXAHEDRON<float>;
template class HEXAHEDRON<double>;
}
