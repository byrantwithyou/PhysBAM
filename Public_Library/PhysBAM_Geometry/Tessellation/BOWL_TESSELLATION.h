//#####################################################################
// Copyright 2011, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __BOWL_TESSELLATION__
#define __BOWL_TESSELLATION__
 
namespace PhysBAM{
template<class T> class BOWL;
template<class T> class TRIANGULATED_SURFACE;

namespace TESSELLATION{
//#####################################################################
    template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const BOWL<T>& ring,const int n=40);
//#####################################################################
}
}
#endif
