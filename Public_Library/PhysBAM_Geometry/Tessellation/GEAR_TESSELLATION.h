//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#ifndef __GEAR_TESSELLATION__
#define __GEAR_TESSELLATION__
 
namespace PhysBAM{
template<class TV> class GEAR;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGULATED_AREA;
template<class T> class SEGMENTED_CURVE_2D;
template<class T,int d> class VECTOR;

namespace TESSELLATION{
//#####################################################################
template<class T> void Boundary_Points(ARRAY<VECTOR<T,2> >& pts,const SMOOTH_GEAR<VECTOR<T,2> >& gear,int n=4);
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<T,3> >& gear,int n=4);
template<class T> TRIANGULATED_AREA<T>* Generate_Triangles(const SMOOTH_GEAR<VECTOR<T,2> >& gear,int n=4);
template<class T> TRIANGULATED_SURFACE<T>* Tessellate_Boundary(const SMOOTH_GEAR<VECTOR<T,3> >& gear,int n=4)
{return Generate_Triangles(gear,n=4);}
template<class T> SEGMENTED_CURVE_2D<T>* Tessellate_Boundary(const SMOOTH_GEAR<VECTOR<T,2> >& gear,int n=4);
//#####################################################################
}
}
#endif
