//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function maxmag  
//#####################################################################
//
// finds the maximum value in magnitude and returns it with the sign
//               
//#####################################################################
#ifndef __maxmag__
#define __maxmag__    

namespace PhysBAM{

template<class T>
inline T maxmag(const T a,const T b)
{return abs(a)>abs(b)?a:b;}

template<class T,class ...Args>
inline T maxmag(const T a,const T b,Args&&... c)
{return maxmag(a,maxmag(b,c...));}

}
#endif

