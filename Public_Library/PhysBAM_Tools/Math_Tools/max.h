//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Sergey Koltakov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function max  
//#####################################################################
#ifndef __max__
#define __max__

namespace PhysBAM{
#define QQ(x) LOG::cout<<#x<<" "<<(x)<<std::endl; // INDEXING: these macros should go away when we are done with them.
#define ZZ(x) LOG::cout<<"Z "<<#x<<" "<<(x+1)<<std::endl;
#define RA(x) LOG::cout<<"R "<<#x<<" ("<<((x).min_corner+1)<<" "<<(x).max_corner<<")"<<std::endl;
template<class T>
inline T max(const T a, const T b)
{return (a>b)?a:b;}

template<class T>
inline T max(const T a,const T b,const T c)
{return max(a,max(b,c));}

template<class T>
inline T max(const T a,const T b,const T c,const T d)
{return max(a,max(b,c,d));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e)
{return max(a,max(b,c,d,e));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f)
{return max(a,max(b,c,d,e,f));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f,const T g)
{return max(a,max(b,c,d,e,f,g));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f,const T g,const T h)
{return max(a,max(b,c,d,e,f,g,h));}

template<class T>
inline T max(const T a,const T b,const T c,const T d,const T e,const T f,const T g,const T h,const T i)
{return max(a,max(b,c,d,e,f,g,h,i));}

}
#endif
