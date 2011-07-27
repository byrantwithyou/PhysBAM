//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FIRE_INTERPOLATION_POLICY__
#define __FIRE_INTERPOLATION_POLICY__

#include <PhysBAM_Dynamics/Interpolation/FIRE_INTERPOLATION_FORWARD.h>
namespace PhysBAM{

template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;
template<class T,int d> class VECTOR;

template<class T_GRID>
struct FIRE_INTERPOLATION_POLICY
{
private:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;
public:
    // multiphase fire
    typedef FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<T_GRID> FACE_LOOKUP_FIRE_MULTIPHASE;
    typedef AVERAGING_UNIFORM<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> AVERAGING_FIRE_MULTIPHASE;
    typedef LINEAR_INTERPOLATION_UNIFORM<T_GRID,T,FACE_LOOKUP_FIRE_MULTIPHASE> LINEAR_INTERPOLATION_FIRE_MULTIPHASE;
    // multiphase fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE;
};

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T>
struct FIRE_INTERPOLATION_POLICY<OCTREE_GRID<T> >
{
private:
    typedef OCTREE_GRID<T> T_GRID;
public:
    typedef FACE_LOOKUP_FIRE_DYADIC<T_GRID> FACE_LOOKUP_FIRE;
    typedef AVERAGING_DYADIC<T_GRID,FACE_LOOKUP_FIRE> AVERAGING_FIRE;
    typedef LINEAR_INTERPOLATION_DYADIC<T_GRID,T,FACE_LOOKUP_FIRE> LINEAR_INTERPOLATION_FIRE;
    // fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID,FACE_LOOKUP_FIRE> FACE_LOOKUP_FIRE_COLLIDABLE;
    // multiphase fire
    typedef FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC<T_GRID> FACE_LOOKUP_FIRE_MULTIPHASE;
    typedef AVERAGING_DYADIC<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> AVERAGING_FIRE_MULTIPHASE;
    typedef LINEAR_INTERPOLATION_DYADIC<T_GRID,T,FACE_LOOKUP_FIRE_MULTIPHASE> LINEAR_INTERPOLATION_FIRE_MULTIPHASE;
    // multiphase fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE;
};

template<class T>
struct FIRE_INTERPOLATION_POLICY<QUADTREE_GRID<T> >
{
private:
    typedef QUADTREE_GRID<T> T_GRID;
public:
    typedef FACE_LOOKUP_FIRE_DYADIC<T_GRID> FACE_LOOKUP_FIRE;
    typedef AVERAGING_DYADIC<T_GRID,FACE_LOOKUP_FIRE> AVERAGING_FIRE;
    typedef LINEAR_INTERPOLATION_DYADIC<T_GRID,T,FACE_LOOKUP_FIRE> LINEAR_INTERPOLATION_FIRE;
    // fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID,FACE_LOOKUP_FIRE> FACE_LOOKUP_FIRE_COLLIDABLE;
    // multiphase fire
    typedef FACE_LOOKUP_FIRE_MULTIPHASE_DYADIC<T_GRID> FACE_LOOKUP_FIRE_MULTIPHASE;
    typedef AVERAGING_DYADIC<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> AVERAGING_FIRE_MULTIPHASE;
    typedef LINEAR_INTERPOLATION_DYADIC<T_GRID,T,FACE_LOOKUP_FIRE_MULTIPHASE> LINEAR_INTERPOLATION_FIRE_MULTIPHASE;
    // multiphase fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_DYADIC<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE;
};
#endif

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T>
struct FIRE_INTERPOLATION_POLICY<RLE_GRID_2D<T> >
{
private:
    typedef RLE_GRID_2D<T> T_GRID;
    typedef VECTOR<T,2> TV;
public:
    // fire
    typedef FACE_LOOKUP_FIRE_RLE<T_GRID> FACE_LOOKUP_FIRE;
    typedef AVERAGING_RLE<T_GRID,FACE_LOOKUP_FIRE> AVERAGING_FIRE;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T,FACE_LOOKUP_FIRE> LINEAR_INTERPOLATION_FIRE;
    // fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID,FACE_LOOKUP_FIRE> FACE_LOOKUP_FIRE_COLLIDABLE;
    // multiphase fire
    typedef FACE_LOOKUP_FIRE_MULTIPHASE_RLE<T_GRID> FACE_LOOKUP_FIRE_MULTIPHASE;
    typedef AVERAGING_RLE<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> AVERAGING_FIRE_MULTIPHASE;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T,FACE_LOOKUP_FIRE_MULTIPHASE> LINEAR_INTERPOLATION_FIRE_MULTIPHASE;
    // multiphase fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE;
};

template<class T>
struct FIRE_INTERPOLATION_POLICY<RLE_GRID_3D<T> >
{
private:
    typedef RLE_GRID_3D<T> T_GRID;
    typedef VECTOR<T,3> TV;
public:
    // fire
    typedef FACE_LOOKUP_FIRE_RLE<T_GRID> FACE_LOOKUP_FIRE;
    typedef AVERAGING_RLE<T_GRID,FACE_LOOKUP_FIRE> AVERAGING_FIRE;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T,FACE_LOOKUP_FIRE> LINEAR_INTERPOLATION_FIRE;
    // fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID,FACE_LOOKUP_FIRE> FACE_LOOKUP_FIRE_COLLIDABLE;
    // multiphase fire
    typedef FACE_LOOKUP_FIRE_MULTIPHASE_RLE<T_GRID> FACE_LOOKUP_FIRE_MULTIPHASE;
    typedef AVERAGING_RLE<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> AVERAGING_FIRE_MULTIPHASE;
    typedef LINEAR_INTERPOLATION_RLE<T_GRID,T,FACE_LOOKUP_FIRE_MULTIPHASE> LINEAR_INTERPOLATION_FIRE_MULTIPHASE;
    // multiphase fire collidable
    typedef FACE_LOOKUP_COLLIDABLE_RLE<T_GRID,FACE_LOOKUP_FIRE_MULTIPHASE> FACE_LOOKUP_FIRE_MULTIPHASE_COLLIDABLE;
};
#endif

//#####################################################################
}
#endif
