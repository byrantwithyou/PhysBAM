//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PARTICLE_GRID_FACE_ITERATOR__
#define __PARTICLE_GRID_FACE_ITERATOR__
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_INDEX.h>
namespace PhysBAM{

template<class TV> class PARTICLE_GRID_WEIGHTS;

template<class TV>
class PARTICLE_GRID_FACE_ITERATOR
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef VECTOR<typename PARTICLE_GRID_WEIGHTS<TV>::SCRATCH,TV::m> SCRATCH;
    SCRATCH& scratch;
    int axis;
    int i;

    PARTICLE_GRID_FACE_ITERATOR(const VECTOR<PARTICLE_GRID_WEIGHTS<TV>*,TV::m> weights,
        int p,bool want_gradient,SCRATCH& scratch);
    FACE_INDEX<TV::m> Index() const {return FACE_INDEX<TV::m>(axis,scratch(axis).index(i));}
    const T& Weight() const {return scratch(axis).weight(i);}
    const TV& Gradient() const {return scratch(axis).gradient(i);}

    void Next() {i++;if(i>=scratch(axis).index.m){axis++;i=0;}}
    bool Valid() const {return axis<TV::m;}
};
}
#endif
