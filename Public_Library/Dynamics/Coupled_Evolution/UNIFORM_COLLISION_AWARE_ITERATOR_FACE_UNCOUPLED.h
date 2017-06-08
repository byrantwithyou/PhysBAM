//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED
//#####################################################################
#ifndef __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED__
#define __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED__

#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Vectors/VECTOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>

namespace PhysBAM{
template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;

template<class TV>
class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED:public FACE_ITERATOR<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;
    const ARRAY<bool,TV_INT>& outside_fluid;
public:
    typedef typename TV::SCALAR T;
    typedef FACE_ITERATOR<TV> BASE;typedef typename GRID<TV>::REGION T_REGION;
    using BASE::grid;using BASE::First_Cell_Index;using BASE::Second_Cell_Index;
    using BASE::Valid;using BASE::face;

    int collision_index;
    const ARRAY<COLLISION_FACE_INFO<TV> >& collision_face_info;
    int scan_end;
    bool use_outside;

    // axis_input==0 means iterate through faces in all dimensions
    explicit UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input=0,bool use_outside_input=true);

    ~UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED();

    void Next()
    {face.index(TV::m-1)++;if(face.index(TV::m-1)>=scan_end) Next_Helper();}

    // TODO: Careful about ghost cells and the validity of outside_fluid.
    void Next_Fluid()
    {Next();if(use_outside) while(Valid() && outside_fluid(First_Cell_Index()) && outside_fluid(Second_Cell_Index())) Next();}

    void Next_Helper();
private:
    int Compare_Collision_Index() const;
};
}
#endif
