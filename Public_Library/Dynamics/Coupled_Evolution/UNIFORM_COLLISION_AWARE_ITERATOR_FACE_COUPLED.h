//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED
//#####################################################################
#ifndef __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED__
#define __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED__

#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/SIDED_FACE_INDEX.h>

namespace PhysBAM{

template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;
template<class TV> struct COLLISION_FACE_INFO;

template<class TV>
class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED
{
public:
    enum WORKAROUND {d=TV::m};
    typedef typename TV::SCALAR T;typedef VECTOR<int,d> TV_INT;typedef typename GRID<TV>::REGION T_REGION;typedef FACE_ITERATOR<TV> BASE;

    const GRID<TV>& grid;
    int collision_index;
    int side;
    FACE_INDEX<d> face;
    const ARRAY<COLLISION_FACE_INFO<TV> >& collision_face_info;

    // TODO: Handle ghost cells, regions, side, and axis properly construction values properly.
    // axis_input==0 means iterate through faces in all dimensions
    explicit UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info);

    void Next()
    {if(++collision_index<collision_face_info.m){face=collision_face_info(collision_index).face;side=collision_face_info(collision_index).side;}}

    bool Valid()
    {return collision_index<collision_face_info.m;}

    const ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >& Get_Simplices() const
    {return collision_face_info(collision_index).simplices;}

    TV_INT First_Cell_Index() const
    {TV_INT i(face.index);i(face.axis)--;return i;}

    TV_INT Second_Cell_Index() const
    {return face.index;}

    TV_INT Real_Cell_Index() const
    {return side==0?First_Cell_Index():Second_Cell_Index();}

    TV_INT Ghost_Cell_Index() const
    {return side==1?First_Cell_Index():Second_Cell_Index();}

    RANGE<TV> Dual_Cell() const
    {TV X=grid.Face(face);return RANGE<TV>(X-(T).5*grid.dX,X+(T).5*grid.dX);}

    TV Location() const
    {return grid.Face(face);}
    
    //void Neighbor_Face_Stencil(VECTOR<SIDED_FACE_INDEX<d>,2*d-2>& neighbor_faces) const
    //{TV_INT real_cell=Real_Cell_Index();
    //for(int i=0;i<d-1;i++){int a=i+(i>=axis);
    //    neighbor_faces(2*i-1)=SIDED_FACE_INDEX<d>(2,a,real_cell);
    //    neighbor_faces(2*i)=SIDED_FACE_INDEX<d>(1,a,real_cell+TV_INT::Axis_Vector(a));}}
};
}
#endif
