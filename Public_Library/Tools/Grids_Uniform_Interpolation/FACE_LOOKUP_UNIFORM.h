//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_LOOKUP_UNIFORM
//#####################################################################
#ifndef __FACE_LOOKUP_UNIFORM__
#define __FACE_LOOKUP_UNIFORM__

#include <Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
namespace PhysBAM{

template<class TV>
class FACE_LOOKUP_UNIFORM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef T ELEMENT;
        
    const ARRAY<T,FACE_INDEX<TV::m> >& V_face;

    FACE_LOOKUP_UNIFORM(const ARRAY<T,FACE_INDEX<TV::m> >& V_face_input)
        :V_face(V_face_input)
    {}

    const ARRAY<T,FACE_INDEX<TV::m> >& Raw_Data() const
    {return V_face;}

    int Number_Of_Ghost_Cells() const
    {return V_face.Number_Of_Ghost_Cells();}
    
    typedef FACE_LOOKUP_UNIFORM LOOKUP;

    const LOOKUP& Starting_Point_Face(const int axis,const TV_INT& face) const
    {return *this;}
    
    const LOOKUP& Starting_Point_Cell(const TV_INT& cell) const
    {return *this;}

    void Set_Reference_Point(const TV& reference_point) const
    {}

    T operator()(const int axis,const TV_INT& face) const
    {return V_face.Component(axis)(face);}

    T operator()(const FACE_INDEX<TV::dimension>& face) const
    {return V_face(face);}

//#####################################################################
};
}
#endif
