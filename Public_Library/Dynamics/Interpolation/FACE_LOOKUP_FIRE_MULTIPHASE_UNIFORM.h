//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM__
#define __FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM__

#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Dynamics/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <cassert>
namespace PhysBAM{

template<class TV> class LEVELSET_MULTIPLE;

template<class TV>
class FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;typedef PROJECTION_DYNAMICS_UNIFORM<TV> T_PROJECTION;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS;typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
public:
    typedef T ELEMENT;
        
    const T_FACE_ARRAYS& V_face;
protected:
    const T_PROJECTION& projection;
    const LEVELSET_MULTIPLE<TV>* levelset_multiple;
public:

    FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM(const T_FACE_ARRAYS& V_face_input,const T_PROJECTION& projection_input,const LEVELSET_MULTIPLE<TV>* levelset_multiple_input=0)
        :V_face(V_face_input),projection(projection_input),levelset_multiple(levelset_multiple_input)
    {}

    const T_FACE_ARRAYS& Raw_Data() const
    {return V_face;}

    class LOOKUP
    {
    private:
        const FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<TV>& face_lookup;
        const int reference_region;
    public:
        typedef T ELEMENT;
        const T_FACE_ARRAYS& V_face;

        LOOKUP(const FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM<TV>& face_lookup_input,const int reference_region_input)
            :face_lookup(face_lookup_input),reference_region(reference_region_input),V_face(face_lookup.V_face)
        {}

        int Number_Of_Ghost_Cells() const
        {return V_face.Number_Of_Ghost_Cells();}

        void Set_Reference_Point(const TV& reference_point) const
        {}
    
        T operator()(const FACE_INDEX<TV::dimension>& face) const
        {return face_lookup.projection.Face_Velocity_With_Ghost_Value_Multiphase(face_lookup.V_face,face.axis,face.index,reference_region);}
        
        T operator()(const int axis,const TV_INT& face) const
        {return face_lookup.projection.Face_Velocity_With_Ghost_Value_Multiphase(face_lookup.V_face,axis,face,reference_region);}
    };

    LOOKUP Starting_Point_Face(const int axis,const TV_INT& face) const
    {assert(levelset_multiple);return LOOKUP(*this,levelset_multiple->Inside_Region_Face(axis,face));}
    
    LOOKUP Starting_Point_Cell(const TV_INT& cell) const
    {assert(levelset_multiple);return LOOKUP(*this,levelset_multiple->Inside_Region(cell));}
        
    LOOKUP Region(const int region) const
    {return LOOKUP(*this,region);}
    
//#####################################################################
};
}
#endif
