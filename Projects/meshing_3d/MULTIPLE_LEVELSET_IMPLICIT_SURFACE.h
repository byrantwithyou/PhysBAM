//#####################################################################
// Copyright 2005, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIPLE_LEVELSET_IMPLICIT_SURFACE
//##################################################################### 
#ifndef __MULTIPLE_LEVELSET_IMPLICIT_SURFACE__
#define __MULTIPLE_LEVELSET_IMPLICIT_SURFACE__

#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>

namespace PhysBAM{

template<class T_input>
class MULTIPLE_LEVELSET_IMPLICIT_SURFACE:public IMPLICIT_OBJECT<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SYMMETRIC_MATRIX<T,3> T_SYMMETRIC_MATRIX;
public:
    using IMPLICIT_OBJECT<TV>::box;

    T inside_threshold;
    LEVELSET_IMPLICIT_OBJECT<TV>* primary_levelset;
    ARRAY<LEVELSET_IMPLICIT_OBJECT<TV>*> secondary_levelsets;
private:
    T minimum_cell_size;
public:

    MULTIPLE_LEVELSET_IMPLICIT_SURFACE()
        :primary_levelset(0)
    {}
    
    void Set_Inside_Threshold(const T inside_threshold_input=(T)1e-3)
    {inside_threshold=inside_threshold_input;}

    LEVELSET_IMPLICIT_OBJECT<TV>* Select(const TV& location) const
    {for(int i=1;i<=secondary_levelsets.m;i++)
        if(secondary_levelsets(i)->levelset.grid.Domain().Inside(location,inside_threshold) ) return secondary_levelsets(i);
    return primary_levelset;}

    TV Surface(const TV& location,const T tolerance=0,const int max_iterations=1) const 
    {return Select(location)->Surface(location,tolerance,max_iterations);}

    T_SYMMETRIC_MATRIX Hessian(const TV& X) const  PHYSBAM_OVERRIDE
    {return Select(X)->Hessian(X);}

    VECTOR<T,2> Principal_Curvatures(const TV& X) const PHYSBAM_OVERRIDE
    {return Select(X)->Principal_Curvatures(X);}

    TV Normal(const TV& location,const int aggregate=-1) const  PHYSBAM_OVERRIDE
    {return Select(location)->Normal(location,aggregate);}
    
    TV Extended_Normal(const TV& location,const int aggregate=-1) const  PHYSBAM_OVERRIDE
    {return Select(location)->Extended_Normal(location,aggregate);}

    T operator()(const TV& location) const  PHYSBAM_OVERRIDE
    {return Select(location)->operator()(location);}

    T Extended_Phi(const TV& location) const  PHYSBAM_OVERRIDE
    {return Select(location)->Extended_Phi(location);}

    void Compute_Normals() PHYSBAM_OVERRIDE
    {primary_levelset->Compute_Normals();for(int i=1;i<=secondary_levelsets.m;i++)secondary_levelsets(i)->Compute_Normals();}

    void Update_Box() PHYSBAM_OVERRIDE
    {primary_levelset->Update_Box();for(int i=1;i<=secondary_levelsets.m;i++)secondary_levelsets(i)->Update_Box();box=primary_levelset->box;}

    void Update_Minimum_Cell_Size(const int maximum_depth=0) PHYSBAM_OVERRIDE
    {primary_levelset->Update_Minimum_Cell_Size();for(int i=1;i<=secondary_levelsets.m;i++)secondary_levelsets(i)->Update_Minimum_Cell_Size();
    minimum_cell_size=primary_levelset->Minimum_Cell_Size();for(int i=1;i<=secondary_levelsets.m;i++)minimum_cell_size=min(minimum_cell_size,secondary_levelsets(i)->Minimum_Cell_Size());}

    T Minimum_Cell_Size_Within_Box(const RANGE<TV>& box)const PHYSBAM_OVERRIDE
    {T result=primary_levelset->Minimum_Cell_Size();
    for(int i=1;i<=secondary_levelsets.m;i++)if(box.Intersection(secondary_levelsets(i)->box,-inside_threshold))result=min(result,secondary_levelsets(i)->Minimum_Cell_Size());
    return result;}

    T Minimum_Cell_Size() const PHYSBAM_OVERRIDE
    {return minimum_cell_size;}

//#####################################################################
};
}
#endif

