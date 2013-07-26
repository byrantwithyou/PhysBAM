//#####################################################################
// Copyright 2004-2006, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KINEMATIC_COLLISION_BODY
//#####################################################################
#ifndef __KINEMATIC_COLLISION_BODY__
#define __KINEMATIC_COLLISION_BODY__

#include <Tools/Grids_Uniform/GRID.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
namespace PhysBAM{

template<class TV>
class KINEMATIC_COLLISION_BODY:public RIGID_BODY<TV>
{
public:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;

    typedef RIGID_BODY<TV> BASE;

    using BASE::implicit_object;using BASE::particle_index;using BASE::Object_Space_Point;using BASE::Add_Structure;

    GRID<TV>* velocity_grid;
    ARRAY<TV,TV_INT>* velocity_field; // a field defined on velocity_grid, in object space

    LINEAR_INTERPOLATION_UNIFORM<TV,TV> interpolation;

    KINEMATIC_COLLISION_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,bool create_collision_geometry,GRID<TV> *initial_velocity_grid=0,ARRAY<TV,TV_INT> *initial_velocity_field=0)
        :BASE(rigid_body_collection,create_collision_geometry),velocity_grid(initial_velocity_grid),velocity_field(initial_velocity_field)
    {
        rigid_body_collection.rigid_body_particles.kinematic(particle_index)=true;
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset=LEVELSET_IMPLICIT_OBJECT<TV>::Create();
        levelset->levelset.grid.Initialize(TV_INT::All_Ones_Vector(),RANGE<TV>::Unit_Box());
        levelset->levelset.phi.Resize(levelset->levelset.grid.Domain_Indices());
        levelset->update_every_frame=true;
        Add_Structure(*levelset);
    }

    ~KINEMATIC_COLLISION_BODY() {}

    void Initialize_Implicit_Object_Levelset(const TV_INT counts,const RANGE<TV>& box)
    {
        LEVELSET_IMPLICIT_OBJECT<TV>* levelset=(LEVELSET_IMPLICIT_OBJECT<TV>*)implicit_object->object_space_implicit_object;
        levelset->levelset.grid.Initialize(counts,box);
        if(levelset->levelset.cell_range)
        {
            delete levelset->levelset.cell_range;
            levelset->levelset.cell_range=0;
        }
        levelset->levelset.phi.Resize(levelset->levelset.grid.Domain_Indices());
        levelset->Update_Box();levelset->Update_Minimum_Cell_Size();
        delete levelset->levelset.cell_range;
    }

    TV Pointwise_Object_Velocity(const TV& X) const
    {return Internal_Pointwise_Object_Velocity(X);}

    TV Pointwise_Object_Velocity_At_Particle(const TV& X,const int particle_index) const
    {return Internal_Pointwise_Object_Velocity(X);}

    TV Pointwise_Object_Velocity(const int simplex_id,const TV& X) const
    {return Internal_Pointwise_Object_Velocity(X);}

    TV Pointwise_Object_Pseudo_Velocity(const int simplex_id,const TV& X,const int state1,const int state2) const
    {PHYSBAM_NOT_IMPLEMENTED();}
private:
    TV Internal_Pointwise_Object_Velocity(const TV& X) const
    {assert(velocity_field); // this version of the function should only be called if velocity field is defined
        TV object_V=interpolation.Clamped_To_Array(*velocity_grid,*velocity_field,Object_Space_Point(X));
        return object_V;}
};
}
#endif
