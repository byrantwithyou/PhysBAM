//#####################################################################
// Copyright 2005-2008, Elliot English, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_COLLISION_GEOMETRY
//#####################################################################
#ifndef __RIGID_COLLISION_GEOMETRY__
#define __RIGID_COLLISION_GEOMETRY__

#include <Tools/Utilities/NONCOPYABLE.h>
#include <Tools/Vectors/VECTOR_1D.h>
#include <Tools/Vectors/VECTOR_2D.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;

template<class TV>
class RIGID_COLLISION_GEOMETRY_BASE:public COLLISION_GEOMETRY<TV>
{
private:
    typedef typename TV::SCALAR T;
    typedef COLLISION_GEOMETRY<TV> BASE;

public:
    using BASE::collision_thickness;

    typedef TV VECTOR_T;
    RIGID_BODY<TV>& rigid_body;
    ARRAY<PAIR<FRAME<TV>,T> > saved_states;

    RIGID_COLLISION_GEOMETRY_BASE(RIGID_BODY<TV>& rigid_body_input);
    virtual ~RIGID_COLLISION_GEOMETRY_BASE();

//#####################################################################
    TV Pointwise_Object_Velocity(const TV& X) const override;
    TV Pointwise_Object_Velocity(const int simplex_id,const TV& X) const override; // extra simplex_id is not used, but for a virtual function in COLLISION GEOMETRY
    virtual TV Pointwise_Object_Velocity_At_Particle(const TV& X,const int particle_index) const; // overridden in DEFORMABLE_OBJECT_WRAPPER where X is ignored
    void Update_Bounding_Box() override;
    const RANGE<TV>& Axis_Aligned_Bounding_Box() const override;
    bool Has_Volumetric_Geometry() const override;

    TV Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2=0,const int max_iterations=1,int* simplex_id=0,T* distance=0) const override;

    T Implicit_Geometry_Extended_Value(const TV& location) const override;
    TV Implicit_Geometry_Normal(const TV& location,const int aggregate=-1) const override;
    TV Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const override;
    TV Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const override;
    TV Implicit_Geometry_Closest_Point_On_Boundary(const TV& location,const T tolerance=0,const int max_iterations=1,T* distance=0) const override;
    bool Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value=0) const override;
    bool Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value=0) const override;
    bool Implicit_Geometry_Lazy_Inside_Extended_Levelset(const TV& location,T contour_value=0) const override;

    TV Simplex_World_Space_Point_From_Barycentric_Coordinates(const int simplex_id,const TV& weights) const override;
    bool Simplex_Intersection(RAY<TV>& ray) const override;
    bool Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const override;
    bool Inside_Any_Simplex(const TV& location,int& simplex_id) const override;
    bool Inside(const TV& location,const T thickness_over_two) const override;
    TV Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2=0,int* simplex_id=0,T* returned_distance=0) const override;
    TV Pointwise_Object_Pseudo_Velocity(const int simplex_id,const TV& X,const int state1,const int state2) const override;
    int Number_Of_Simplices() const override;

    void Save_State(const int state_index,const T time=0) override;
    void Restore_State(const int state_index) override;
    void Average_States(const int state1, const int state2, const int result_state,const T interpolation_distance) override;
    void Delete_State(const int state_index) override;
    void Read_State(TYPED_ISTREAM& input,const int state_index) override;
    void Write_State(TYPED_OSTREAM& output,const int state_index) const override;
//#####################################################################
};
}
#endif
