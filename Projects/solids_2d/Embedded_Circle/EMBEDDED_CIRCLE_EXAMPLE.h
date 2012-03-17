//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDED_CIRCLE_EXAMPLE
//#####################################################################
#ifndef __EMBEDDED_CIRCLE_EXAMPLE__
#define __EMBEDDED_CIRCLE_EXAMPLE__

#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_SPLINE_MODEL_2D.h>
#include <Forces_And_Torques/BODY_FORCES_2D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_2D.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
namespace PhysBAM{

template<class T,class RW>
class EMBEDDED_CIRCLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;using BASE::fluids_parameters;

    CIRCLE<T> circle;
    T initial_height;
    VECTOR_2D<T> initial_velocity;
    T initial_angular_velocity;
    int m_input,n_input;

    EMBEDDED_CIRCLE_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),initial_height(1),initial_velocity(0,0),initial_angular_velocity(0)
    {
        m_input=13;n_input=13;
        last_frame=240;
        restart=false;restart_frame=0;  
        solids_parameters.cfl=(T)1;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Embedded_Circle/output";

        solids_parameters.use_constant_mass=true;
    }

    virtual ~EMBEDDED_CIRCLE_EXAMPLE()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
virtual void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Area();
    EMBEDDED_TRIANGULATED_AREA<T>& embedded_triangulated_area=*solids_parameters.deformable_body_parameters.list(index).embedded_triangulated_area;
    TRIANGLES_OF_MATERIAL_2D<T>& triangles_of_material=*solids_parameters.deformable_body_parameters.list(index).triangles_of_material;
    TRIANGULATED_AREA<T>& triangulated_area=embedded_triangulated_area.triangulated_area;
    DEFORMABLE_PARTICLES<T,VECTOR_2D<T> >& particles=triangulated_area.particles;
    
    RED_GREEN_GRID_2D<T> red_green_grid;
    GRID<TV> red_green_uniform_grid(4*m_input+1,4*n_input+1,circle.Bounding_Box());
    red_green_grid.Initialize(red_green_uniform_grid,1);
    red_green_grid.Build_Triangulated_Area(triangulated_area);
    ARRAY<T> phi(particles.Size());for(int p=0;p<phi.m;p++) phi(p)=circle.Signed_Distance(particles.X(p));
    std::cout<<particles.X.array.m<<" "<<particles.V.array.m<<", s "<<particles.store_velocity<<"\n";
    embedded_triangulated_area.Calculate_Segmented_Curve_From_Levelset_On_Triangle_Nodes(phi);
    triangles_of_material.Create_Material_Area();
    // TODO:  ADD Update Mass to Boundary Surface.
   
    triangulated_area.Set_Density(1000);triangulated_area.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);   
    triangulated_area.Update_Bounding_Box();
    VECTOR_2D<T> center(triangulated_area.bounding_box->Center());T bottom=triangulated_area.bounding_box->ymin;
    for(int i=0;i<particles.Size();i++){
        particles.V(i)=initial_velocity+initial_angular_velocity*(particles.X(i)-center).Rotate_Counterclockwise_90();
        particles.X(i).y+=initial_height-bottom;}

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>("../../Public_Data/Rigid_Bodies_2D/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT_2D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TRIANGULATED_AREA<T>& triangulated_area=deformable_object.embedded_triangulated_area->triangulated_area;

    deformable_object.Delete_Forces();
    solid_body_collection.Add_Force(Create_Body_Forces<T>(triangulated_area));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(triangulated_area,new DIAGONALIZED_SPLINE_MODEL_2D<T>((T)5e4,(T).3,(T).5,7,(T).01)));
}
//#####################################################################
};
}
#endif
