//#####################################################################
// Copyright 2007-2008, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TET_SOUP_EXAMPLE
//##################################################################### 
#ifndef __TET_SOUP_EXAMPLE__
#define __TET_SOUP_EXAMPLE__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
#include <Solids_And_Fluids/SOLIDS_PARAMETERS_3D.h>
namespace PhysBAM{

template<class T,class RW>
class TET_SOUP_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;
    using BASE::write_last_frame;using BASE::data_directory;using BASE::solids_parameters;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;
    int number_embedded_tets;
    RANDOM_NUMBERS random_numbers;
    BOX_3D<T> unit_box;

    TET_SOUP_EXAMPLE(const PARSE_ARGS& parse_args):
        SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<T,GRID<TV> >::NONE),tests(*this,solids_parameters),number_embedded_tets(2000),unit_box(0,1,0,1,0,1)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).2;
        last_frame=10000;
        frame_rate=24;
        output_directory="Tet_Soup/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-6;
        solids_parameters.implicit_solve_parameters.cg_iterations=500;
        solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=true;
        random_numbers.Set_Seed(1234);
    }

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Get_Intersecting_Tetrahedron
//#####################################################################
int Get_Intersecting_Tetrahedron(const DEFORMABLE_PARTICLES<T,TV>& particles,const TV& location,const TETRAHEDRALIZED_VOLUME<T>& dynamic_volume)
{
    ARRAY<int> intersection_list;dynamic_volume.tetrahedron_hierarchy->Intersection_List(location,intersection_list);
    for(int i=0;i<intersection_list.m;i++) if(TETRAHEDRON<T>(particles.X.Subset(dynamic_volume.mesh.elements(intersection_list(i)))).Inside(location)) return intersection_list(i);
    return 0;
}
//#####################################################################
// Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    DEFORMABLE_PARTICLES<T,TV>& particles=deformable_object.particles;

    // read coarse sphere to get the dynamic particles
    TETRAHEDRALIZED_VOLUME<T>& dynamic_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
        RIGID_BODY_STATE<TV>(FRAME_3D<T>(TV(0,(T)3,0))),true,true);
    dynamic_volume.Initialize_Tetrahedron_Hierarchy();
    ARRAY<T>::copy(0,particles.mass.array);

    TETRAHEDRALIZED_VOLUME<T>& tet_soup=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);

    // create random embedded tet soup
    TETRAHEDRON<T> tet(particles.X.Subset(dynamic_volume.mesh.elements(1)));tet.X=tet.X-tet.Center();
    for(int i=0;i<number_embedded_tets;i++){
        int element=random_numbers.Get_Uniform_Integer(1,dynamic_volume.mesh.elements.m);VECTOR<T,3> weights=random_numbers.Get_Direction<TV>();
        VECTOR<T,3> translation=TETRAHEDRON<T>(particles.X.Subset(dynamic_volume.mesh.elements(element))).Point_From_Barycentric_Coordinates(weights);
        T scale=random_numbers.Get_Uniform_Number((T).2,(T)2);
        QUATERNION<T> rotation=QUATERNION<T>::From_Rotation_Vector(random_numbers.Get_Direction<TV>());
        ARRAY<int> intersection_list;ARRAY<TV> X(4);
        for(int p=0;p<4;p++){
            X(p)=translation+(rotation.Rotate(tet.X(p)))*scale;
            int intersecting_tet=Get_Intersecting_Tetrahedron(particles,X(p),dynamic_volume);if(intersecting_tet<0) break;
            intersection_list.Append(intersecting_tet);}
        if(intersection_list.m!=4){i--;continue;}
        int offset=particles.Size();particles.Add_Elements(4);
        tet_soup.mesh.elements.Append(VECTOR<int,4>(1,2,3,4)+offset);
        for(int p=0;p<4;p++){
            VECTOR<int,4> parents=dynamic_volume.mesh.elements(intersection_list(p));
            VECTOR<T,3> weights=TETRAHEDRON<T>(particles.X.Subset(parents)).Barycentric_Coordinates(X(p));
            particles.mass(offset+p)=1;
            LINEAR_BINDING<T,TV,4>* binding=new LINEAR_BINDING<T,TV,4>(particles,offset+p,parents,weights);
            solid_body_collection.deformable_body_collection.binding_list.Add_Binding(binding);
            particles.X(offset+p)=binding->Embedded_Position();}}

    deformable_object.Add_Structure(&tet_soup);
    tests.Add_Ground();

    // correct number nodes
    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    tet_soup.Initialize_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& colliding_surface=tests.Create_Drifted_Surface(*tet_soup.triangulated_surface,solid_body_collection.deformable_body_collection.soft_bindings,true); //TODO: use binding springs?
    colliding_surface.Update_Number_Nodes();
    deformable_object.Add_Structure(&colliding_surface);

    // collisions
//    deformable_object.collisions.collision_structures.Append(&dynamic_volume);
    deformable_object.collisions.collision_structures.Append(&colliding_surface);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_object.structures);
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    // add forces
    solid_body_collection.Add_Force(new GRAVITY<T,TV>(tet_soup));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tet_soup,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));

    solid_body_collection.Update_Fragments();

    LOG::cout<<"NUMBER of FRAGMENTS="<<solid_body_collection.deformable_object.particles_of_fragment.m<<std::endl;
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
};
}
#endif
