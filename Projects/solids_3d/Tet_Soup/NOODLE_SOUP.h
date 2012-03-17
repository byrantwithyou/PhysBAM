//#####################################################################
// Copyright 2007-2008, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NOODLE_SOUP
//##################################################################### 
#ifndef __NOODLE_SOUP__
#define __NOODLE_SOUP__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
#include <Solids_And_Fluids/SOLIDS_PARAMETERS_3D.h>
namespace PhysBAM{

template<class T,class RW>
class NOODLE_SOUP:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;
    using BASE::write_last_frame;using BASE::data_directory;using BASE::solids_parameters;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;
    T ribbon_length,ribbon_width;
    T bowl_scale;
    std::string noodle_data_directory;
    DEFORMABLE_OBJECT<T,TV> noodle_deformable_object;
    PARTICLE_SUBSET<DEFORMABLE_PARTICLES<T,TV> > gravity_particles;
    ARRAY<TRIANGULATED_SURFACE<T>*> noodles;

    NOODLE_SOUP(const PARSE_ARGS& parse_args):
        SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<T,GRID<TV> >::NONE),tests(*this,solids_parameters),ribbon_length(250),ribbon_width(2.5),bowl_scale(10),
       noodle_data_directory("/data/shinar/PhysBAM/Projects/solids_3d/Tet_Soup/noodle/output/"),gravity_particles(solid_body_collection.deformable_object.particles)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T)5;
        last_frame=10000;
        frame_rate=24;
        output_directory="Tet_Soup/noodle_soup/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-6;
        solids_parameters.implicit_solve_parameters.cg_iterations=500;
        solids_parameters.perform_self_collision=true;
        solids_parameters.collide_with_interior=false;
    }

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    DEFORMABLE_PARTICLES<T,TV>& particles=deformable_object.particles;

    // read coarse sphere to get the dynamic particles
    TETRAHEDRALIZED_VOLUME<T>& dynamic_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
        RIGID_BODY_STATE<TV>(FRAME_3D<T>(TV(0,0,0))),true,true);
    T radius=0;for(int p=0;p<particles.Size();p++) radius=max(radius,particles.X(p).Magnitude());
    ARRAY<T>::copy(0,particles.mass.array);

    // add the noodles and compute their bounding box
    int noodle_particles_start=particles.Size()+1;
    Add_Noodles();
    BOX_3D<T> box;box.Enlarge_To_Include_Points(particles.X.array);
    T scale=1.1*maxabs(box.xmin,box.xmax,box.ymin,box.ymax,box.xmin,box.xmax);

    // enlarge the embedding sphere
    for(int p=1;p<noodle_particles_start;p++) particles.X(p)*=scale;
    dynamic_volume.Initialize_Tetrahedron_Hierarchy();

    // add bindings
    for(int p=noodle_particles_start;p<=particles.Size();p++){
        gravity_particles.active_indices.Append(p);
        int parent_tet=Get_Intersecting_Tetrahedron(particles,particles.X(p),dynamic_volume);if(parent_tet<0) PHYSBAM_FATAL_ERROR();
        VECTOR<int,4> parents=dynamic_volume.mesh.elements(parent_tet);
        VECTOR<T,3> weights=TETRAHEDRON<T>(particles.X.Subset(parents)).Barycentric_Coordinates(particles.X(p));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new LINEAR_BINDING<T,TV,4>(particles,p,parents,weights));}
    gravity_particles.Update_Subset_Index_From_Element_Index();

    particles.X.array+=TV(1.5*scale,(T)15*scale,0);

    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();
    for(int i=0;i<noodles.m;i++) noodles(i)->Update_Number_Nodes();

    for(int n=0;n<noodles.m;n++){
        TRIANGULATED_SURFACE<T>& colliding_surface=tests.Create_Drifted_Surface(*noodles(n),solid_body_collection.deformable_body_collection.soft_bindings,true);
        colliding_surface.Update_Number_Nodes();
        deformable_object.Add_Structure(&colliding_surface);
        deformable_object.collisions.collision_structures.Append(&colliding_surface);}

    // correct number nodes
    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    tests.Add_Ground(.3,20);

    // add planes for the noodle soup to collide with
    RIGID_BODY<TV>& plane1=tests.Add_Rigid_Body("ground",.2,.8);plane1.frame.t+=TV(2*scale,12*scale,0);plane1.frame.r=QUATERNION<T>::From_Euler_Angles(0,0,pi/4);
    RIGID_BODY<TV>& plane2=tests.Add_Rigid_Body("ground",.2,.8);plane2.frame.t+=TV(-2*scale,8*scale,0);plane2.frame.r=QUATERNION<T>::From_Euler_Angles(0,0,-pi/4);
    RIGID_BODY<TV>& plane3=tests.Add_Rigid_Body("ground",.2,.8);plane3.frame.t+=TV(2*scale,4*scale,0);plane3.frame.r=QUATERNION<T>::From_Euler_Angles(0,0,pi/4);

    // collisions
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents(particles.mass.array);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    // add forces
    for(int n=0;n<noodles.m;n++){
        TRIANGULATED_SURFACE<T>& noodle=*noodles(n);
        solid_body_collection.Add_Force(Create_Edge_Springs<T>(noodle,(T)1e4,1));
        TRIANGLE_BENDING_ELEMENTS<T>* bend=Create_Bending_Elements(noodle,(T)1e3);
        bend->Set_Area_Cutoff_From_Triangulated_Surface(noodle,.1);bend->use_force_differential=false;
        solid_body_collection.Add_Force(bend);}
    solid_body_collection.Add_Force(new GRAVITY<T,TV>(gravity_particles));

    solid_body_collection.Update_Fragments();

    LOG::cout<<"NUMBER of FRAGMENTS="<<solid_body_collection.deformable_object.particles_of_fragment.m<<std::endl;
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Get_Intersecting_Tetrahedron
//#####################################################################
int Get_Intersecting_Tetrahedron(const DEFORMABLE_PARTICLES<T,TV>& particles,const TV& location,const TETRAHEDRALIZED_VOLUME<T>& dynamic_volume)
{
    ARRAY<int> intersection_list;dynamic_volume.tetrahedron_hierarchy->Intersection_List(location,intersection_list);
    for(int i=0;i<intersection_list.m;i++) if(TETRAHEDRON<T>(particles.X.Subset(dynamic_volume.mesh.elements(intersection_list(i)))).Inside(location)) return intersection_list(i);
    return 0;
}
//#####################################################################
// Function Read_Noodle_Triangulated_Surface
//#####################################################################
void Read_Noodle_Triangulated_Surface(const int frame_number,const FRAME_3D<T>& frame,const bool read_static_variables=false)
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    noodle_deformable_object.Read(stream_type,noodle_data_directory,noodle_data_directory,frame_number,-1,read_static_variables);
    for(int p=0;p<noodle_deformable_object.particles.Size();p++) noodle_deformable_object.particles.X(p)=frame*noodle_deformable_object.particles.X(p);
    noodles.Append((TRIANGULATED_SURFACE<T>*)noodle_deformable_object.structures(1)->Append_Particles_And_Create_Copy(deformable_object.particles));
}
//#####################################################################
// Function Add_Noodles
//#####################################################################
void Add_Noodles()
{
    Read_Noodle_Triangulated_Surface(412,FRAME_3D<T>(TV(0,-10,0)),true);
    Read_Noodle_Triangulated_Surface(398,FRAME_3D<T>(TV(0,10,0),QUATERNION<T>::From_Euler_Angles(pi,.15*pi,0)));
    Read_Noodle_Triangulated_Surface(375,FRAME_3D<T>(TV(-10,0,0),QUATERNION<T>::From_Euler_Angles(.5*pi,.5*pi,0)));
    Read_Noodle_Triangulated_Surface(385,FRAME_3D<T>(TV(10,0,0),QUATERNION<T>::From_Euler_Angles(.5*pi,-.5*pi,0)));
    Read_Noodle_Triangulated_Surface(403,FRAME_3D<T>(TV(0,0,-10),QUATERNION<T>::From_Euler_Angles(.5*pi,0,0)));
    Read_Noodle_Triangulated_Surface(392,FRAME_3D<T>(TV(0,0,10),QUATERNION<T>::From_Euler_Angles(-.5*pi,0,0)));
}
//#####################################################################
};
}
#endif
