//#####################################################################
// Copyright 2006-2008, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHOWER_CURTAIN_EXAMPLE
//##################################################################### 
#ifndef __SHOWER_CURTAIN_EXAMPLE__
#define __SHOWER_CURTAIN_EXAMPLE__

#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Solids_And_Fluids/SOLIDS_PARAMETERS_3D.h>
namespace PhysBAM{

template<class T,class RW>
class SHOWER_CURTAIN_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>
{
public:
    typedef VECTOR<T,3> TV;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::last_frame;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::write_last_frame;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::data_directory;

    SOLIDS_STANDARD_TESTS<TV,RW> tests;
    SEGMENT_MESH segment_mesh;
    T stiffness,overdamping_fraction;
    bool velocity_damping,use_blocks;
    ARRAY<int> velocity_damping_list;
    T time_to_close,close_start_time;
    VECTOR<T,3> goal_position;
    TRIANGULATED_SURFACE<T>* cloth_panel;
    int number_of_tori;
    int side_length;
    T fudge_adjustment;
    T count,torus_increment,start_position,height;
    ARRAY<FRAME_3D<T> > stored_frames;
    ARRAY<PAIR<int,TV> > stored_bindings;
    int pole_id;
    INTERPOLATION_CURVE<T,TV> motion_curve;
    ARRAY<int> rigid_body_bound_particles;
    bool hires;

    SHOWER_CURTAIN_EXAMPLE(const PARSE_ARGS& parse_args):
        SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<T,GRID<TV> >::NONE),tests(*this,solids_parameters),stiffness(1e4),overdamping_fraction(1),
        velocity_damping(false),use_blocks(false),time_to_close(20),close_start_time(5),number_of_tori(10),side_length(20),count(1),pole_id(0),hires(false)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T)1;
        last_frame=10000;
        frame_rate=24;
        output_directory="Shower_Curtain/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;
        solids_parameters.implicit_solve_parameters.cg_tolerance=1e-3;
        solids_parameters.implicit_solve_parameters.cg_iterations=500;
        solids_parameters.min_collision_loops=1;
        solids_parameters.max_collision_loops=1;
        solids_parameters.perform_per_collision_step_repulsions=false;
        solids_parameters.perform_per_time_step_repulsions=true;
        solids_parameters.perform_self_collision=true;solids_parameters.perform_collision_body_collisions=true;
        if(parse_args.Is_Value_Set("-dampen")) overdamping_fraction=(T)parse_args.Get_Double_Value("-dampen");
        if(parse_args.Is_Value_Set("-stiffen")) stiffness=(T)parse_args.Get_Double_Value("-stiffen");
        if(parse_args.Is_Value_Set("-o")) output_directory=parse_args.Get_String_Value("-o");
        if(parse_args.Is_Value_Set("-velocityDamping")) velocity_damping=true;
        if(parse_args.Is_Value_Set("-useBlocks")) use_blocks=true;
        if(parse_args.Is_Value_Set("-fudge")) fudge_adjustment=(T)parse_args.Get_Double_Value("-fudge");
        if(parse_args.Is_Value_Set("-restart")) {restart_frame=parse_args.Get_Integer_Value("-restart");restart=true;std::cout<<"restart frame found: "<<restart_frame<<std::endl;}
        std::cout<<"restart: "<<restart<<" on frame "<<restart_frame<<std::endl;
    }

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
int Add_Rigid_Body(const std::string& rigid_body_name,const std::string& filename,const FRAME_3D<T>& frame,const T scale,const bool with_phi)
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    DEFORMABLE_PARTICLES<T,TV>& particles=deformable_object.particles;
    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(filename,scale,0,with_phi);
    rigid_body.frame=frame;
    rigid_body.Set_Coefficient_Of_Restitution(0);
    rigid_body.Set_Name(rigid_body_name);
    int rigid_body_particle=tests.Add_Rigid_Body_Particle(rigid_body);
    int gravity_particle=particles.Add_Element();
    particles.mass(gravity_particle)=rigid_body.mass.mass;
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<T,TV>(particles,gravity_particle,deformable_object.rigid_body_particles,rigid_body_particle,TV()));
    stored_bindings.Append(PAIR<int,TV>(rigid_body_particle,TV()));
    rigid_body_bound_particles.Append(gravity_particle);
    return rigid_body.id_number;
}
//#####################################################################
// Function Add_Plate_Binding
//#####################################################################
void Add_Torus_Binding(const int torus_id,const int particle_id)
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    DEFORMABLE_PARTICLES<T,TV>& particles=deformable_object.particles;

    // add deformable particles for each of the sphere and plank
    VECTOR<T,3> plate_object_space_position=solids_parameters.rigid_body_parameters.list(torus_id)->frame.Inverse()*particles.X(particle_id);
    
    // add bindings to their rigid body particles
    solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<T,TV>(particles,particle_id,deformable_object.rigid_body_particles,torus_id,plate_object_space_position));
    stored_bindings.Append(PAIR<int,TV>(torus_id,plate_object_space_position));
    rigid_body_bound_particles.Append(particle_id);
}
//#####################################################################
// Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT<T,TV>& deformable_object=solid_body_collection.deformable_object;
    DEFORMABLE_PARTICLES<T,TV>& particles=deformable_object.particles;

    T x_shift=0,y_shift=0,z_shift=0;

    RIGID_BODY_STATE<TV> initial_state;initial_state.frame=FRAME_3D<T>(VECTOR<T,3>(x_shift,y_shift,z_shift),QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0)));
    if(!hires) cloth_panel=&tests.Create_Cloth_Panel(4*side_length,side_length,1,&initial_state);
    else cloth_panel=&tests.Create_Cloth_Panel(8*side_length,side_length,1,&initial_state);

    torus_increment=side_length/number_of_tori;
    start_position=-(side_length-torus_increment)*.5;
    height=(side_length-torus_increment)*.5;
    for(int torus=0;torus<number_of_tori;torus++){
        if(!hires) Add_Rigid_Body("ring","torus",FRAME_3D<T>(VECTOR<T,3>(start_position+torus*torus_increment,height+.37+fudge_adjustment,0),QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0))),.51,true);
        else Add_Rigid_Body("ring","torus",FRAME_3D<T>(VECTOR<T,3>(start_position+torus*torus_increment,height+.37+fudge_adjustment,0),QUATERNION<T>(pi/2,VECTOR<T,3>(1,0,0))),.51,true);}

    // add the rigid body bindings
    //int particle_indices[14]={328,329,330, 408,412,489, 493,570,571, 573,574,652, 653,654}; old
       int particle_indices[22]={5,84,85,86,87, 88,165,166,167,168, 169,246,247,248,249, 250,327,328,329,330, 331,410};
    //int particle_indices[54]={331,489,490,491,492,      493,494,495,649,650,      651,652,653,654,655,      656,657,810,811,817,      818,971,972,978,979,
    //                          1131,1132,1140,1141,1293, 1294,1300,1301,1454,1455, 1461,1462,1615,1616,1617, 1618,1620,1621,1622,1623, 1777,1778,1779,1780,1781,
    //                          1782,1783,1941,1942}; old
    //int particle_indices[57]={9,167,168,169,170, 171,172,173,327,328, 329,330,331,332,333, 334,335,488,489,495, 496,649,650,656,657, 
    //                          809,810,818,819,970, 971,972,978,979,980, 1132,1133,1139,1140,1293, 1294,1295,1296,1298,1299, 1300,1301,1455,1456,1457,
    //                          1458,1459,1460,1461,1618, 1619,1620};
    if(hires){
        for(int t=0;t<number_of_tori;t++) for(int p=0;p<57;p++){
            int torus_particle_number=deformable_object.rigid_body_particles.id.array.Find(t);Add_Torus_Binding(torus_particle_number,particle_indices[p]+(t-1)*16);}}
    else{
        for(int t=0;t<number_of_tori;t++) for(int p=0;p<22;p++){
            int torus_particle_number=deformable_object.rigid_body_particles.id.array.Find(t);Add_Torus_Binding(torus_particle_number,particle_indices[p]+(t-1)*8);}}
    
    // remove the unwanted triangles
    int triangle_indices[24]={324,325,326,327,328, 483,484,485,486,487, 488,489,643,644,645, 646,647,648,649,804, 805,806,807,808};    
    //int triangle_indices[24]={330,331,332,333,334, 489,490,491,492,493, 494,495,649,650,651, 652,653,654,655,810, 811,812,813,814}; old
    //int triangle_indices[79]={1294,1295,1296,1609,1610, 1611,1612,1613,1614,1615, 1616,1617,1618,1619,1929, 1930,1931,1932,1933,1934, 1935,1936,1937,1938,1939,
    //                          1940,2249,2250,2251,2252, 2253,2254,2255,2256,2257, 2258,2259,2260,2261,2568, 2569,2570,2571,2572,2573, 2574,2575,2576,2577,2578,
    //                          2579,2580,2581,2889,2890, 2891,2892,2893,2894,2895, 2896,2897,2898,2899,2900, 3209,3210,3211,3212,3213, 3214,3215,3216,3217,3218,
    //                          3219,3534,3535,3536}; old
    //int triangle_indices[79]={1290,1291,1292,1605,1606,1607,1608,1609, 1610,1611,1612,1613,1614, 1615,1925,1926,1927,1928, 1929,1930,1931,1932,1933, 1934,1935,1936,2245,2246,
    //                          2247,2248,2249,2250,2251, 2252,2253,2254,2255,2256,2257,2564,2565,2566,2567, 2568,2569,2570,2571,2572, 2573,2574,2575,2576,2577, 2885,2886,2887,2888,2889,
    //                          2890,2891,2892,2893,2894, 2895,2896,3205,3206,3207, 3208,3209,3210,3211,3212, 3213,3214,3215,3530,3531, 3532};
    ARRAY<int> triangles_to_remove;
    if(hires){for(int t=0;t<number_of_tori;t++) for(int e=0;e<79;e++) triangles_to_remove.Append(triangle_indices[e]+(t-1)*5120);}
    else{for(int t=0;t<number_of_tori;t++) for(int e=0;e<24;e++) triangles_to_remove.Append(triangle_indices[e]+(t-1)*1280);}
    Sort(triangles_to_remove);for(int i=triangles_to_remove.m;i>=1;i--) cloth_panel->mesh.elements.Remove_Index(triangles_to_remove(i));
 
    deformable_object.collisions.collision_structures.Append(cloth_panel);

   // add the rod
    if(restart){
        if(!hires) pole_id=Add_Rigid_Body("pole","Rings_Test/tall_cylinder",FRAME_3D<T>(VECTOR<T,3>(0,height+.37+fudge_adjustment,0),QUATERNION<T>(pi/2,VECTOR<T,3>(0,0,1))),.5,true);
        else pole_id=Add_Rigid_Body("pole","Rings_Test/tall_cylinder",FRAME_3D<T>(VECTOR<T,3>(0,height+.37+fudge_adjustment,0),QUATERNION<T>(pi/2,VECTOR<T,3>(0,0,1))),.5,true);
        solids_parameters.rigid_body_parameters.list(pole_id)->Update_Bounding_Box();
        solids_parameters.rigid_body_parameters.list(pole_id)->is_kinematic=true;}
    else{
        if(!hires) Add_Rigid_Body("pole","Rings_Test/tall_cylinder",FRAME_3D<T>(VECTOR<T,3>(-30,height+.37+fudge_adjustment,0),QUATERNION<T>(pi/2,VECTOR<T,3>(0,0,1))),.5,true);
        else Add_Rigid_Body("pole","Rings_Test/tall_cylinder",FRAME_3D<T>(VECTOR<T,3>(-30,height+.37+fudge_adjustment,0),QUATERNION<T>(pi/2,VECTOR<T,3>(0,0,1))),.5,true);}

    //tests.Add_Ground();
    //solids_parameters.rigid_body_parameters.list(number_of_tori+2)->frame.t.y=-13;
    // correct number nodes
    for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();

    // make self_collision surface
    TRIANGULATED_SURFACE<T>* self_collision_surface=TRIANGULATED_SURFACE<T>::Create(particles);
    self_collision_surface->mesh.Initialize_Mesh(cloth_panel->mesh);self_collision_surface->Update_Number_Nodes();
    for(int e=self_collision_surface->mesh.elements.m;e>=1;e--){
        // remove elements with any particle being in the rigid_body_bound_particles list
        int i,j,k;self_collision_surface->mesh.elements.Get(e,i,j,k);
        int i_found=rigid_body_bound_particles.Find(i);
        int j_found=rigid_body_bound_particles.Find(j);
        int k_found=rigid_body_bound_particles.Find(k);
        if(i_found || j_found || k_found) self_collision_surface->mesh.elements.Remove_Index(e);}

    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(self_collision_surface);

    // correct mass
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    // add forces
    if(restart) solid_body_collection.Add_Force(new GRAVITY<T,TV>(particles));
    solid_body_collection.Add_Force(Create_Edge_Springs<T>(*cloth_panel,stiffness*1/(1+sqrt((T)2)),overdamping_fraction*5)); // were *2 and *10
    //solid_body_collection.Add_Force(Create_Bending_Elements<T>(*cloth_panel,1e-4));
    solid_body_collection.Add_Force(Create_Bending_Springs<T>(*cloth_panel,2/(1+sqrt((T)2)),8));

    //solid_body_collection.Add_Force(Create_Altitude_Springs<T>(*cloth_panel,stiffness*20/(1+sqrt((T)2)),overdamping_fraction*15)); // were *2 and *10
    //solid_body_collection.Add_Force(Create_Bending_Springs<T>(*cloth_panel,2/(1+sqrt((T)2)),8));

//    solid_body_collection.Add_Force(Create_Edge_Springs<T>(*cloth_panel,1e2,2));
//    solid_body_collection.Add_Force(Create_Altitude_Springs<T>(*cloth_panel,1e2,2));
//    solid_body_collection.Add_Force(Create_Bending_Elements<T>(*cloth_panel,3e-3));


    solids_parameters.collision_body_list.Add_Body(solids_parameters.rigid_body_parameters.list(number_of_tori+1));
    solids_parameters.collision_body_list.Add_Body(solids_parameters.rigid_body_parameters.list(number_of_tori+2));

    solid_body_collection.Update_Fragments();
    LOG::cout<<"NUMBER of FRAGMENTS="<<solid_body_collection.deformable_object.particles_of_fragment.m<<std::endl;
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Initialize_Bodies();

    int time_adjustment=56;

    motion_curve.Add_Control_Point(0,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(1,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(3,VECTOR<T,3>(0,height+.37+fudge_adjustment,-2));
    motion_curve.Add_Control_Point(7,VECTOR<T,3>(0,height+.37+fudge_adjustment,2));
    motion_curve.Add_Control_Point(9,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(11,VECTOR<T,3>(0,height+.37+fudge_adjustment+2,0));
    motion_curve.Add_Control_Point(15,VECTOR<T,3>(0,height+.37+fudge_adjustment-2,0));
    motion_curve.Add_Control_Point(17,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(30,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));

    motion_curve.Add_Control_Point(0+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(1+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(3+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment,-2));
    motion_curve.Add_Control_Point(7+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment,2));
    motion_curve.Add_Control_Point(9+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(11+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment+2,0));
    motion_curve.Add_Control_Point(15+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment-2,0));
    motion_curve.Add_Control_Point(17+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));
    motion_curve.Add_Control_Point(30+time_adjustment,VECTOR<T,3>(0,height+.37+fudge_adjustment,0));

}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
//    if(current_position_time < 56 && current_position_time > 55) for(int i=0;i<V.m;i++) V(i)=TV();
//    if(current_position_time > 17+55 && current_position_time < 11) V(6561)=TV(.5,.1,.2); // pulled 2 uses .5 for Y, pulled 3 uses .1 for Y
    if(current_position_time<1) for(int i=0;i<V.m;i++) V(i)=TV();
    if(current_position_time>17 && current_position_time<56) V(6561)=TV(.5,.1,.2); // pulled 2 uses .5 for Y, pulled 3 uses .1 for Y
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    //if(current_position_time < 56 && current_position_time > 55) for(int i=0;i<V.m;i++) V(i)=TV();
    //if(current_position_time > 55+17) V(6561)=TV();
    if(current_position_time < 1) for(int i=0;i<V.m;i++) V(i)=TV();
    if(current_position_time > 17) V(6561)=TV();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{ 
    if(restart)    twist(number_of_tori+1)=TWIST<TV>();
    else{
        for(int t=0;t<number_of_tori;t++) twist(t).angular=TV(0,pow(-1.,t)*.1,0);
        for(int t=0;t<number_of_tori;t++) twist(t).linear=t<=number_of_tori/2?TV(.4,0,0):TV(-.4,0,0);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(restart) twist(number_of_tori+1)=TWIST<TV>();
    else for(int t=0;t<number_of_tori;t++) twist(t)=TWIST<TV>();
}
//#####################################################################
// Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    LOG::cout<<" * * * * * * * * * * * * * * preprocess frame called for frame "<<frame<<std::endl;
    // rearrange tori
    if(frame==12){
        std::cout<<"putting bar in correct position"<<std::endl;
        solids_parameters.rigid_body_parameters.list(number_of_tori+1)->frame.t.x=0;
        std::cout<<"adding gravity to particles"<<std::endl;}
    T change=(number_of_tori+1)/2;
    if(frame<12){
        for(int t=0;t<number_of_tori;t++){
            solids_parameters.rigid_body_parameters.list(t)->frame.t+=TV((change-t)*.05,0,0);
            solids_parameters.rigid_body_parameters.list(t)->frame.r=QUATERNION<T>(pow(-1.,t)*.1,VECTOR<T,3>(0,1,0))*solids_parameters.rigid_body_parameters.list(t)->frame.r;}}
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    frame=FRAME_3D<T>(motion_curve.Value(time),QUATERNION<T>(pi/2,VECTOR<T,3>(0,0,1)));
}
//#####################################################################
// Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    int count=1;
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T,GRID<TV>,RW>::Read_Output_Files_Solids(frame);

    for(int i=0;i<solid_body_collection.deformable_body_collection.binding_list.bindings.m;i++)
        if(RIGID_BODY_BINDING<T,TV>* binding=dynamic_cast<RIGID_BODY_BINDING<T,TV>*>(solid_body_collection.deformable_body_collection.binding_list.bindings(i))){
            binding->rigid_body_particles=&solid_body_collection.deformable_object.rigid_body_particles;
            binding->particle_index=stored_bindings(count).x;
            binding->object_space_position=stored_bindings(count).y;
            count++;}
}
};
}
#endif
