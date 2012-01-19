//#####################################################################
// Copyright 2007, Andrew Selle, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Simulation of Hair
//#####################################################################
#include "MOCAP_FLESH.h"
using namespace PhysBAM;
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include "MOCAP_FLESH.h"

using namespace PhysBAM;
//#####################################################################
// Function MOCAP_FLESH
//#####################################################################
template<class T_input> MOCAP_FLESH<T_input>::
MOCAP_FLESH(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),steps_per_frame(1)
{
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Register_Options()
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-steps",1,"steps per frame");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Parse_Options()
{
    BASE::Parse_Options();
    if(parse_args->Is_Value_Set("-steps")) steps_per_frame=parse_args->Get_Integer_Value("-steps");
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Limit_Solids_Dt(T& dt,const T time)
{
    dt=clamp(dt,(T)0,one_over_frame_rate/(T)steps_per_frame);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Initialize_Bodies()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    ARRAY<float> masses;

    output_directory=STRING_UTILITIES::string_sprintf("Mocap_Flesh/output_%d",test_number);

    solids_evolution=new QUASISTATIC_EVOLUTION<TV>(solids_parameters,solid_body_collection);
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
    solids_parameters.implicit_solve_parameters.cg_iterations=500;
    solids_parameters.newton_iterations=20;
    solids_parameters.newton_tolerance=(T)1e-2;
    solid_body_collection.print_residuals=true;
    solids_parameters.use_partially_converged_result=true;
    one_over_frame_rate=(T)1/frame_rate;
    LOG::cout<<"MONITOR begin_frame="<<this->first_frame<<std::endl;
    LOG::cout<<"MONITOR output_directory="<<(FILE_UTILITIES::Get_Working_Directory()+"/"+output_directory)<<std::endl;

    // read data
    TETRAHEDRALIZED_VOLUME<T>& body_original=*TETRAHEDRALIZED_VOLUME<T>::Create();
    TETRAHEDRALIZED_VOLUME<T>& body=*(TETRAHEDRALIZED_VOLUME<T>*)body_original.Append_Particles_And_Create_Copy(particles);
    FILE_UTILITIES::Read_From_File<T>("body.tet",body);
    FILE_UTILITIES::Read_From_File<T>("attachments",attachments);
    FILE_UTILITIES::Read_From_File<T>("body_motion.mtn",body_motion);

    last_frame=body_motion.trajectories(1).counts.x;
    LOG::cout<<"MONITOR end_frame="<<last_frame<<std::endl;

    //int offset=particles.array_collection->Size();
    TRIANGULATED_SURFACE<T>& surface_original=*TRIANGULATED_SURFACE<T>::Create();
    FILE_UTILITIES::Read_From_File<T>("body.tri",surface_original);
    TRIANGULATED_SURFACE<T>& surface=*(TRIANGULATED_SURFACE<T>*)surface_original.Append_Particles_And_Create_Copy(particles);
    surface.Update_Number_Nodes();
    body.Update_Number_Nodes();

    ARRAY<int> tets;ARRAY<PAIR<int,TV> > bindings;const T tolerance=1;
    body.Initialize_Hierarchy();
    for(int p=0;p<surface_original.particles.array_collection->Size();p++){
        tets.Remove_All();body.hierarchy->Intersection_List(surface_original.particles.X(p),tets,1e-4);bool got_bind=false;
        for(int tt=0;tt<tets.m;tt++){int t=tets(tt);
            TV bary=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(surface_original.particles.X(p),body.particles.X.Subset(body.mesh.elements(t)));
            if(bary.x>-tolerance && bary.y>-tolerance && bary.z>-tolerance && bary.x+bary.y+bary.z<(T)1+tolerance){bindings.Append(PAIR<int,TV>(t,bary));got_bind=true;break;}}
        if(!got_bind){LOG::cout<<"no binding on particle "<<p<<std::endl;bindings.Append(PAIR<int,TV>(0,TV(0,0,0)));}}

    /*for(int i=0;i<bindings.m;i++){int p=offset+i;
        VECTOR<int,4> nodes=body.mesh.elements(bindings(i).x);
        binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,p,nodes,bindings(i).y));}*/
    particles.Store_Velocity();

    // mass
    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(body,density,true);
    //################################################################
    // Mass Update Phase
    //################################################################
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    solid_body_collection.Add_Force(Create_Quasistatic_Finite_Volume(body,new NEO_HOOKEAN<T,3>((T)5e3,(T).4,(T).01,(T).25),true,true));
    //solid_body_collection.Add_Force(Create_Finite_Volume(body,new NEO_HOOKEAN<T,3>((T)5e2,(T).4,(T).01,(T).25),true,(T).1));
    deformable_body_collection.deformable_geometry.Add_Structure(&body);
    deformable_body_collection.deformable_geometry.Add_Structure(&surface);

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        solid_body_collection.deformable_body_collection.mpi_solids->KD_Tree_Partition(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,particles.X);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    for(int joint=0;joint<attachments.m;joint++){
        const ARRAY<PAIR<int,TV> >& joint_attachment=attachments(joint);
        for(int i=0;i<joint_attachment.m;i++) V(joint_attachment(i).x)=TV();}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    T frame_float=(velocity_time*frame_rate+1);
    int frame_lower=(int)frame_float;
    for(int joint=0;joint<attachments.m;joint++){
        const ARRAY<PAIR<int,TV> >& joint_attachment=attachments(joint);
        for(int i=0;i<joint_attachment.m;i++){
            int p=joint_attachment(i).x;
            const TV& joint_relative_position=joint_attachment(i).y;
            V(p)=one_over_frame_rate*(body_motion.trajectories(joint)(frame_lower+1).targeted_transform*joint_relative_position-body_motion.trajectories(joint)(frame_lower).targeted_transform*joint_relative_position);}}   
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time)
{
    for(int joint=0;joint<attachments.m;joint++){
        const ARRAY<PAIR<int,TV> >& joint_attachment=attachments(joint);
        for(int i=0;i<joint_attachment.m;i++) X(joint_attachment(i).x)=TV();}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class T_input> void MOCAP_FLESH<T_input>::
Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    T frame_float=(time*frame_rate+1);
    int frame_lower=(int)frame_float;
    T alpha=frame_float-(T)frame_lower;
    for(int joint=0;joint<attachments.m;joint++){
        const ARRAY<PAIR<int,TV> >& joint_attachment=attachments(joint);
        for(int i=0;i<joint_attachment.m;i++){
            int p=joint_attachment(i).x;
            const TV& joint_relative_position=joint_attachment(i).y;
            X(p)=alpha*(body_motion.trajectories(joint)(frame_lower+1).targeted_transform*joint_relative_position)+(1-alpha)*(body_motion.trajectories(joint)(frame_lower).targeted_transform*joint_relative_position);}}   
}
//#####################################################################
template class MOCAP_FLESH<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MOCAP_FLESH<double>;
#endif
