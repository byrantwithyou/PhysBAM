//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Unnur Gretarsdottir, Jon Gretarsson, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Math_Tools/cyclic_shift.h>
#include <Core/Math_Tools/Robust_Arithmetic.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SEGMENT_BENDING_SPRINGS<TV>::
SEGMENT_BENDING_SPRINGS(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh)
    :LINEAR_SPRINGS<TV>(particles,bending_segment_mesh)
{
    Initialize(segment_mesh);
    Set_Stiffness(0);Set_Damping(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SEGMENT_BENDING_SPRINGS<TV>::
~SEGMENT_BENDING_SPRINGS()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SEGMENT_BENDING_SPRINGS<TV>::
Initialize(SEGMENT_MESH& segment_mesh_input)
{
   if(!segment_mesh_input.adjacent_elements) segment_mesh_input.Initialize_Adjacent_Elements();
    ARRAY<VECTOR<int,2> > segment_list;
    for(int s1=0;s1<segment_mesh_input.elements.m;s1++){
        VECTOR<int,2> s1_nodes=segment_mesh_input.elements(s1);
        const ARRAY<int>& adjacent=(*segment_mesh_input.adjacent_elements)(s1);
        for(int a=0;a<adjacent.m;a++){int s2=adjacent(a);
            if(s2>s1){
                VECTOR<int,2> s2_nodes=segment_mesh_input.elements(s2);
                for(int i=0;i<2;i++){
                    int j=s2_nodes.Find(s1_nodes[i]);
                    if(j>=0) segment_list.Append(s1_nodes.Remove_Index(i).Append(s2_nodes.Remove_Index(j).x));}}}}
    segment_mesh.Initialize_Mesh(segment_mesh_input.number_nodes,segment_list);
}
//#####################################################################
// Function Create_Segment_Bending_Springs
//#####################################################################
template<class TV> SEGMENT_BENDING_SPRINGS<TV>* PhysBAM::
Create_Segment_Bending_Springs(DEFORMABLE_PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const typename TV::SCALAR stiffness,const typename TV::SCALAR overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const typename TV::SCALAR max_strain_per_time_step,const bool use_rest_state_for_strain_rate,
    const typename TV::SCALAR restlength_enlargement_fraction,const bool verbose)
{
    SEGMENT_BENDING_SPRINGS<TV>* bend=new SEGMENT_BENDING_SPRINGS<TV>(particles,segment_mesh);
    bend->Set_Restlength_From_Particles();
    bend->Set_Stiffness(stiffness);
    bend->Set_Overdamping_Fraction(overdamping_fraction);
    bend->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    bend->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(restlength_enlargement_fraction) bend->Clamp_Restlength_With_Fraction_Of_Springs(restlength_enlargement_fraction);
    if(verbose) bend->Print_Restlength_Statistics();
    return bend;
}
//#####################################################################
// Function Create_Segment_Bending_Springs
//#####################################################################
template<class TV> SEGMENT_BENDING_SPRINGS<TV>* PhysBAM::
Create_Segment_Bending_Springs(SEGMENTED_CURVE<TV>& segmented_curve,const typename TV::SCALAR stiffness,const typename TV::SCALAR overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const typename TV::SCALAR max_strain_per_time_step,const bool use_rest_state_for_strain_rate,
    const typename TV::SCALAR restlength_enlargement_fraction,const bool verbose)
{
    return Create_Segment_Bending_Springs(dynamic_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles),segmented_curve.mesh,stiffness,overdamping_fraction,limit_time_step_by_strain_rate,max_strain_per_time_step,
        use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}
//#####################################################################
namespace PhysBAM{
template class SEGMENT_BENDING_SPRINGS<VECTOR<float,2> >;
template class SEGMENT_BENDING_SPRINGS<VECTOR<float,3> >;
template SEGMENT_BENDING_SPRINGS<VECTOR<float,2> >* Create_Segment_Bending_Springs<VECTOR<float,2> >(SEGMENTED_CURVE<VECTOR<float,2> >&,float,float,bool,float,bool,
    float,bool);
template SEGMENT_BENDING_SPRINGS<VECTOR<float,3> >* Create_Segment_Bending_Springs<VECTOR<float,3> >(SEGMENTED_CURVE<VECTOR<float,3> >&,float,float,bool,float,bool,
    float,bool);
template class SEGMENT_BENDING_SPRINGS<VECTOR<double,2> >;
template class SEGMENT_BENDING_SPRINGS<VECTOR<double,3> >;
template SEGMENT_BENDING_SPRINGS<VECTOR<double,2> >* Create_Segment_Bending_Springs<VECTOR<double,2> >(SEGMENTED_CURVE<VECTOR<double,2> >&,double,double,bool,double,bool,
    double,bool);
template SEGMENT_BENDING_SPRINGS<VECTOR<double,3> >* Create_Segment_Bending_Springs<VECTOR<double,3> >(SEGMENTED_CURVE<VECTOR<double,3> >&,double,double,bool,double,bool,
    double,bool);
}
