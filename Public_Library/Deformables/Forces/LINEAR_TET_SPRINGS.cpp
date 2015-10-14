//#####################################################################
// Copyright 2007, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/constants.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Particles/PARTICLES.h>
#include <Geometry/Basic_Geometry/PLANE.h>
#include <Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <fstream>
#include <iostream>
using namespace PhysBAM;
//#####################################################################
// LINEAR_TET_SPRINGS
//#####################################################################
template<class T> LINEAR_TET_SPRINGS<T>::
LINEAR_TET_SPRINGS(DEFORMABLE_PARTICLES<TV>& particles,TETRAHEDRON_MESH& mesh)
    :DEFORMABLES_FORCES<TV>(particles),mesh(mesh),spring_parameters(mesh.elements.Size(),false),
    edge_restlength_squared(mesh.elements.Size(),false),minimum_edge_compression_squared((T)1e-6),minimum_sin((T)pi/(T)180)
{
    Use_Springs_Compressed_Beyond_Threshold_Only(false);
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Set_Restlength_From_Particles()
{
    Set_Restlength_From_Material_Coordinates(particles.X);
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> material_coordinates)
{
    Invalidate_CFL();
    for(int t=0;t<mesh.elements.m;t++){
        const VECTOR<int,4>& element_nodes=mesh.elements(t);
        for(int s=0;s<4;s++){
            const VECTOR<int,4> spring_nodes=Spring_Nodes(s,element_nodes);
            SPRING_PARAMETER& param=spring_parameters(t)(s);
            param.restlength=param.visual_restlength=
                PLANE<T>(material_coordinates(spring_nodes[1]),material_coordinates(spring_nodes[2]),material_coordinates(spring_nodes[3])).Signed_Distance(material_coordinates(spring_nodes[0]));
            assert(param.restlength);}
        for(int s=4;s<7;s++){
            SPRING_PARAMETER& param=spring_parameters(t)(s);
            const VECTOR<int,4> spring_nodes=Spring_Nodes(s,element_nodes);
            SEGMENT_3D<T> segment_1(material_coordinates(spring_nodes[0]),material_coordinates(spring_nodes[1])),
                segment_2(material_coordinates(spring_nodes[2]),material_coordinates(spring_nodes[3]));
            VECTOR<T,2> weights;
            TV normal=segment_1.Shortest_Vector_Between_Lines(segment_2,weights);
            param.visual_restlength=param.restlength=normal.Magnitude();}
        int count=0;
        for(int i=0;i<3;i++) for(int j=i+1;j<4;j++) edge_restlength_squared(t)(count++)=(material_coordinates(element_nodes[j])-material_coordinates(element_nodes[i])).Magnitude_Squared();}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T LINEAR_TET_SPRINGS<T>::
CFL_Strain_Rate() const
{
    T dx=0,max_strain_rate=0;
    ARRAY_VIEW<const TV> V(particles.V);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const VECTOR<int,4>& element_nodes=mesh.elements(t);
        T minimum_signed_distance=FLT_MAX;TV minimum_normal;TV weights;VECTOR<int,4> spring_nodes;
        int min_spring_index=Find_Shortest_Spring(t,element_nodes,spring_nodes,minimum_signed_distance,minimum_normal,weights);
        if(min_spring_index==0) continue; // TODO: in coincident case, no spring present
        const SPRING_PARAMETER& params=spring_parameters(t)(min_spring_index);
        if(use_rest_state_for_strain_rate) dx=params.restlength;else dx=minimum_signed_distance;
        if(use_springs_compressed_beyond_threshold){
            T length;if(use_rest_state_for_strain_rate) length=minimum_signed_distance;else length=dx;
            if(length-params.restlength>-spring_compression_fraction_threshold*params.restlength) continue;}
        if(min_spring_index<5){
            T strain_rate=TV::Dot_Product(V(spring_nodes[0])-(weights.x*V(spring_nodes[1])+weights.y*V(spring_nodes[2])+weights.z*V(spring_nodes[3])),minimum_normal)/dx;
            max_strain_rate=max(max_strain_rate,abs(strain_rate));}
        else{
            T strain_rate=TV::Dot_Product((1-weights.x)*V(spring_nodes[0])+weights.x*V(spring_nodes[1])-
                (1-weights.y)*V(spring_nodes[2])-weights.y*V(spring_nodes[3]),minimum_normal)/dx;
            max_strain_rate=max(max_strain_rate,abs(strain_rate));}}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    // TODO: add CFL component for edge edges
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const VECTOR<int,4>& element_nodes=mesh.elements(t);
        FREQUENCY_DATA element_frequency;
        for(int s=0;s<4;s++){ // point face
            const VECTOR<int,4> spring_nodes=Spring_Nodes(s,element_nodes);const SPRING_PARAMETER& param=spring_parameters(t)(s);
            T one_over_mass_times_restlength=particles.one_over_effective_mass(spring_nodes[0]);
            for(int j=1;j<4;j++) one_over_mass_times_restlength+=((T)1/9)*particles.one_over_effective_mass(spring_nodes[j]);
            one_over_mass_times_restlength/=param.restlength;
            T elastic_squared=4*param.youngs_modulus*one_over_mass_times_restlength*one_over_cfl_number_squared;
            T damping=2*param.damping*one_over_mass_times_restlength*one_over_cfl_number;
            if(elastic_squared>element_frequency.elastic_squared) element_frequency.elastic_squared=elastic_squared;
            if(damping>element_frequency.damping) element_frequency.damping=damping;}
        for(int s=4;s<7;s++){ // edge edge
            const VECTOR<int,4> spring_nodes=Spring_Nodes(s,element_nodes);const SPRING_PARAMETER& param=spring_parameters(t)(s);
            T one_over_mass_times_restlength=(T).25*VECTOR<T,4>(particles.one_over_effective_mass.Subset(spring_nodes)).Sum()/param.restlength; // VECTOR<T,4> dodges a gcc 4.1.1 bug
            T elastic_squared=4*param.youngs_modulus*one_over_mass_times_restlength*one_over_cfl_number_squared;
            T damping=2*param.damping*one_over_mass_times_restlength*one_over_cfl_number;
            if(elastic_squared>element_frequency.elastic_squared) element_frequency.elastic_squared=elastic_squared;
            if(damping>element_frequency.damping) element_frequency.damping=damping;}
        for(int i=0;i<4;i++){
            frequency(element_nodes[i]).elastic_squared+=element_frequency.elastic_squared;
            frequency(element_nodes[i]).damping+=element_frequency.damping;}}
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_elements.Update(mesh.elements,particle_is_simulated);
    // TODO: implement me
    //if(cache_strain) strains_of_spring.Resize(mesh.elements.m,false,false);
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian)
{
    spring_states.Resize(mesh.elements.m);

    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const VECTOR<int,4>& element_nodes=mesh.elements(t);
        SPRING_STATE& state=spring_states(t);
        state.id=Find_Shortest_Spring(t,element_nodes,state.spring_nodes,state.current_length,state.direction,state.weights); // TODO: minimum weights needs to be computed
        if(state.id!=0){
            const SPRING_PARAMETER& params=spring_parameters(t)(state.id);
            state.coefficient=params.damping/params.restlength;}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& spring_state=spring_states(t);
        if(spring_state.id!=0){
            const SPRING_PARAMETER& spring_param=spring_parameters(t)(spring_state.id);
            const VECTOR<int,4>& spring_nodes=spring_state.spring_nodes;
            
            TV force=spring_param.youngs_modulus/spring_param.restlength*(spring_state.current_length-spring_param.visual_restlength)*spring_state.direction;
            if(spring_state.id<5){
                F(spring_nodes[0])-=force;
                F(spring_nodes[1])+=spring_state.weights.x*force;F(spring_nodes[2])+=spring_state.weights.y*force;F(spring_nodes[3])+=spring_state.weights.z*force;}
            else if(spring_state.id<8){
                F(spring_nodes[0])-=(1-spring_state.weights.x)*force;F(spring_nodes[1])-=spring_state.weights.x*force;
                F(spring_nodes[2])+=(1-spring_state.weights.y)*force;F(spring_nodes[3])+=spring_state.weights.y*force;}}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& spring_state=spring_states(t);
        if(spring_state.id!=0){

            const VECTOR<int,4>& spring_nodes=spring_state.spring_nodes;
            if(spring_state.id<5){
                TV force=(spring_state.coefficient*TV::Dot_Product(V(spring_nodes[0])
                        -spring_state.weights.x*V(spring_nodes[1])-spring_state.weights.y*V(spring_nodes[2])-spring_state.weights.z*V(spring_nodes[3]),spring_state.direction))*spring_state.direction;
                F(spring_nodes[0])-=force;
                F(spring_nodes[1])+=spring_state.weights.x*force;F(spring_nodes[2])+=spring_state.weights.y*force;F(spring_nodes[3])+=spring_state.weights.z*force;}
            else if(spring_state.id<8){
                TV force=(spring_state.coefficient*TV::Dot_Product((1-spring_state.weights.x)*V(spring_nodes[0])+spring_state.weights.x*V(spring_nodes[1])
                        -(1-spring_state.weights.y)*V(spring_nodes[2])-spring_state.weights.y*V(spring_nodes[3]),spring_state.direction))*spring_state.direction;
                F(spring_nodes[0])-=(1-spring_state.weights.x)*force;F(spring_nodes[1])-=spring_state.weights.x*force;
                F(spring_nodes[2])+=(1-spring_state.weights.y)*force;F(spring_nodes[3])+=spring_state.weights.y*force;}}}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& spring_state=spring_states(t);
        if(spring_state.id!=0){
            const SPRING_PARAMETER& spring_param=spring_parameters(t)(spring_state.id);
            const VECTOR<int,4>& spring_nodes=spring_state.spring_nodes;
            if(spring_state.id<5){
                TV force=(spring_param.youngs_modulus/spring_param.restlength*TV::Dot_Product(V(spring_nodes[0])
                        -spring_state.weights.x*V(spring_nodes[1])-spring_state.weights.y*V(spring_nodes[2])-spring_state.weights.z*V(spring_nodes[3]),spring_state.direction))*spring_state.direction;
                F(spring_nodes[0])-=force;
                F(spring_nodes[1])+=spring_state.weights.x*force;F(spring_nodes[2])+=spring_state.weights.y*force;F(spring_nodes[3])+=spring_state.weights.z*force;}
            else if(spring_state.id<8){
                TV force=(spring_param.youngs_modulus/spring_param.restlength*TV::Dot_Product((1-spring_state.weights.x)*V(spring_nodes[0])+spring_state.weights.x*V(spring_nodes[1])
                        -(1-spring_state.weights.y)*V(spring_nodes[2])-spring_state.weights.y*V(spring_nodes[3]),spring_state.direction))*spring_state.direction;
                F(spring_nodes[0])-=(1-spring_state.weights.x)*force;F(spring_nodes[1])-=spring_state.weights.x*force;
                F(spring_nodes[2])+=(1-spring_state.weights.y)*force;F(spring_nodes[3])+=spring_state.weights.y*force;}}}
}
//#####################################################################
// Function Use_Springs_Compressed_Beyond_Threshold_Only
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Use_Springs_Compressed_Beyond_Threshold_Only(const bool use_springs_compressed_beyond_threshold_input,const T threshold_fraction)
{
    use_springs_compressed_beyond_threshold=use_springs_compressed_beyond_threshold_input;spring_compression_fraction_threshold=threshold_fraction;
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    ARRAY_VIEW<T>& one_over_effective_mass=particles.one_over_effective_mass;
    for(int t=0;t<mesh.elements.m;t++){
        const VECTOR<int,4>& element_nodes=mesh.elements(t);
        for(int s=0;s<4;s++){
            VECTOR<int,4> spring_nodes=Spring_Nodes(s,element_nodes);
            SPRING_PARAMETER& spring_param=spring_parameters(t)(s);
            T harmonic_mass=Pseudo_Inverse(one_over_effective_mass(spring_nodes[0])
                +((T)1/9)*(one_over_effective_mass(spring_nodes[1])+one_over_effective_mass(spring_nodes[2])+one_over_effective_mass(spring_nodes[3])));
            spring_param.damping=overdamping_fraction*2*sqrt(spring_param.youngs_modulus*spring_param.restlength*harmonic_mass);}
        for(int s=4;s<7;s++){
            VECTOR<int,4> spring_nodes=Spring_Nodes(s,element_nodes);
            SPRING_PARAMETER& spring_param=spring_parameters(t)(s);
            T harmonic_mass=Pseudo_Inverse((T).25*(one_over_effective_mass(spring_nodes[0])+one_over_effective_mass(spring_nodes[1])
                    +one_over_effective_mass(spring_nodes[2])+one_over_effective_mass(spring_nodes[3])));
            spring_param.damping=overdamping_fraction*2*sqrt(spring_param.youngs_modulus*spring_param.restlength*harmonic_mass);}}
    Invalidate_CFL();
}
//#####################################################################
// Function Clamp_Restlength_With_Fraction_Of_Springs
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Clamp_Restlength_With_Fraction_Of_Springs(const T fraction)
{
    {ARRAY<T> length(spring_count*spring_parameters.m,false);
    for(int s=0;s<spring_parameters.m;s++) for(int k=0;k<spring_count;k++) length(spring_count*(s-1)+k)=spring_parameters(s)(k).restlength;
    length.Sort();
    T minimum_restlength=length(min((int)(fraction*length.m)+1,length.m));LOG::cout<<"Enlarging the restlength of all altitude springs below "<<minimum_restlength<<std::endl;
    for(int i=0;i<spring_parameters.m;i++) for(int k=0;k<spring_count;k++) spring_parameters(i)(k).restlength=max(minimum_restlength,spring_parameters(i)(k).restlength);}
    {ARRAY<T> edge_length(6*edge_restlength_squared.m,false);
    for(int s=0;s<edge_restlength_squared.m;s++) for(int k=0;k<6;k++) edge_length(6*(s-1)+k)=edge_restlength_squared(s)(k);
    edge_length.Sort();
    T minimum_edge_restlength=edge_length(min((int)(fraction*edge_length.m)+1,edge_length.m));
    LOG::cout<<"Enlarging the restlength of all altitude spring edges below "<<minimum_edge_restlength<<std::endl;
    for(int i=0;i<edge_restlength_squared.m;i++) for(int k=0;k<6;k++) edge_restlength_squared(i)(k)=max(minimum_edge_restlength,edge_restlength_squared(i)(k));}
}
//#####################################################################
// Function Print_Restlength_Statistics
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Print_Restlength_Statistics()
{
    LOG::cout<<"Tetrahedron Springs - Total Springs = "<<spring_count*spring_parameters.m<<std::endl;
    ARRAY<T> length(spring_count*spring_parameters.m,false),visual_restlength(spring_count*spring_parameters.m,false);
    for(int s=0;s<spring_parameters.m;s++) for(int k=0;k<spring_count;k++){
            length(spring_count*(s-1)+k)=spring_parameters(s)(k).restlength;visual_restlength(spring_count*(s-1)+k)=spring_parameters(s)(k).visual_restlength;}
    length.Sort();
    visual_restlength.Sort();
    int one_percent=(int)(.01*length.m),ten_percent=(int)(.1*length.m),median=(int)(.5*length.m);
    LOG::cout<<"Tetrahedron Springs - Smallest Restlength = "<<length(0)<<", Visual Restlength = "<<visual_restlength(0)<<std::endl;
    LOG::cout<<"Tetrahedron Springs - One Percent Restlength = "<<length(one_percent)<<", Visual Restlength = "<<visual_restlength(one_percent)<<std::endl;
    LOG::cout<<"Tetrahedron Springs - Ten Percent Restlength = "<<length(ten_percent)<<", Visual Restlength = "<<visual_restlength(ten_percent)<<std::endl;
    LOG::cout<<"Tetrahedron Springs - Median Restlength = "<<length(median)<<", Visual Restlength = "<<visual_restlength(median)<<std::endl;
}
//#####################################################################
// Function Find_Shortest_Spring
//#####################################################################
template<class T> int LINEAR_TET_SPRINGS<T>::
Find_Shortest_Spring(const int tet,const VECTOR<int,4> element_nodes,VECTOR<int,4>& spring_nodes_output,T& minimum_signed_distance,TV& minimum_normal,TV& weights) const
{
    ARRAY_VIEW<const TV> X(particles.X);
    int maximum_cross_squared_index=-1;T maximum_cross_squared=(T)-FLT_MAX;TV maximum_cross;
    for(unsigned char h=0;h<4;h++){
        VECTOR<int,4> spring_nodes=Spring_Nodes(h,element_nodes);
        TV u_cross_v=TV::Cross_Product(X(spring_nodes[2])-X(spring_nodes[1]),X(spring_nodes[3])-X(spring_nodes[1]));
        T u_cross_v_squared=u_cross_v.Magnitude_Squared();
        if(u_cross_v_squared>maximum_cross_squared){maximum_cross_squared_index=h;maximum_cross_squared=u_cross_v_squared;maximum_cross=u_cross_v;}}
    for(unsigned char h=4;h<7;h++){
        VECTOR<int,4> spring_nodes=Spring_Nodes(h,element_nodes);
        TV u_cross_v=TV::Cross_Product(X(spring_nodes[1])-X(spring_nodes[0]),X(spring_nodes[3])-X(spring_nodes[2]));
        T u_cross_v_squared=u_cross_v.Magnitude_Squared();
        if(u_cross_v_squared>maximum_cross_squared){maximum_cross_squared_index=h;maximum_cross_squared=u_cross_v_squared;maximum_cross=u_cross_v;}}
    // check robustness
    TV u,v;
    VECTOR<int,4> spring_nodes=Spring_Nodes(maximum_cross_squared_index,element_nodes);
    VECTOR<int,2> edge_indices=Edge_Indices(maximum_cross_squared_index);
    if(maximum_cross_squared_index<5){u=X(spring_nodes[2])-X(spring_nodes[1]);v=X(spring_nodes[3])-X(spring_nodes[1]);}
    else{u=X(spring_nodes[1])-X(spring_nodes[0]);v=X(spring_nodes[3])-X(spring_nodes[2]);}
    T u_length_squared=u.Magnitude_Squared(),v_length_squared=v.Magnitude_Squared();
    if(u_length_squared<minimum_edge_compression_squared*edge_restlength_squared(tet)(edge_indices[0])
        || v_length_squared<minimum_edge_compression_squared*edge_restlength_squared(tet)(edge_indices[1])) return 0;
    if(abs(maximum_cross_squared)<minimum_sin*u_length_squared*v_length_squared){
        VECTOR<int,2> edge_index;T max_distance_squared=0;
        for(int i=0;i<3;i++) for(int j=i+1;j<4;j++){
            T distance_squared=(X(element_nodes[i])-X(element_nodes[j])).Magnitude_Squared();
            if(distance_squared>max_distance_squared){edge_index=VECTOR<int,2>(i,j);max_distance_squared=distance_squared;}}
        VECTOR<int,2> other_nodes=element_nodes.Remove_Index(edge_index[1]).Remove_Index(edge_index[0]);
        if((X(element_nodes[edge_index[0]])-X(other_nodes[1])).Magnitude_Squared()<(X(element_nodes[edge_index[0]])-X(other_nodes[0])).Magnitude_Squared())
            other_nodes=other_nodes.Reversed();
        VECTOR<int,2> edge1(element_nodes[edge_index[0]],other_nodes[1]),edge2(other_nodes[0],element_nodes[edge_index[1]]);
        bool found=false;
        for(maximum_cross_squared_index=4;maximum_cross_squared_index<7;maximum_cross_squared_index++){
            spring_nodes=Spring_Nodes(maximum_cross_squared_index,element_nodes);
            if(spring_nodes.Slice<0,1>()==edge1 || spring_nodes.Slice<2,3>()==edge1 || spring_nodes.Slice<0,1>().Reversed()==edge1 || spring_nodes.Slice<2,3>().Reversed()==edge1){
                found=true;break;}}
        PHYSBAM_ASSERT(found);
        minimum_normal=(X(edge1[0])-X(edge2[1])).Orthogonal_Vector().Normalized();
        minimum_signed_distance=0;
        TV midpoint=(T).5*(X(edge1[1])+X(edge2[0]));
        weights=VECTOR<T,3>(SEGMENT_3D<T>(X(spring_nodes[0]),X(spring_nodes[1])).Interpolation_Fraction(midpoint),
            SEGMENT_3D<T>(X(spring_nodes[2]),X(spring_nodes[3])).Interpolation_Fraction(midpoint),0);}
    else{
        minimum_normal=maximum_cross.Normalized();
        minimum_signed_distance=TV::Dot_Product(minimum_normal,X(spring_nodes[0])-X(spring_nodes[2]));
        if(maximum_cross_squared_index<5){
            weights=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(spring_nodes[0]),X(spring_nodes[1]),X(spring_nodes[2]),X(spring_nodes[3]));}
        else if(maximum_cross_squared_index<8){
            VECTOR<T,2> dummy_weights;
            SEGMENT_3D<T>(X(spring_nodes[0]),X(spring_nodes[1])).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(X(spring_nodes[2]),X(spring_nodes[3])),dummy_weights);
            weights=dummy_weights.Append(0);}}
    // finish
    spring_nodes_output=spring_nodes;
    return maximum_cross_squared_index;
}
//#####################################################################
// Function Spring_Nodes
//#####################################################################
template<class T> VECTOR<int,4> LINEAR_TET_SPRINGS<T>::
Spring_Nodes(unsigned char pair_id,const VECTOR<int,4>& n)
{
    switch(pair_id){
        case 0: return VECTOR<int,4>(n[0],n[1],n[3],n[2]);case 1: return VECTOR<int,4>(n[1],n[0],n[2],n[3]); // point face
        case 2: return VECTOR<int,4>(n[2],n[0],n[3],n[1]);case 3: return VECTOR<int,4>(n[3],n[0],n[1],n[2]); // point face
        case 4: return VECTOR<int,4>(n[0],n[1],n[2],n[3]);case 5: return VECTOR<int,4>(n[1],n[2],n[0],n[3]);case 6: return VECTOR<int,4>(n[0],n[2],n[3],n[1]); // edge edge
        default: PHYSBAM_FATAL_ERROR();}
}
template<class T> VECTOR<int,2> LINEAR_TET_SPRINGS<T>::
Edge_Indices(unsigned char pair_id)
{
    switch(pair_id){
        case 0: return VECTOR<int,2>(4,3);case 1: return VECTOR<int,2>(1,2);case 2: return VECTOR<int,2>(2,0);case 3: return VECTOR<int,2>(0,1);
        case 4: return VECTOR<int,2>(0,5);case 5: return VECTOR<int,2>(3,2);case 6: return VECTOR<int,2>(1,4);
        default: PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Create_Tet_Springs
//#####################################################################
template<class T> LINEAR_TET_SPRINGS<T>* PhysBAM::
Create_Tet_Springs(DEFORMABLE_PARTICLES<VECTOR<T,3> >& particles,TETRAHEDRON_MESH& mesh,const T stiffness,const T overdamping_fraction,
    const bool use_compressed_by_threshold_only/*=true*/,const T fraction_compression/*=(T).1*/,const bool limit_time_step_by_strain_rate/*=true*/,const T max_strain_per_time_step/*=(T).1*/,
    const bool use_rest_state_for_strain_rate/*=true*/,const T restlength_enlargement_fraction/*=0*/,const bool verbose/*=true*/)
{
    LINEAR_TET_SPRINGS<T>* lts=new LINEAR_TET_SPRINGS<T>(particles,mesh);
    lts->Set_Restlength_From_Particles();
    if(restlength_enlargement_fraction) lts->Clamp_Restlength_With_Fraction_Of_Springs(restlength_enlargement_fraction);
    lts->Set_Stiffness(stiffness);
    lts->Set_Overdamping_Fraction(overdamping_fraction);
    lts->Use_Springs_Compressed_Beyond_Threshold_Only(use_compressed_by_threshold_only,fraction_compression);
    lts->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    lts->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(verbose) lts->Print_Restlength_Statistics();
    return lts;
}
//#####################################################################
// Function Create_Tet_Springs
//#####################################################################
template<class T> LINEAR_TET_SPRINGS<T>* PhysBAM::
Create_Tet_Springs(TETRAHEDRALIZED_VOLUME<T>& volume,const T stiffness,
    const T overdamping_fraction,const bool use_compressed_by_threshold_only/*=true*/,const T fraction_compression/*=(T).1*/,const bool limit_time_step_by_strain_rate/*=true*/,
    const T max_strain_per_time_step/*=(T).1*/,const bool use_rest_state_for_strain_rate/*=true*/,const T restlength_enlargement_fraction/*=0*/,const bool verbose/*=true*/)
{
    return Create_Tet_Springs(dynamic_cast<DEFORMABLE_PARTICLES<VECTOR<T,3> >&>(volume.particles),volume.mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,fraction_compression,limit_time_step_by_strain_rate,
        max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Set_Stiffness(const T youngs_modulus_input)
{
    Invalidate_CFL();
    for(int i=0;i<spring_parameters.m;i++)
        for(int s=0;s<spring_count;s++)
            spring_parameters(i)(s).youngs_modulus=youngs_modulus_input;
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Set_Stiffness(const ARRAY<VECTOR<T,spring_count> >& youngs_modulus_input)
{
    Invalidate_CFL();
    for(int i=0;i<spring_parameters.m;i++)
        for(int s=0;s<spring_count;s++)
            spring_parameters(i)(s).youngs_modulus=youngs_modulus_input(i)(s);
}
//#####################################################################
// Function Set_Restlength
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Set_Restlength(const ARRAY<VECTOR<T,spring_count> >& restlength_input)
{Invalidate_CFL();
    for(int i=0;i<spring_parameters.m;i++)
        for(int s=0;s<spring_count;s++)
            spring_parameters(i)(s).restlength=restlength_input(i)(s);
}
//#####################################################################
// Function Clamp_Restlength
//#####################################################################
template<class T> void LINEAR_TET_SPRINGS<T>::
Clamp_Restlength(const T clamped_restlength)
{Invalidate_CFL();
    for(int i=0;i<spring_parameters.m;i++)
        for(int s=0;s<spring_count;s++)
            spring_parameters(i)(s).restlength=max(spring_parameters(i)(s).visual_restlength,clamped_restlength);
}
//#####################################################################
namespace PhysBAM{
template LINEAR_TET_SPRINGS<float>* Create_Tet_Springs<float>(TETRAHEDRALIZED_VOLUME<float>&,float,float,bool,float,bool,float,bool,float,bool);
template class LINEAR_TET_SPRINGS<float>;
template LINEAR_TET_SPRINGS<double>* Create_Tet_Springs<double>(TETRAHEDRALIZED_VOLUME<double>&,double,double,bool,double,bool,double,bool,double,bool);
template class LINEAR_TET_SPRINGS<double>;
}
