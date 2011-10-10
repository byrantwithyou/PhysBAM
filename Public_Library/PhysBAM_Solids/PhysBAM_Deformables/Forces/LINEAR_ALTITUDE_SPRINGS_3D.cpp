//#####################################################################
// Copyright 2002-2010, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>::
LINEAR_ALTITUDE_SPRINGS_3D(PARTICLES<TV>& particles,TETRAHEDRON_MESH& tetrahedron_mesh)
    :LINEAR_ALTITUDE_SPRINGS<VECTOR<T,3>,3>(particles,tetrahedron_mesh)
{
    use_shortest_spring_only=true;
}
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>::
~LINEAR_ALTITUDE_SPRINGS_3D()
{
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient) // assumes mass and restlength are already defined
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        for(int s=1;s<=4;s++){int node1,node2,node3,node4;
            Fill_Node_Indices(i,j,k,l,s,node1,node2,node3,node4);
            T one_over_triangle_mass=(T)one_ninth*(particles.one_over_effective_mass(node2)+particles.one_over_effective_mass(node3)+particles.one_over_effective_mass(node4)),
                  one_over_particle_mass=particles.one_over_effective_mass(node1),harmonic_mass=Pseudo_Inverse(one_over_particle_mass+one_over_triangle_mass);
            parameters(t)(s).youngs_modulus=scaling_coefficient*harmonic_mass/parameters(t)(s).restlength;}}
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Restlength_From_Particles()
{
    Set_Restlength_From_Material_Coordinates(particles.X);
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> material_coordinates)
{
    Invalidate_CFL();
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        // TODO: use a for loop instead
        TV barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(i),material_coordinates(j),material_coordinates(l),material_coordinates(k));
        parameters(t)(1).restlength=(material_coordinates(i)-barycentric.x*material_coordinates(j)-barycentric.y*material_coordinates(l)-barycentric.z*material_coordinates(k)).Normalize();
        barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(j),material_coordinates(i),material_coordinates(k),material_coordinates(l));
        parameters(t)(2).restlength=(material_coordinates(j)-barycentric.x*material_coordinates(i)-barycentric.y*material_coordinates(k)-barycentric.z*material_coordinates(l)).Normalize();
        barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(k),material_coordinates(i),material_coordinates(l),material_coordinates(j));
        parameters(t)(3).restlength=(material_coordinates(k)-barycentric.x*material_coordinates(i)-barycentric.y*material_coordinates(l)-barycentric.z*material_coordinates(j)).Normalize();
        barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(l),material_coordinates(i),material_coordinates(j),material_coordinates(k));
        parameters(t)(4).restlength=(material_coordinates(l)-barycentric.x*material_coordinates(i)-barycentric.y*material_coordinates(j)-barycentric.z*material_coordinates(k)).Normalize();
        assert(parameters(t)(1).restlength>0);assert(parameters(t)(2).restlength>0);assert(parameters(t)(3).restlength>0);assert(parameters(t)(4).restlength>0);
        for(int k=1;k<=4;k++) parameters(t)(k).visual_restlength=parameters(t)(k).restlength;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        // TODO: use a for loop instead
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(i)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(1).damping=overdamping_fraction*2*sqrt(parameters(t)(1).youngs_modulus*parameters(t)(1).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(j)+(T)one_ninth*(particles.one_over_effective_mass(i)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(2).damping=overdamping_fraction*2*sqrt(parameters(t)(2).youngs_modulus*parameters(t)(2).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(k)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(i)+particles.one_over_effective_mass(l)));
        parameters(t)(3).damping=overdamping_fraction*2*sqrt(parameters(t)(3).youngs_modulus*parameters(t)(3).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(l)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(i)));
        parameters(t)(4).damping=overdamping_fraction*2*sqrt(parameters(t)(4).youngs_modulus*parameters(t)(4).restlength*harmonic_mass);}
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Overdamping_Fraction(const ARRAY<VECTOR<T,4> >& overdamping_fraction) // 1 is critically damped
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        // TODO: use a for loop instead
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(i)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(1).damping=overdamping_fraction(t)(1)*2*sqrt(parameters(t)(1).youngs_modulus*parameters(t)(1).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(j)+(T)one_ninth*(particles.one_over_effective_mass(i)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(2).damping=overdamping_fraction(t)(2)*2*sqrt(parameters(t)(2).youngs_modulus*parameters(t)(2).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(k)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(i)+particles.one_over_effective_mass(l)));
        parameters(t)(3).damping=overdamping_fraction(t)(3)*2*sqrt(parameters(t)(3).youngs_modulus*parameters(t)(3).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(l)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(i)));
        parameters(t)(4).damping=overdamping_fraction(t)(4)*2*sqrt(parameters(t)(4).youngs_modulus*parameters(t)(4).restlength*harmonic_mass);}
    Invalidate_CFL();
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    ARRAY<VECTOR<T,4> > save_damping(parameters.m);for(int i=1;i<=parameters.m;i++) for(int k=1;k<=4;k++) save_damping(i)(k)=parameters(i)(k).damping;
    Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=1;k<=parameters.m;k++) for(int i=1;i<=4;i++) parameters(k)(i).damping=max(parameters(k)(i).damping,save_damping(k)(i));
}
//#####################################################################
// Function Fill_Node_Indices
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Fill_Node_Indices(int i,int j,int k,int l,int isolated_node_number,int& node1,int& node2,int& node3,int& node4) const
{
    // node1 is the isolated vertex and nodes2,3,4 are the triangle
    switch(isolated_node_number){
        case 1:node1=i;node2=j;node3=l;node4=k;break;
        case 2:node1=j;node2=i;node3=k;node4=l;break;
        case 3:node1=k,node2=i,node3=l,node4=j;break;
        default:node1=l;node2=i;node3=j;node4=k;}
}
//#####################################################################
// Function Fill_Spring_State
//#####################################################################
template<class T> bool LINEAR_ALTITUDE_SPRINGS_3D<T>::
Fill_Spring_State(int t,int isolated_node_number,int node1,int node2,int node3,int node4,SPRING_STATE& state)
{
    ARRAY_VIEW<const TV> X=particles.X;
    state.barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3),X(node4));
    state.direction=X(node1)-state.barycentric.x*X(node2)-state.barycentric.y*X(node3)-state.barycentric.z*X(node4);
    state.current_length=state.direction.Normalize();
    TV normal=TRIANGLE_3D<T>::Normal(X(node2),X(node3),X(node4));
    if(!state.current_length) state.direction=normal;
    else if(TV::Dot_Product(state.direction,normal)<0){
        state.direction=-state.direction;
        state.current_length=-state.current_length;}
    T rl=parameters(t)(isolated_node_number).restlength;
    if(use_springs_compressed_beyond_threshold && state.current_length-rl>-spring_compression_fraction_threshold*rl)
    {state.node=0;return false;}
    else{
        state.coefficient=parameters(t)(isolated_node_number).damping/rl;
        if(use_plasticity) Compute_Plasticity(isolated_node_number,t,state.current_length);
        state.node=isolated_node_number;
        return true;}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    if(use_shortest_spring_only) spring_states.Resize(mesh.elements.m);
    else spring_states_all_springs.Resize(mesh.elements.m);
    int used_springs=0,total_elements=0;
    ARRAY_VIEW<const TV> X=particles.X;int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data(); // use shortest spring only
        total_elements++;
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        if(use_shortest_spring_only){
            int hmin=0;T cross_area_max=-(T)FLT_MAX;
            if(is_position_update){
                for(int h=1;h<=4;h++){
                    Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                    T cross_area=TV::Cross_Product(X(node3)-X(node2),X(node4)-X(node2)).Magnitude_Squared();
                    if(cross_area>cross_area_max){hmin=h;cross_area_max=cross_area;}}}
            else hmin=spring_states(t).node;
            Fill_Node_Indices(i,j,k,l,hmin,node1,node2,node3,node4);
            SPRING_STATE& state=spring_states(t);
            if(Fill_Spring_State(t,hmin,node1,node2,node3,node4,state)) used_springs++;}
        else{
            for(int h=1;h<=4;h++){
                Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                SPRING_STATE& state=spring_states_all_springs(t)[h];
                if(Fill_Spring_State(t,h,node1,node2,node3,node4,state)) used_springs++;}}}

    residual_PE.Resize(mesh.elements.m);
    if(!total_PE.m){
        total_PE.Resize(mesh.elements.m);
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            total_PE(t)=Potential_Energy(t,time);}}

    if(print_number_used) LOG::cout<<"using "<<used_springs<<" of "<<total_elements<<" altitude springs"<<std::endl;
    if(compute_half_forces){
        if(use_shortest_spring_only) for(int i=1;i<=spring_states.m;i++){
            spring_states(i).sqrt_coefficient=sqrt(spring_states(i).coefficient);}
        else for(int i=1;i<=spring_states_all_springs.m;i++) for(int h=1;h<=4;h++){
                spring_states_all_springs(i)(h).sqrt_coefficient=sqrt(spring_states_all_springs(i)(h).coefficient);}}
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                const SPRING_PARAMETER& parameter=parameters(t)(state.node);
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                T rl=parameter.restlength,vrl=parameter.visual_restlength;
                TV force=parameter.youngs_modulus/rl*(state.current_length-vrl)*state.direction;
                F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;F(node4)+=state.barycentric.z*force;}}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                TV force=(state.coefficient*TV::Dot_Product(V(node1)-state.barycentric.x*V(node2)-state.barycentric.y*V(node3)-state.barycentric.z*V(node4),
                        state.direction))*state.direction;
                F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;F(node4)+=state.barycentric.z*force;}}}
}
//#####################################################################
// Function Velocity_Dependent_Size
//#####################################################################
template<class T> int LINEAR_ALTITUDE_SPRINGS_3D<T>::
Velocity_Dependent_Forces_Size() const
{
    int aggregate_id=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node) aggregate_id++;}}
    return aggregate_id;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    int aggregate_id=1;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                aggregate(aggregate_id++)+=state.sqrt_coefficient*TV::Dot_Product(V(node1)-state.barycentric.x*V(node2)-state.barycentric.y*V(node3)-state.barycentric.z*V(node4),state.direction);}}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    int aggregate_id=1;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                TV force=state.sqrt_coefficient*aggregate(aggregate_id++)*state.direction;
                F(node1)+=force;F(node2)-=state.barycentric.x*force;F(node3)-=state.barycentric.y*force;F(node4)-=state.barycentric.z*force;}}}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                const SPRING_PARAMETER& parameter=parameters(t)(state.node);
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                T rl=parameter.restlength;
                TV force=parameter.youngs_modulus/rl*TV::Dot_Product(V(node1)-state.barycentric.x*V(node2)-state.barycentric.y*V(node3)-state.barycentric.z*V(node4),state.direction)*state.direction;
                F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;F(node4)+=state.barycentric.z*force;}}}
}
template<class T> bool LINEAR_ALTITUDE_SPRINGS_3D<T>::
Compute_Strain_Rate_And_Strain(int t,int isolated_node_number,int node1,int node2,int node3,int node4,T& strain_rate,T& strain) const
{
    ARRAY_VIEW<const TV> X(particles.X),V(particles.V);
    TV direction=TRIANGLE_3D<T>::Normal(X(node2),X(node3),X(node4));
    T rl=parameters(t)(isolated_node_number).restlength;
    T dx;if(use_rest_state_for_strain_rate) dx=rl;else dx=TV::Dot_Product(direction,X(node1)-X(node2));
    if(use_springs_compressed_beyond_threshold){
        T length;if(use_rest_state_for_strain_rate) length=TV::Dot_Product(direction,X(node1)-X(node2));else length=dx;
        if(length-rl>-spring_compression_fraction_threshold*rl) return false;}
    TV barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3),X(node4));
    TV v_interpolated=barycentric.x*V(node2)+barycentric.y*V(node3)+barycentric.z*V(node4);
    strain_rate=TV::Dot_Product(V(node1)-v_interpolated,direction)/dx;
    strain=(dx-parameters(t)(isolated_node_number).visual_restlength)/rl;
    return true;
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0;int node1,node2,node3,node4; // node1 is the vertex and nodes2,3,4 are the triangle
    T strain_rate,strain;
    ARRAY<VECTOR<int,4> >& elements=mesh.elements;ARRAY_VIEW<const TV> X(particles.X);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int i,j,k,l;elements(t).Get(i,j,k,l);
        if(use_shortest_spring_only){
            int hmin=0;T cross_area_max=(T)-FLT_MAX;
            for(int h=1;h<=4;h++){
                Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                T cross_area=TV::Cross_Product(X(node3)-X(node2),X(node4)-X(node2)).Magnitude_Squared();
                if(cross_area>cross_area_max){hmin=h;cross_area_max=cross_area;}}
            Fill_Node_Indices(i,j,k,l,hmin,node1,node2,node3,node4);
            if(!Compute_Strain_Rate_And_Strain(t,hmin,node1,node2,node3,node4,strain_rate,strain)) continue;
            max_strain_rate=max(max_strain_rate,abs(strain_rate));
            if(cache_strain){strains_of_spring(t)=VECTOR<T,2>(abs(strain_rate),abs(strain));}}
        else{
            for(int h=1;h<=4;h++){
                Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                if(!Compute_Strain_Rate_And_Strain(t,h,node1,node2,node3,node4,strain_rate,strain)) continue;
                Compute_Strain_Rate_And_Strain(t,h,node1,node2,node3,node4,strain_rate,strain);
                max_strain_rate=max(max_strain_rate,abs(strain_rate));
                if(cache_strain){strains_of_spring_all_springs(t)[h]=VECTOR<T,2>(abs(strain_rate),abs(strain));}}}}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void LINEAR_ALTITUDE_SPRINGS_3D<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
    if((use_shortest_spring_only && !spring_states.m) || (!use_shortest_spring_only && !spring_states_all_springs.m)) return;
    ARRAY_VIEW<const TV> X=particles.X;
    FORCE_DATA<TV> force_data;
    if(force_name.empty()) force_data.name="LINEAR_ALTITUDE_SPRINGS_3D";
    else force_data.name=force_name;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                const SPRING_PARAMETER& parameter=parameters(t)(state.node);
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);

                force_data.state=(state.current_length-parameter.visual_restlength)/parameter.visual_restlength;
                force_data.first_action_point=X(node1);
                force_data.second_action_point=state.barycentric.x*X(node2)+state.barycentric.y*X(node3)+state.barycentric.z*X(node4);
                force_data_list.Append(force_data);}}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Potential_Energy(const int t,const T time) const
{
    const SPRING_STATE& state=spring_states(t);
    if(state.node){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
        const SPRING_PARAMETER& parameter=parameters(t)(state.node);
        Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
        T rl=parameter.restlength,vrl=parameter.visual_restlength;

        ARRAY_VIEW<const TV> X=particles.X;
        VECTOR<T,3> barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3),X(node4));
        TV direction=X(node1)-barycentric.x*X(node2)-barycentric.y*X(node3)-barycentric.z*X(node4);
        T current_length=direction.Normalize();
        TV normal=TRIANGLE_3D<T>::Normal(X(node2),X(node3),X(node4));
        if(!current_length) direction=normal;
        else if(TV::Dot_Product(direction,normal)<0){
            direction=-direction;
            current_length=-current_length;}
        return (T).5*parameter.youngs_modulus/rl*sqr(current_length-vrl);}
    else return 0;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        potential_energy+=Potential_Energy(t,time);}
    return potential_energy;
}
//#####################################################################
// Function Residual_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Residual_Energy(const T time) const
{
    T residual_energy=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        residual_energy+=residual_PE(t);}
    return residual_energy;
}

template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>* PhysBAM::
Create_Altitude_Springs(PARTICLES<VECTOR<T,3> >& particles,TETRAHEDRON_MESH& mesh,
    const T stiffness,const T overdamping_fraction,const bool use_compressed_by_threshold_only,const T fraction_compression,const bool limit_time_step_by_strain_rate,
    const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const T restlength_enlargement_fraction,const bool verbose)
{
    return Create_Altitude_Springs_Base(particles,mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,fraction_compression,
        limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>* PhysBAM::
Create_Altitude_Springs(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,
    const T stiffness,const T overdamping_fraction,const bool use_compressed_by_threshold_only,const T fraction_compression,
    const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const T restlength_enlargement_fraction,
    const bool verbose)
{
    return Create_Altitude_Springs(dynamic_cast<PARTICLES<VECTOR<T,3> >&>(tetrahedralized_volume.particles),tetrahedralized_volume.mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,
        fraction_compression,limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}
//#####################################################################
template class LINEAR_ALTITUDE_SPRINGS_3D<float>;
template LINEAR_ALTITUDE_SPRINGS_3D<float>* PhysBAM::Create_Altitude_Springs<float>(TETRAHEDRALIZED_VOLUME<float>&,float,float,bool,float,bool,float,bool,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_ALTITUDE_SPRINGS_3D<double>;
template LINEAR_ALTITUDE_SPRINGS_3D<double>* PhysBAM::Create_Altitude_Springs<double>(TETRAHEDRALIZED_VOLUME<double>&,double,double,bool,double,bool,double,bool,double,bool);
#endif
