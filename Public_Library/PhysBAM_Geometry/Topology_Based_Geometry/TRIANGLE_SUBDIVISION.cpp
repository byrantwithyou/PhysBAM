//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eilene Hao, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
namespace PhysBAM{
//#####################################################################
// Function Clean_Memory
//#####################################################################
void TRIANGLE_SUBDIVISION::
Clean_Memory()
{
    if(delete_segment_mesh){delete triangle_mesh.segment_mesh;triangle_mesh.segment_mesh=0;}
    if(delete_neighbor_nodes){delete triangle_mesh.neighbor_nodes;triangle_mesh.neighbor_nodes=0;}
    if(delete_topologically_sorted_neighbor_nodes){delete triangle_mesh.topologically_sorted_neighbor_nodes;triangle_mesh.topologically_sorted_neighbor_nodes=0;}
    if(delete_boundary_mesh){delete triangle_mesh.boundary_mesh;triangle_mesh.boundary_mesh=0;}
    if(delete_incident_elements){delete triangle_mesh.incident_elements;triangle_mesh.incident_elements=0;}
    delete_segment_mesh=delete_neighbor_nodes=delete_topologically_sorted_neighbor_nodes=delete_boundary_mesh=delete_incident_elements=false;
}
//#####################################################################
// Function Refine_Mesh
//#####################################################################
void TRIANGLE_SUBDIVISION::
Refine_Mesh(TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input)
{
    start_index_for_new_nodes=start_index_for_new_nodes_input?start_index_for_new_nodes_input:triangle_mesh.number_nodes+1;
    if(!triangle_mesh.segment_mesh) delete_segment_mesh=true; // make it below, but delete it when this class is deleted
    bool triangle_edges_defined=triangle_mesh.element_edges!=0;
    if(!triangle_edges_defined)triangle_mesh.Initialize_Element_Edges(); // this makes a segment mesh as well

    refined_triangle_mesh.elements.Exact_Resize(4*triangle_mesh.elements.m);
    
    int new_t=0;
    for(int t=0;t<triangle_mesh.elements.m;t++){
        int i,j,k;triangle_mesh.elements(t).Get(i,j,k);
        int edge_ij,edge_jk,edge_ki;(*triangle_mesh.element_edges)(t).Get(edge_ij,edge_jk,edge_ki);
        int ij=start_index_for_new_nodes-1+edge_ij,jk=start_index_for_new_nodes-1+edge_jk,ki=start_index_for_new_nodes-1+edge_ki;
        refined_triangle_mesh.elements(new_t++).Set(i,ij,ki);
        refined_triangle_mesh.elements(new_t++).Set(j,jk,ij);
        refined_triangle_mesh.elements(new_t++).Set(k,ki,jk);
        refined_triangle_mesh.elements(new_t++).Set(jk,ki,ij);
        refined_triangle_mesh.number_nodes=max(refined_triangle_mesh.number_nodes,i,j,k,ij,jk,ki);} // update the number of nodes
    
    if(!triangle_edges_defined){delete triangle_mesh.element_edges;triangle_mesh.element_edges=0;} 
    // we keep segment_mesh for later use, i.e. it's not deleted!
}
//#####################################################################
// Function Refine_Mesh_Dual
//#####################################################################
void TRIANGLE_SUBDIVISION::
Refine_Mesh_Dual(TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input)
{
    start_index_for_new_nodes=start_index_for_new_nodes_input?start_index_for_new_nodes_input:triangle_mesh.number_nodes+1;
    bool adjacent_elements_defined=triangle_mesh.adjacent_elements!=0;if(!adjacent_elements_defined)triangle_mesh.Initialize_Adjacent_Elements();

    refined_triangle_mesh.elements.Exact_Resize(3*triangle_mesh.elements.m);
    
    int new_t=0;
    for(int t=0;t<triangle_mesh.elements.m;t++){
        assert((*triangle_mesh.adjacent_elements)(t).m==3);  // boundary subdivision unimplemented
        int tv=start_index_for_new_nodes-1+t;
        int ti,tj,tk;triangle_mesh.elements(t).Get(ti,tj,tk);
        for(int a=0;a<3;a++){
            int s=(*triangle_mesh.adjacent_elements)(t)(a);if(s<t)continue;
            int si,sj,sk;triangle_mesh.elements(s).Get(si,sj,sk);
            if(ti==si||ti==sj||ti==sk){cyclic_shift(ti,tj,tk);if(ti==si||ti==sj||ti==sk)cyclic_shift(ti,tj,tk);}
            int sv=start_index_for_new_nodes-1+s;
            refined_triangle_mesh.elements(new_t++).Set(tv,tj,sv);
            refined_triangle_mesh.elements(new_t++).Set(tv,sv,tk);}}
    refined_triangle_mesh.number_nodes=triangle_mesh.number_nodes+triangle_mesh.elements.m;
    
    if(!adjacent_elements_defined){delete triangle_mesh.adjacent_elements;triangle_mesh.adjacent_elements=0;} 
}
//#####################################################################
// Function Apply_Linear_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Linear_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==start_index_for_new_nodes-1+triangle_mesh.segment_mesh->elements.m);
    // neighbor_nodes are used to identify which nodes are in the base mesh - construct here and delete when this class is deleted
    if(!triangle_mesh.neighbor_nodes){delete_neighbor_nodes=true;triangle_mesh.Initialize_Neighbor_Nodes();}
    for(int i=0;i<triangle_mesh.number_nodes;i++)if((*triangle_mesh.neighbor_nodes)(i).m) subdivided_values(i)=base_values(i);
    // interpolate values on edges
    for(int k=0;k<triangle_mesh.segment_mesh->elements.m;k++)
        subdivided_values(start_index_for_new_nodes-1+k)=(T).5*(base_values(triangle_mesh.segment_mesh->elements(k)(0))+base_values(triangle_mesh.segment_mesh->elements(k)(1)));
}
//#####################################################################
// Function Apply_Fractal_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Fractal_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values,const float power)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==start_index_for_new_nodes-1+triangle_mesh.segment_mesh->elements.m);
    // neighbor_nodes are used to identify which nodes are in the base mesh - construct here and delete when this class is deleted
    if(!triangle_mesh.neighbor_nodes){delete_neighbor_nodes=true;triangle_mesh.Initialize_Neighbor_Nodes();}
    if(!triangle_mesh.incident_elements){delete_incident_elements=true;triangle_mesh.Initialize_Incident_Elements();}
    for(int i=0;i<triangle_mesh.number_nodes;i++)if((*triangle_mesh.neighbor_nodes)(i).m) subdivided_values(i)=base_values(i);
    // interpolate values on edges
    RANDOM_NUMBERS<T> random;
    for(int k=0;k<triangle_mesh.segment_mesh->elements.m;k++){
        int node1=triangle_mesh.segment_mesh->elements(k)(0);
        int node2=triangle_mesh.segment_mesh->elements(k)(1);
        TV midpoint=(T).5*(base_values(node1)+base_values(node2));
        TV offset=base_values(node1)-base_values(node2);
        T modulus=(T)power*random.Get_Gaussian();
        ARRAY<int> adjacent_triangles;
        triangle_mesh.Triangles_On_Edge(node1,node2,&adjacent_triangles);
        TV normal;
        for(int i=0;i<adjacent_triangles.m;i++){
            int n1,n2,n3;triangle_mesh.elements(adjacent_triangles(i)).Get(n1,n2,n3);
            normal+=PLANE<T>::Normal(base_values(n1),base_values(n2),base_values(n3));}
        normal.Normalize();
        //subdivided_values(start_index_for_new_nodes-1+k)=midpoint+normal*edge_magnitude*modulus;
        subdivided_values(start_index_for_new_nodes-1+k)=midpoint+normal*modulus;}
}
//#####################################################################
// Function Apply_Loop_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Loop_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==start_index_for_new_nodes-1+triangle_mesh.segment_mesh->elements.m);
    if(!triangle_mesh.topologically_sorted_neighbor_nodes){delete_topologically_sorted_neighbor_nodes=true;triangle_mesh.Initialize_Topologically_Sorted_Neighbor_Nodes();}
    if(!triangle_mesh.boundary_mesh){delete_boundary_mesh=true;triangle_mesh.Initialize_Boundary_Mesh();}
    if(!triangle_mesh.boundary_mesh->neighbor_nodes) triangle_mesh.boundary_mesh->Initialize_Neighbor_Nodes();
    ARRAY<ARRAY<int> > &neighbors=*triangle_mesh.topologically_sorted_neighbor_nodes,&boundary_neighbors=*triangle_mesh.boundary_mesh->neighbor_nodes;
    // node values
    for(int i=0;i<triangle_mesh.number_nodes;i++)if(neighbors(i).m){
        if(Node_Is_A_Corner(i) || boundary_neighbors(i).m>2 || boundary_neighbors(i).m==1)subdivided_values(i)=base_values(i);
        else if(boundary_neighbors(i).m==2) // if this is a regular boundary node
            subdivided_values(i)=(T).75*base_values(i)+(T).125*(base_values(boundary_neighbors(i)(0))+base_values(boundary_neighbors(i)(1)));
        else{ // interior node
            T alpha;switch(neighbors(i).m){
                case 3:alpha=(T).4375;break;case 4:alpha=(T).5;break;case 5:alpha=(T).54546609462891;break;case 6:alpha=(T).625;break;
                default:{T lambda=(T).375+(T).25*cos(T(2*pi)/neighbors(i).m);alpha=1-lambda*(4+lambda*(5*lambda-8))/(2*(1-lambda))+sqr(lambda);}}
            VECTOR<T,3> neighbor_sum=base_values(neighbors(i)(0));for(int j=1;j<neighbors(i).m;j++)neighbor_sum+=base_values(neighbors(i)(j));
            subdivided_values(i)=alpha*base_values(i)+(1-alpha)/neighbors(i).m*neighbor_sum;}}
    // edge values
    for(int i=0;i<triangle_mesh.segment_mesh->elements.m;i++){
        int index=start_index_for_new_nodes-1+i;
        int j,end1,end2;triangle_mesh.segment_mesh->elements(i).Get(end1,end2);
        if(boundary_neighbors(end1).m && boundary_neighbors(end2).m && boundary_neighbors(end1).Find(end2,j)) // if boundary edge
            subdivided_values(index)=(T).5*(base_values(end1)+base_values(end2));
        else if(neighbors(end1).m==6 && neighbors(end2).m==6){ // if next to regular vertices (the most common situation)
            j=0;neighbors(end1).Find(end2,j);assert(j);
            int common1=neighbors(end1)(j==0?neighbors(end1).m:j-1),common2=neighbors(end1)(j==neighbors(end1).m?1:j+1);
            subdivided_values(index)=(T).375*(base_values(end1)+base_values(end2))+(T).125*(base_values(common1)+base_values(common2));}
        else{
            if(neighbors(end1).m != 6){
                T beta;switch(neighbors(end1).m){
                    case 3:beta=(T).625;break;case 4:beta=(T).640625;break;case 5:beta=(T).65906781074217;break;
                    default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end1).m);beta=lambda*(4+lambda*(5*lambda-8))/(2*(1-lambda));}}
                subdivided_values(index)=(1-beta)*base_values(end1);
                int loop_start=0;neighbors(end1).Find(end2,loop_start);assert(loop_start);
                int j=loop_start;
                do{
                    T u=cos((T)(2*pi)*(j-loop_start)/neighbors(end1).m);
                    T w;switch(neighbors(end1).m){
                        case 3:w=(T)one_sixth*((T)1.25+u);break;case 4:w=(T).5*sqr((T).5+(T).375*u);break;
                        case 5:w=(T).16362712429686842801278667714785*sqr((T).55278640450004206071816526625374+u);break;
                        default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end1).m);
                            w=2*cube(lambda)/(neighbors(end1).m*(1-lambda))*(1+u)*sqr(1/lambda-(T)1.5+u);}}
                    subdivided_values(index)+=w*base_values(neighbors(end1)(j));
                    j++;if(j>neighbors(end1).m)j=1;}
                while(j!=loop_start);}
            if(neighbors(end2).m != 6){
                T beta;switch(neighbors(end2).m){
                    case 3:beta=(T).625;break;case 4:beta=(T).640625;break;case 5:beta=(T).65906781074217;break;
                    default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end2).m);beta=lambda*(4+lambda*(5*lambda-8))/(2*(1-lambda));}}
                if(neighbors(end1).m==6) subdivided_values(index)=(1-beta)*base_values(end2);
                else{subdivided_values(index)*=(T).5;subdivided_values(index)+=(T).5*(1-beta)*base_values(end2);} // need to average in this case
                int loop_start=0;neighbors(end2).Find(end1,loop_start);assert(loop_start);
                int j=loop_start;
                do{
                    T u=cos((T)(2*pi)*(j-loop_start)/neighbors(end2).m);
                    T w;switch(neighbors(end2).m){
                        case 3:w=(T)one_sixth*((T)1.25+u);break;case 4:w=(T).5*sqr((T).5+(T).375*u);break;
                        case 5:w=(T).16362712429686842801278667714785*sqr((T).55278640450004206071816526625374+u);break;
                        default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end2).m);
                            w=2*cube(lambda)/(neighbors(end2).m*(1-lambda))*(1+u)*sqr(1/lambda-(T)1.5+u);}}
                    if(neighbors(end1).m!=6)w*=(T).5; // need to average in this case
                    subdivided_values(index)+=w*base_values(neighbors(end2)(j));
                    j++;if(j>neighbors(end2).m)j=1;}
                while(j!=loop_start);}}}
}
//#####################################################################
// Function Apply_Root_Three_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Root_Three_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==start_index_for_new_nodes-1+triangle_mesh.elements.m);
    if(!triangle_mesh.neighbor_nodes){delete_neighbor_nodes=true;triangle_mesh.Initialize_Neighbor_Nodes();}
    for(int p=0;p<triangle_mesh.number_nodes;p++){
        int n=(*triangle_mesh.neighbor_nodes)(p).m;if(n<3)continue;T one_over_n=(T)1/n;
        T alpha=T(2./9)*(2-cos(T(2*pi)*one_over_n));
        TV sum=base_values((*triangle_mesh.neighbor_nodes)(p)(0));for(int i=1;i<n;i++)sum+=base_values((*triangle_mesh.neighbor_nodes)(p)(i));
        subdivided_values(p)=(1-alpha)*base_values(p)+alpha*one_over_n*sum;}
    for(int t=0;t<triangle_mesh.elements.m;t++){
        int i,j,k;triangle_mesh.elements(t).Get(i,j,k);
        subdivided_values(start_index_for_new_nodes-1+t)=(T)one_third*(base_values(i)+base_values(j)+base_values(k));}
}
//#####################################################################
template void TRIANGLE_SUBDIVISION::Apply_Linear_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Fractal_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values,const float power);
template void TRIANGLE_SUBDIVISION::Apply_Loop_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Root_Three_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void TRIANGLE_SUBDIVISION::Apply_Linear_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Fractal_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values,const float power);
template void TRIANGLE_SUBDIVISION::Apply_Loop_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Root_Three_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values);
#endif
}
