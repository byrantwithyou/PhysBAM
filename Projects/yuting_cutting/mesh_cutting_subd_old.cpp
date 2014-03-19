//
//  mesh_cutting.cpp
//
//
//  Created by Yuting Wang on 4/5/12.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.
//

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry_Intersections/SEGMENT_3D_TRIANGLE_3D_INTERSECTION.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>

#include <iostream>
#include "mesh_cutting_subd_old.h"
#include "DEFORMABLE_OBJECTS.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <set>

//need modiry:
//Translate_cc: translate sim_mesh first, than update cutting particles.
//split and cut
//dragging
    
    
void read_tsc(){__asm__("rdtsc");}
struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);read_tsc();}
void stop_timer(){gettimeofday(&stoptime,NULL);read_tsc();}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}
#define DE cout<<"file "<<__FILE__<<"   line "<<__LINE__<<"  "<<&volume->particles.X<<"   "<<volume->particles.X<<endl;

using namespace PhysBAM;
using namespace std;

//won't use memory
static const int NumNodesPerTriangle = 3;
static const int NumEdgesPerTriangle = 3;
static const int NumFacesPerTet = 4;
static const int NumNodesPerTet = 4;
static const int NumEdgesPerTet = 6;
static const int NumMaterialsPerFace = 6;
static const int NumSplitTurnOns = 36;
static const int MaxPieces = 24;
    
//index transformations
int face2node[NumFacesPerTet][NumNodesPerTriangle] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};

int face2edge[NumFacesPerTet][NumEdgesPerTriangle] = {{0, 1, 2}, {0, 4, 3}, {1, 5, 4}, {2, 3, 5}};

int edge2node[NumEdgesPerTet][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 0}, {3, 1}, {3, 2}};
int edge2opposite_node[NumEdgesPerTet][2] = {{2, 3}, {0, 3}, {1, 3}, {1, 2}, {0, 2}, {0, 1}};

int material2face[MaxPieces][NumFacesPerTet-1] = {{0,1,24}, {1,2,25}, {2,3,26},{3,4,27},{4,5,28},{5,0,29},
                                                  {6,7,24},{7,8,25},{8,9,32},{9,10,33},{10,11,31},{11,6,30},
                                                  {12,13,26},{13,14,27},{14,15,34},{15,16,35},{16,17,33},{17,12,32},
                                                  {18,19,28},{19,20,29},{20,21,30},{21,22,31},{22,23,35},{23,18,34}};

int material2material[MaxPieces][NumFacesPerTet-1] = {{5,1,6}, {0,2,7}, {1,3,12},{2,4,13},{3,5,18},{4,0,19},
                                                      {11,7,0},{6,8,1},{7,9,17},{8,10,16},{9,11,21},{10,6,20},
                                                      {17,13,2},{12,14,3},{13,15,23},{14,16,22},{15,17,9},{16,12,8},
                                                      {23,19,4},{18,20,5},{19,21,11},{20,22,10},{21,23,15},{22,18,14}};

int material2edge[MaxPieces] = {0, 0, 1, 1, 2, 2, 
                                0, 0, 4, 4, 3, 3,
                                1, 1, 5, 5, 4, 4,
                                2, 2, 3, 3, 5, 5};

int material2node[MaxPieces] = {0,1,1,2,2,0,
                                0,1,1,3,3,0,
                                1,2,2,3,3,1,
                                2,0,0,3,3,2};

int face2material[NumSplitTurnOns][2] = {{5,0},{0,1},{1,2},{2,3},{3,4},{4,5},
                                         {11,6},{6,7},{7,8},{8,9},{9,10},{10,11},
                                         {17,12},{12,13},{13,14},{14,15},{15,16},{16,17},
                                         {23,18},{18,19},{19,20},{20,21},{21,22},{22,23},
                                         {6,0},{7,1},{12,2},{13,3},{18,4},{19,5},
                                         {11,20},{10,21},{8,17},{9,16},{14,23},{15,22},
};

int face2edge_center[NumSplitTurnOns] = {-1,0,-1,1,-1,2,
                                         -1,0,-1,4,-1,3,
                                         -1,1,-1,5,-1,4,
                                         -1,2,-1,3,-1,5,
                                         0,0,1,1,2,2,3,3,4,4,5,5};

int face2face_center[NumSplitTurnOns] = {0,0,0,0,0,0,
                                         1,1,1,1,1,1,
                                         2,2,2,2,2,2,
                                         3,3,3,3,3,3,
                                         -1,-1,-1,-1,-1,-1,
                                         -1,-1,-1,-1,-1,-1};


int face2nodes[NumFacesPerTet][NumNodesPerTriangle] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
int face2opposite_node[NumFacesPerTet] = {3, 2, 0, 1};

int node2face[NumNodesPerTet][3] = {{0,1,3},{0,1,2},{0,3,2},{1,2,3}};
int node2edge[NumNodesPerTet][3][2] = {{{0,2},{0,3},{2,3}},{{0,1},{0,4},{1,4}},{{1,2},{2,5},{1,5}},{{3,4},{4,5},{3,5}}};

int edge2turnon[NumEdgesPerTet][2] = {{1,7},{3,13},{5,19},{11,21},{9,17},{15,23}};
int node2turnon[NumNodesPerTet][3] = {{0,6,20},{2,8,12},{4,14,18},{10,16,22}};

double coefmin = 1;
double coefmax = 0;
template<class T> 
MESH_CUTTING<T>::TET_CUTTING::TET_CUTTING()
{
    material_ids = IDENTITY_ARRAY<>(MaxPieces);
    for (int l = 0; l < material_ids.m; l++) { has_material(material_ids(l)) = 1; }
    tet_center.value = CENTER::Constant_Vector(1.0/NumNodesPerTet);
    for (int i = 0; i < NumFacesPerTet; i++){
        for (int j = 0; j < NumNodesPerTriangle; j++){
            face_centers(i).value(face2node[i][j]) = 1.0/3;
    }}
    for (int i = 0; i < NumEdgesPerTet; i++){
        edge_centers(i).value(edge2node[i][0]) = 0.5;
        edge_centers(i).value(edge2node[i][1]) = 0.5;
    }
}

template<class T> 
typename MESH_CUTTING<T>::TET_CUTTING MESH_CUTTING<T>::TET_CUTTING::Generate_Sub_Tet(const ARRAY<int>& material_ids_input)
{
    //set centers
    for (int j = 0; j < NumSplitTurnOns; j++){
        //has split along this face
        if(material_ids_input.Contains(face2material[j][0]) ^ material_ids_input.Contains(face2material[j][1])){
            if(tet_center.Is_Not_Set()) tet_center.Set();
            //these two lines might return assert error if the split is along a full edge, where the center is never appended
            if(face_centers[face2face_center[j]].Is_Not_Set()) face_centers[face2face_center[j]].Set();
            if(edge_centers[face2edge_center[j]].Is_Not_Set()) edge_centers[face2edge_center[j]].Set();
    }}
    TET_CUTTING tc = *this;
    tc.material_ids = material_ids_input;
    tc.has_material = VECTOR<bool, MaxPieces>();
    for (int l = 0; l < material_ids_input.m; l++) {
         tc.has_material(material_ids_input(l)) = 1;
    }
    return tc;
}

template<class T> 
ARRAY<int> MESH_CUTTING<T>::TET_CUTTING::Find_CC(const int& material_id, VECTOR<bool, MaxPieces>& picked)
{
    ARRAY<int> cc;
    cc.Append(material_id);
    picked[material_id] = 1;
    for (int i = 0; i < cc.m; i++){
        int r = cc(i);
        for (int j = 0; j < NumFacesPerTet-1; j++){
            int face_id = material2face[r][j];
            if (!(turned_on[face_id])){
                for (int k = 0; k < 2; k++){
                    int neighbor = face2material[face_id][k];
                    if (neighbor!=r){
                        if (has_material[neighbor] && !(picked[neighbor])){
                            cc.Append(neighbor);
                            picked[neighbor] = 1;
    }}}}}}
    return cc;
}

template<class T> 
MESH_CUTTING<T>::MESH_CUTTING(){}

template<class T> 
MESH_CUTTING<T>::MESH_CUTTING(TETRAHEDRALIZED_VOLUME<T>* sim_volume_input, T timestep_input, int ratio_input): sim_volume(sim_volume_input), timestep(timestep_input), ratio(ratio_input)
{
    Initialize_Cutting_Volume();
    Update_For_Draw();
}

template<class T>
T MESH_CUTTING<T>::levelset_eyeball(const TV& p)
{
    return (p-TV(0,0,-0.6)).Magnitude_Squared() - 0.09;
}

template<class T>
void MESH_CUTTING<T>::Subdivide_Cutting_Mesh_Into_Eyeball()
{
    //initialize intersections and edge2tet mapping
    HASHTABLE<EDGE, T> intersections;
    HASHTABLE<EDGE, ARRAY<int> > edge2tet;
    T b = -1;
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        for (int j = 0; j < NumEdgesPerTet; j++) {
            int i1 = volume->mesh.elements(i)(edge2node[j][0]);
            int i2 = volume->mesh.elements(i)(edge2node[j][1]);
            EDGE e(i1, i2);
            e.Sort();
            ARRAY<int>& t = edge2tet.Get_Or_Insert(e);
            t.Append(i);
            if (!intersections.Get(e, b)) {
                TV p1 = volume->particles.X(e[0]);
                TV p2 = volume->particles.X(e[1]);
                T l1 = levelset_eyeball(p1);
                T l2 = levelset_eyeball(p2);
                if ((l1 * l2) < 0) {
                    intersections.Set(e, fabs(l2)/(fabs(l1)+fabs(l2)));
                }
            }
        }
    }
    
    //subdivide
    for (HASHTABLE_ITERATOR<EDGE, T> it(intersections); it.Valid(); it.Next()) {
        T w = it.Data();
        if (w > 0) {
            EDGE e = it.Key();
            TV p1 = volume->particles.X(e[0]);
            TV p2 = volume->particles.X(e[1]);
            ARRAY<int> tets;
            edge2tet.Get(e, tets);
            int new_node_id = weights_in_sim.m;
            for (int i = 0; i < tets.m; i++) {
                //subdivide a tet into 2 along the edge intersection
                int old_element_id = tets(i);
                TET t1 = volume->mesh.elements(old_element_id);
                TET t2 = t1;
                
                //new node
                if (i == 0) {
                    CENTER c;
                    for (int i = 0; i < NumNodesPerTet; i++) {
                        if (t1[i] == e[0]) {
                            c[i] = w;
                        }
                        else if (t1[i] == e[1]) {
                            c[i] = 1 - w;
                        }
                    }
                    weights_in_sim.Append(PARENT(ctet2stet(old_element_id), c));
                }
                
                //new elements
                int ii[2];
                for (int j = 0, k = 0; j < NumNodesPerTet; j++) {
                    if (t1[j]!=e[0] && t1[j]!=e[1]) {
                        ii[k] = t1[j];
                        ++k;
                    }
                }
                for (int j = 0; j < NumNodesPerTet; j++) {
                    if (t1[j] == e[1]) {
                        t1[j] = new_node_id;
                        break;
                    }
                }
                volume->mesh.elements(tets(i)) = t1;
                for (int j = 0; j < NumNodesPerTet; j++) {
                    if (t2[j] == e[0]) {
                        t2[j] = new_node_id;
                        break;
                    }
                }
                volume->mesh.elements.Append(t2);
                ctet2stet.Append(ctet2stet(tets(i)));
                int new_element_id = volume->mesh.elements.m-1;
                
                //update edge2tet
                if (ARRAY<int>* ts = edge2tet.Get_Pointer(EDGE(ii[0], ii[1]).Sorted())) {
                    ts->Append(new_element_id);
                }
                for (int k = 0; k < 2; k++) {
                    if (ARRAY<int>* ts = edge2tet.Get_Pointer(EDGE(e[1], ii[k]).Sorted())){
                        int id = -1;
                        if (!ts->Find(old_element_id, id)){
                            cout << "should be able to find!!!\n";
                        }
                        else {
                            ts->operator()(id) = new_element_id;
                        }
                    }
                }
            }
        }
    }
    tet_cuttings.Resize(volume->mesh.elements.m);
    Update_Cutting_Particles();
    
    //eyeball color information 
    is_blue.Resize(volume->mesh.elements.m);
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        is_blue(i) = 1;
        for (int j = 0; j < NumNodesPerTet; j++) {
            if (levelset_eyeball(volume->particles.X(volume->mesh.elements(i)(j))) > 0) {
                is_blue(i) = 0;
                break;
            }
        }
    }
}

template<class T>
void MESH_CUTTING<T>::Set_Dirichlet_Nodes_For_Eyeball()
{
    for (int i = 0; i < sim_volume->particles.X.m; i++) {
        if (sim_volume->particles.X(i)(2) > 0 && sim_volume->particles.X(i).Magnitude_Squared() < 0.25) {
            diri_nodes.Set(i);
        }
    }
}

template<class T>
void MESH_CUTTING<T>::Initialize_Cutting_Volume()
{    
    //create cutting volume, initially copied from sim_volume
    volume = TETRAHEDRALIZED_VOLUME<T>::Create();    
    int initial_num_nodes = sim_volume->particles.X.m;
    weights_in_sim.Resize(initial_num_nodes);
    ARRAY<bool> node_visited(initial_num_nodes);
    for (int i = 0; i < initial_num_nodes; i++){
        node_visited(i) = 0;
    }
    
    for (int i = 0; i < sim_volume->mesh.elements.m; i++){
        TET e = sim_volume->mesh.elements(i);
        volume->mesh.elements.Append(e);
        ctet2stet.Append(i);
        for (int j = 0; j < NumNodesPerTet; j++){
            CENTER c; c[j] = 1;
            if(!node_visited(e[j])){
               weights_in_sim(e[j]).id = i;
               weights_in_sim(e[j]).weight = c;
               node_visited(e[j]) = 1;
            }
        }
    }
    tet_cuttings.Resize(sim_volume->mesh.elements.m);
    Update_Cutting_Particles();
    
    //refine cutting volume
    int num_refinements = 0;
    for (int i = 0; i < num_refinements; i++){
        Refine_Cutting_Volume();
    }
    
    cutting_particle_material_space.Resize(volume->particles.X.m);
    for(int i=0;i<volume->particles.X.m;++i)
        cutting_particle_material_space(i)=volume->particles.X(i);
        
    //subdevide cutting volume so it can be an eyeball...
    is_blue.Resize(volume->mesh.elements.m);
//    Subdivide_Cutting_Mesh_Into_Eyeball();
//    Set_Dirichlet_Nodes_For_Eyeball();
    
}

template<class T> 
typename MESH_CUTTING<T>::CENTER MESH_CUTTING<T>::Weight_In_Sim_Tet(const CENTER& weight_in_element, const TET& element, const int sim_tet_id)
{
    CENTER c;
    for (int k = 0; k < NumNodesPerTet; k++){
        CENTER cc;
        int parent_id = weights_in_sim(element[k]).id;
        if (parent_id == sim_tet_id){
            cc = weights_in_sim(element[k]).weight;
        }
        else {
            for (int i = 0; i < NumNodesPerTet; i++) {
                for (int j = 0; j < NumNodesPerTet; j++) {
                    if (sim_volume->mesh.elements(sim_tet_id)(i) == sim_volume->mesh.elements(parent_id)(j)) {
                        cc[i] = weights_in_sim(element[k]).weight[j];
                        break;
                    }
                }
            }
        }
        c += (weight_in_element[k] * cc);
    }
    return c;
}

template<class T> 
typename MESH_CUTTING<T>::CENTER MESH_CUTTING<T>::Weight_In_Sim_Tet(const CENTER& weight_in_element, const TET& element, const int sim_tet_id, const ARRAY<TET>& original_mesh)
{
    CENTER c;
    for (int k = 0; k < NumNodesPerTet; k++){
        CENTER cc;
        int parent_id = weights_in_sim(element[k]).id;
        if (parent_id == sim_tet_id){
            cc = weights_in_sim(element[k]).weight;
        }
        else {
            for (int i = 0; i < NumNodesPerTet; i++) {
                for (int j = 0; j < NumNodesPerTet; j++) {
                    if (original_mesh(sim_tet_id)(i) == original_mesh(parent_id)(j)) {
                        cc[i] = weights_in_sim(element[k]).weight[j];
                        break;
                    }
                }
            }
        }
        c += (weight_in_element[k] * cc);
    }
    return c;
}

template<class T> 
void MESH_CUTTING<T>::Refine_Cutting_Volume()
{
    HASHTABLE<VECTOR<int,2>,int> new_nodes;
    HASHTABLE<VECTOR<int,3>,int> new_nodes2;
    ARRAY<TET> new_elements;
    ARRAY<TET_CUTTING> new_cuttings;
    ARRAY<int> new_ctet2stet;
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        const int parent_sim_tet = ctet2stet(i);
        VECTOR<int,4> &element = volume->mesh.elements(i);
        TET_CUTTING& tc = tet_cuttings(i);
        
        //tet center
        int cid = weights_in_sim.m;
        {
            CENTER c = Weight_In_Sim_Tet(tc.tet_center.value, element, parent_sim_tet);
            weights_in_sim.Append(PARENT(parent_sim_tet, c));
        }
        
        //edge centers
        VECTOR<int, NumEdgesPerTet> ec;
        for (int j = 0; j < NumEdgesPerTet; j++){
            CENTER c = Weight_In_Sim_Tet(tc.edge_centers[j].value, element, parent_sim_tet);
            VECTOR<int,2> edge = VECTOR<int,2>(element(edge2node[j][0]),element(edge2node[j][1])).Sorted();
            //store the weight
            if (!(new_nodes.Get(edge,ec[j]))){
                int eid = weights_in_sim.m;
                ec[j] = eid;
                new_nodes.Set(edge,eid);
                weights_in_sim.Append(PARENT(parent_sim_tet, c));
            }
        }
        
        VECTOR<int, NumFacesPerTet> fc;
        for (int j = 0; j < NumFacesPerTet; j++){
            CENTER c = Weight_In_Sim_Tet(tc.face_centers[j].value, element, parent_sim_tet);
            TV_INT face = TV_INT(element(face2node[j][0]),element(face2node[j][1]),element(face2node[j][2])).Sorted();
            if (!(new_nodes2.Get(face,fc[j]))){
                int fid = weights_in_sim.m;
                fc[j] = fid;
                new_nodes2.Set(face,fid);
                weights_in_sim.Append(PARENT(parent_sim_tet, c));
            }
        }
        
        VECTOR<bool,4> has_node;
        for (int j = 0; j < tc.material_ids.m; j++) {
            has_node(material2node[tc.material_ids(j)]) = 1;
        }
        for (int j = 0; j < 4; j++){
            if(has_node(j)) {
                for (int k = 0; k < 3; k++){
                    new_elements.Append(VECTOR<int,4>(element(j), cid, fc[node2face[j][k]], ec[node2edge[j][k][0]]));
                    new_elements.Append(VECTOR<int,4>(element(j), cid, fc[node2face[j][k]], ec[node2edge[j][k][1]]));
                    new_cuttings.Append(TET_CUTTING());
                    new_cuttings.Append(TET_CUTTING());
                    new_ctet2stet.Append(parent_sim_tet);
                    new_ctet2stet.Append(parent_sim_tet);
                }
            }
        }    
    }
    
    tet_cuttings = new_cuttings;
    volume->mesh.elements = new_elements;
    ctet2stet = new_ctet2stet;
    Update_Cutting_Particles();
}

template<class T>
void MESH_CUTTING<T>::Update_Cutting_Particles()
{
    start_timer();
    volume->particles.Resize(weights_in_sim.m);
    for (int i = 0; i < weights_in_sim.m; i++){
        volume->particles.X(i) = weight2vec_sim(weights_in_sim(i).id, weights_in_sim(i).weight);
    }
    stop_timer();
    //printf("update_cutting_particles:    %f\n",get_time());
}

template<class T> 
MESH_CUTTING<T>::~MESH_CUTTING()
{


}


template<class T>
void MESH_CUTTING<T>::Initialize_Elasticity()
{
    start_timer();
    int SIZE_MESH = sim_volume->mesh.elements.m;
    int NUM_NODES = sim_volume->particles.X.m;

    deformable_object = new DEFORMABLE_OBJECT_3D<T>(SIZE_MESH, NUM_NODES, 0);
    our_mesh = new GEOMETRY::TETRAHEDRON_MESH(SIZE_MESH,NUM_NODES);

    for (int i = 0; i < NUM_NODES; i++) {
        for (int k = 0; k < 3; k++){
            deformable_object->Positions()(i*3+k) = sim_volume->particles.X(i)(k);
            deformable_object->Velocities()(i*3+k) = sim_volume->particles.V(i)(k);
        }
    }
        
    for (int i = 0; i < SIZE_MESH; i++){
        for (int k = 0; k<4; k++){
            our_mesh->Nodes_Of_Element(i)(k) = sim_volume->mesh.elements(i)(k);
        }
    }
    deformable_object->Tetrahedron_Mesh() = *our_mesh;

    deformable_object->Tetrahedron_Mesh().Initialize_Oriented_Boundary_Triangles();
    //set up physics
    //elastic constitutive model
    fem = new FEM_HYPERELASTICITY_3D<T>(deformable_object->Tetrahedron_Mesh(),deformable_object->Positions());
    le = new FIXED_COROTATED_ELASTICITY_3D<T>((T)100000,(T).3,fem->F());
    fem->Set_Constitutive_Model(*le);
    fem->Initialize_Undeformed_Configuration();
    
    int undeformed_size = fem->Dm_inverse.Size();
    undeformed_config_copy.Resize(undeformed_size);
    for (int i = 0; i < undeformed_size; i++){
        undeformed_config_copy(i) = fem->Dm_inverse(i);
    }
    
    //Choose time stepping: backward Euler
    T dt = (timestep/ratio)/1000;
    int number_steps = 1000;
    T start_time = 0;
    T end_time = number_steps*dt;
    be = new BACKWARD_EULER_TIME_STEPPING_3D<T>(dt,end_time,start_time,*deformable_object);
    be->Set_Elastic_Forces(*fem);
    nodal_volumes = new ALGEBRA::VECTOR<T>(deformable_object->Positions().Size());
    fem->Nodal_Volume_Fractions(*nodal_volumes);
    
    be->Initialize_BE_Matrix(*nodal_volumes);
    //be->Initialize_CG();
    be->Initialize_MINRES();
    
    my_constrained = new ALGEBRA::VECTOR<int>(3*diri_nodes.Size());
    my_constrained_locations = new ALGEBRA::VECTOR<T>(3*diri_nodes.Size());
    int i = 0;
    for (HASHTABLE_ITERATOR<int> it(diri_nodes); it.Valid(); it.Next()) {
        int fixed_node = it.Key();
        my_constrained->Set_Value(3*i,3*fixed_node); my_constrained->Set_Value(3*i+1,3*fixed_node+1); my_constrained->Set_Value(3*i+2,3*fixed_node+2); // first node constrained
        my_constrained_locations->Set_Value(3*i,sim_volume->particles.X(fixed_node)(0));
        my_constrained_locations->Set_Value(3*i+1,sim_volume->particles.X(fixed_node)(1));
        my_constrained_locations->Set_Value(3*i+2,sim_volume->particles.X(fixed_node)(2));
        i++;
    }
    
    //z component of every node fixed, x,y component of all diri nodes fixed.
//    my_constrained = new ALGEBRA::VECTOR<int>(2*diri_nodes.Size()+sim_volume->particles.X.m);
//    my_constrained_locations = new ALGEBRA::VECTOR<T>(2*diri_nodes.Size()+sim_volume->particles.X.m);
//    int i = 0;
//    for (HASHTABLE_ITERATOR<int> it(diri_nodes); it.Valid(); it.Next()) {
//        int fixed_node = it.Key();
//        my_constrained->Set_Value(i,3*fixed_node);
//        my_constrained->Set_Value(i+1,3*fixed_node+1);
//        my_constrained_locations->Set_Value(i,sim_volume->particles.X(fixed_node)(0));
//        my_constrained_locations->Set_Value(i+1,sim_volume->particles.X(fixed_node)(1));
//        i += 2;
//    }
//    for (int j = 0; j < sim_volume->particles.X.m; ++j) {
//        my_constrained->Set_Value(i,3*j+2);
//        my_constrained_locations->Set_Value(i,sim_volume->particles.X(j)(2));
//        i++;
//    }
    
    be->Set_Boundary_Conditions(*my_constrained, *my_constrained_locations);
    stop_timer();
    printf("elasticity initialization time:    %f\n",get_time());
}

template<class T>
void MESH_CUTTING<T>::Reinitialize_Elasticity()
{
    delete deformable_object;
    delete our_mesh;
    delete be;
    delete le; 
    delete fem;
    delete nodal_volumes;
    
    int SIZE_MESH = sim_volume->mesh.elements.m;
    int NUM_NODES = sim_volume->particles.X.m;

    deformable_object = new DEFORMABLE_OBJECT_3D<T>(SIZE_MESH, NUM_NODES, 0);
    our_mesh = new GEOMETRY::TETRAHEDRON_MESH(SIZE_MESH,NUM_NODES);

    for (int i = 0; i < NUM_NODES; i++) {
        for (int k = 0; k < 3; k++){
            deformable_object->Positions()(i*3+k) = sim_volume->particles.X(i)(k);
            deformable_object->Velocities()(i*3+k) = sim_volume->particles.V(i)(k);
        }
    }
        
    for (int i = 0; i < SIZE_MESH; i++){
        for (int k = 0; k<4; k++){
            our_mesh->Nodes_Of_Element(i)(k) = sim_volume->mesh.elements(i)(k);
        }
    }
    deformable_object->Tetrahedron_Mesh() = *our_mesh;

    deformable_object->Tetrahedron_Mesh().Initialize_Oriented_Boundary_Triangles();
    //set up physics
    //elastic constitutive model
    fem = new FEM_HYPERELASTICITY_3D<T>(deformable_object->Tetrahedron_Mesh(),deformable_object->Positions());
    le = new FIXED_COROTATED_ELASTICITY_3D<T>((T)100000,(T).3,fem->F());
    fem->Set_Constitutive_Model(*le);
    for (int i = 0; i < undeformed_config_copy.m; i++){
         fem->Dm_inverse(i) = undeformed_config_copy(i);
         //cout << fem->Dm_inverse(i).Determinant() << endl;
    }
    
    //Choose time stepping: backward Euler
    T dt = (timestep/ratio)/1000;
    int number_steps = 1000;
    T start_time = 0;
    T end_time = number_steps*dt;
    be = new BACKWARD_EULER_TIME_STEPPING_3D<T>(dt,end_time,start_time,*deformable_object);
    be->Set_Elastic_Forces(*fem);
    nodal_volumes = new ALGEBRA::VECTOR<T>(deformable_object->Positions().Size());
    fem->Nodal_Volume_Fractions(*nodal_volumes);
    be->Initialize_BE_Matrix(*nodal_volumes);
    //be->Initialize_CG();
    be->Initialize_MINRES();
    
    delete my_constrained;
    delete my_constrained_locations;
//    my_constrained = new ALGEBRA::VECTOR<int>(3*diri_nodes.Size());
//    my_constrained_locations = new ALGEBRA::VECTOR<T>(3*diri_nodes.Size());
//    int i = 0;
//    for (HASHTABLE_ITERATOR<int> it(diri_nodes); it.Valid(); it.Next()) {
//        int fixed_node = it.Key();
//        my_constrained->Set_Value(3*i,3*fixed_node); my_constrained->Set_Value(3*i+1,3*fixed_node+1); my_constrained->Set_Value(3*i+2,3*fixed_node+2); // first node constrained
//        my_constrained_locations->Set_Value(3*i,sim_volume->particles.X(fixed_node)(0));
//        my_constrained_locations->Set_Value(3*i+1,sim_volume->particles.X(fixed_node)(1));
//        my_constrained_locations->Set_Value(3*i+2,sim_volume->particles.X(fixed_node)(2));
//        i++;
//    }
    
    
    //z component of every node fixed, x,y component of all diri nodes fixed.
    my_constrained = new ALGEBRA::VECTOR<int>(2*diri_nodes.Size()+sim_volume->particles.X.m);
    my_constrained_locations = new ALGEBRA::VECTOR<T>(2*diri_nodes.Size()+sim_volume->particles.X.m);
    int i = 0;
    for (HASHTABLE_ITERATOR<int> it(diri_nodes); it.Valid(); it.Next()) {
        int fixed_node = it.Key();
        my_constrained->Set_Value(i,3*fixed_node);
        my_constrained->Set_Value(i+1,3*fixed_node+1);
        my_constrained_locations->Set_Value(i,sim_volume->particles.X(fixed_node)(0));
        my_constrained_locations->Set_Value(i+1,sim_volume->particles.X(fixed_node)(1));
        i += 2;
    }
    for (int j = 0; j < sim_volume->particles.X.m; ++j) {
        my_constrained->Set_Value(i,3*j+2);
        my_constrained_locations->Set_Value(i,sim_volume->particles.X(j)(2));
        i++;
    }
    
    
    be->Set_Boundary_Conditions(*my_constrained, *my_constrained_locations);
    
    
}

template<class T>
void MESH_CUTTING<T>::Split(const int& tet_id, const ARRAY<int>& tri_ids, const TRIANGULATED_SURFACE<T>& cutting_surface, ARRAY<int>& sim_node_from, ARRAY<bool>& sim_tet_split)
{
    T eps = 1e-6;
    TET_CUTTING tc = tet_cuttings(tet_id);
    VECTOR<bool, NumSplitTurnOns> can_split;
    for (int i = 0; i < NumSplitTurnOns; i++){
        if(tc.has_material(face2material[i][0]) && tc.has_material(face2material[i][1])) {
            can_split(i) = 1;
    }}

    const VECTOR<int, NumNodesPerTet> tet_element = volume->mesh.elements(tet_id);
    int NumTris = tri_ids.m;
    bool overlap = 0;
    for (int i = 0; i < NumTris; i++) {
        //cout << endl;
        VECTOR<int, NumNodesPerTriangle> tri_element = cutting_surface.mesh.elements(tri_ids(i));
        TV& x1 = cutting_surface.particles.X(tri_element(0));
        TV& x2 = cutting_surface.particles.X(tri_element(1));
        TV& x3 = cutting_surface.particles.X(tri_element(2));
        TRIANGLE_3D<T> triangle(x1, x2, x3);
        //if triangle overlap with any tet face, mark it as non-mergable
        bool overlap1 = 0;
        for (int j = 0; j < NumFacesPerTet; j++){
            TV p1, p2, p3;
            int node_num1 = face2node[j][0];
            int node_num2 = face2node[j][1];
            int node_num3 = face2node[j][2];
            p1 = volume->particles.X(tet_element(node_num1));
            p2 = volume->particles.X(tet_element(node_num2));
            p3 = volume->particles.X(tet_element(node_num3));
            MATRIX<T,3> m1(p1-p3, p2-p3, x2-x1);
            MATRIX<T,3> m2(p1-p3, p2-p3, x3-x1);
            
            if (!m1.Determinant() && !m2.Determinant() && (triangle.Inside(p1, 1e-5) || triangle.Inside(p2, 1e-5) || triangle.Inside(p3, 1e-5))) {
                cutFaces.Set(VECTOR<int,3>(tet_element(node_num1), tet_element(node_num2), tet_element(node_num3)).Sorted());
                overlap = 1;
                break;
            }
        }
        if (overlap1) {
            overlap = 1;
            continue;
        }
        
        //tet centers
        ARRAY<CENTER> nodes_inside_tet;
        for (int j = 0; j < NumNodesPerTriangle; j++){
            MATRIX<T,3> ls(volume->particles.X(tet_element(0))-volume->particles.X(tet_element(3)), volume->particles.X(tet_element(1))-volume->particles.X(tet_element(3)), volume->particles.X(tet_element(2))-volume->particles.X(tet_element(3)));
            if (ls.Determinant() != 0) {
                TV w = ls.Inverse() * (cutting_surface.particles.X(tri_element(j))-volume->particles.X(tet_element(3)));
                if (w[0] > 0 && w[1] > 0 && w[2] > 0 
                    && (w[0]+w[1]+w[2]) < 1){
                    CENTER c(w[0], w[1], w[2], 1-w[0]-w[1]-w[2]);
                    nodes_inside_tet.Append(c);            
                }
            }
//            else {
//                cout << "bad element!" << endl;
//                exit(1);
//            }
        }
        if (nodes_inside_tet.m == NumNodesPerTriangle){ continue; }
        
        //otherwise, might have something to turn on unless it's outside the tet. To decide, the following intersection information is needed. We turn on turn off face material when getting the information. After that, we turn on segments on faces and edges.
        ARRAY<CENTER> intersections;//for computing potential center if no node of the cutting triangle is inside the tet
        ARRAY<CENTER> vcenters[NumNodesPerTet];
        ARRAY<CENTER> ecenters[NumEdgesPerTet];
        ARRAY<CENTER> fcenters[NumFacesPerTet];

        //for each face(excluding boundary) of the tet, compute its intersection with the cutting triangle's edges
        for (int j = 0; j < NumFacesPerTet; j++){
            TV p1, p2, p3;
            int node_num1 = face2node[j][0];
            int node_num2 = face2node[j][1];
            int node_num3 = face2node[j][2];
            p1 = volume->particles.X(tet_element(node_num1));
            p2 = volume->particles.X(tet_element(node_num2));
            p3 = volume->particles.X(tet_element(node_num3));
            //for each edge of the cutting triangle
            for (int k=0; k < NumEdgesPerTriangle; k++){
                TV& q1 = cutting_surface.particles.X(tri_element(k));
                TV& q2 = cutting_surface.particles.X(tri_element((k+1)%3));
                MATRIX<T,3> ls(p1-p3, p2-p3, q2-q1);
                if (ls.Determinant() != 0) {
                    TV coef = ls.Inverse() * (q2-p3);
                    //cout << coef << endl;
                    if (coef[0] > 0 && coef[1] > 0
                        && coef[0] + coef[1] < 1
                        && coef[2] >= 0 && coef[2] <= 1) {
                        //intersection
                        CENTER c(0,0,0,0);
                        c[node_num1] = coef[0];
                        c[node_num2] = coef[1];
                        c[node_num3] = 1 - coef[0] - coef[1];
                        fcenters[j].Append(c);
                        intersections.Append(c);
                    }
                }
                else {
                    //cout << "cutting triangle's edge parallel to the tet face\n";
                    //exit(1);
                }
            }
        }
  
        //for each edge of the tet
        //cutting triangle intersecting the edge?
        for (int j = 0; j < NumEdgesPerTet; j++){
            int node_num1 = edge2node[j][0];
            TV y1 = volume->particles.X(tet_element(node_num1));
            
            int node_num2 = edge2node[j][1];
            TV y2 = volume->particles.X(tet_element(node_num2));
        
                        
            MATRIX<T,3> ls(x1-x3, x2-x3, y2-y1);
            if (ls.Determinant() != 0) {
                TV coef = ls.Inverse() * (y2-x3);
        
                if (coef[0] >= 0 && coef[1] >= 0
                    && coef[0] + coef[1] <= 1){
                    if (coef[2] > 0 && coef[2] < 1) {
                        if (coef[2] > coefmax) coefmax = coef[2]; 
                        if (coef[2] < coefmin) coefmin = coef[2]; 
                        //intersection on the edge
                        CENTER c(0, 0, 0, 0);
                        c[node_num1] = coef[2];
                        c[node_num2] = 1 - coef[2];
                        ecenters[j].Append(c);
                        intersections.Append(c);
                        //cout << "intersect with edge " << j << endl;
                    }
                    else if (fabs(coef[2]-1)<eps) {
                        cout << "cut on a node " << node_num1 << endl;
//                        x1[0] += 0.001;
//                        x2[0] += 0.001;
//                        x3[0] += 0.001;
                        //return(-1);
                        //exit (1);q
                        if (vcenters[node_num1].m == 0){
                            CENTER c(0, 0, 0, 0);
                            c[node_num1] = 1;
                            vcenters[node_num1].Append(c);
                            intersections.Append(c);
                        }    
                    }
                    else if (fabs(coef[2])<eps) {
                        cout << "cut on a node " << node_num2 << endl;
                        //                        x1[0] += 0.001;
                        //                        x2[0] += 0.001;
                        //                        x3[0] += 0.001;
                        //return(-1);
                        //exit (1);q
                        if (vcenters[node_num2].m == 0){
                            CENTER c(0, 0, 0, 0);
                            c[node_num2] = 1;
                            vcenters[node_num2].Append(c);
                            intersections.Append(c);
                        }    
                    }
                }
            }
//            else {//seems never fail to split, always extra splits due to inconsistency
//                cout << "parallel************************************************************\n";
//            }    
        }
        //cutting triangle cut along an edge?
        for (int j = 0; j < NumEdgesPerTet; j++){
            int node_num1 = edge2node[j][0];
            TV y1 = volume->particles.X(tet_element(node_num1));
            
            int node_num2 = edge2node[j][1];
            TV y2 = volume->particles.X(tet_element(node_num2));
            
            if (MATRIX<T,3>(x1-x3, x2-x3, y2-y1).Determinant() == 0 && (vcenters[node_num1].m || vcenters[node_num2].m) && (ecenters[0].m || ecenters[1].m || ecenters[2].m || ecenters[3].m || ecenters[4].m || ecenters[5].m || fcenters[0].m || fcenters[1].m || fcenters[2].m || fcenters[3].m)) {
                cout << "cut along edge " << j << ": ";
                CENTER c;
                c[node_num1] = 0.5;
                c[node_num2] = 0.5;
                intersections.Append(c);
                if (vcenters[node_num1].m) {
                    //cout << node_num1 << "    ";
                    tc.turned_on[24+2*j] = 1;
                }
                if (vcenters[node_num2].m) {
                    //cout << node_num2;
                    tc.turned_on[24+2*j+1] = 1;
                }
                //cout << endl;
                break;
            }
        }
        
        //now I know whethe they intersect or not. if so, add the potential tet center
        bool they_do_intersect = 0;
        if (nodes_inside_tet.m == 1 || nodes_inside_tet.m == 2){
            they_do_intersect = 1;
            if (tc.tet_center.Is_Not_Set()){
                tc.tet_center.Append(nodes_inside_tet.Average());
            }
        }
        else if (intersections.m >= 3) {
            they_do_intersect = 1;
            if (tc.tet_center.Is_Not_Set()){
                tc.tet_center.Append(intersections.Average());
            }
        }

        //turn on segments on the tet faces
        if (they_do_intersect) {
            
            for (int j = 0; j < NumEdgesPerTet; j++){
                if (tc.edge_centers[j].Is_Not_Set() && ecenters[j].m) {
                    tc.edge_centers[j].Append(ecenters[j](0));
                }
            }
            //turn on segments on the faces
            //for each face
            for (int j = 0; j < NumFacesPerTet; j++){
                int e1 = face2edge[j][0];
                int e2 = face2edge[j][1];
                int e3 = face2edge[j][2];
                int node_num1 = face2node[j][0];
                int node_num2 = face2node[j][1];
                int node_num3 = face2node[j][2];
                //face interior and a node or an edge cut
                if (vcenters[node_num1].m+vcenters[node_num2].m+vcenters[node_num3].m+ecenters[e1].m+ecenters[e2].m+ecenters[e3].m==1 && fcenters[j].m == 0) {cout << volume->particles.X(tet_element(face2node[j][0])) << volume->particles.X(tet_element(face2node[j][1])) << volume->particles.X(tet_element(face2node[j][2])) << endl; cout << x1 << x2 << x3 << endl << endl;}
                if (fcenters[j].m == 1) {
                    //decide which one to turn on
                    int turn_on_index = -1;
                    if (vcenters[node_num1].m) {
                        turn_on_index = 6 * j;
                    }
                    else if (vcenters[node_num2].m) {
                        turn_on_index = 6 * j + 2;
                    }
                    else if (vcenters[node_num3].m) {
                        turn_on_index = 6 * j + 4;
                    }
                    else if (ecenters[e1].m) {
                        turn_on_index = 6 * j + 1;
                    }
                    else if (ecenters[e2].m>0) {
                        turn_on_index = 6 * j + 3;
                    }
                    else if (ecenters[e3].m) {
                        turn_on_index = 6 * j + 5;
                    }
                    //turn on
                    if (turn_on_index!=-1){
                        if (can_split[turn_on_index]){
                            tc.turned_on[turn_on_index] = 1;
                            //push back potential centers if not set
                            if (tc.face_centers[j].Is_Not_Set()) {
                                tc.face_centers[j].Append(fcenters[j](0));
                            }
                        }
                    }
                }
                //two edges cut
                else if(ecenters[e1].m > 0 && ecenters[e2].m > 0) {
                    CENTER c = ecenters[e1](0)+ecenters[e2](0);
                    fcenters[j].Append(c*0.5);
                    bool effective_cut = 0;
                    int turn_on_index = 6 * j + 1;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    turn_on_index = 6 * j + 3;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    if (effective_cut && tc.face_centers[j].Is_Not_Set()) {
                        tc.face_centers[j].Append(fcenters[j](0));
                    }
                }
                else if(ecenters[e1].m && ecenters[e3].m) {
                    CENTER c = ecenters[e1](0)+ecenters[e3](0);
                    fcenters[j].Append(c*0.5);
                    bool effective_cut = 0;
                    int turn_on_index = 6 * j + 1;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    turn_on_index = 6 * j + 5;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    if (effective_cut && (tc.face_centers[j].Is_Not_Set())) {
                        tc.face_centers[j].Append(fcenters[j](0));
                    }
                }
                else if(ecenters[e2].m && ecenters[e3].m) {
                    CENTER c = ecenters[e2](0)+ecenters[e3](0);
                    fcenters[j].Append(c*0.5);
                    bool effective_cut = 0;
                    int turn_on_index = 6 * j + 3;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    turn_on_index = 6 * j + 5;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    if (effective_cut && (tc.face_centers[j].Is_Not_Set())) {
                        tc.face_centers[j].Append(fcenters[j](0));
                    }
                }
                //an edge and a node cut
                else if(ecenters[e1].m && vcenters[node_num3].m) {
                    bool effective_cut = 0;
                    CENTER c = ecenters[e1](0)+vcenters[node_num3](0);
                    fcenters[j].Append(c*0.5);
                    int turn_on_index = 6 * j + 1;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    turn_on_index = 6 * j + 4;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    if (effective_cut && tc.face_centers[j].Is_Not_Set()) {
                        tc.face_centers[j].Append(fcenters[j](0));}
                }
                else if(ecenters[e2].m && vcenters[node_num1].m) {
                    bool effective_cut = 0;
                    CENTER c = ecenters[e2](0)+vcenters[node_num1](0);
                    fcenters[j].Append(c*0.5);
                    int turn_on_index = 6 * j + 3;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    turn_on_index = 6 * j;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    if (effective_cut && tc.face_centers[j].Is_Not_Set()) {
                        tc.face_centers[j].Append(fcenters[j](0));}
                }
                else if(ecenters[e3].m && vcenters[node_num2].m) {
                    bool effective_cut = 0;
                    CENTER c = ecenters[e3](0)+vcenters[node_num2](0);
                    fcenters[j].Append(c*0.5);
                    int turn_on_index = 6 * j + 5;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    turn_on_index = 6 * j + 2;
                    if (can_split[turn_on_index]) {
                        tc.turned_on[turn_on_index] = 1;
                        effective_cut = 1;
                    }
                    if (effective_cut && tc.face_centers[j].Is_Not_Set()) {
                        tc.face_centers[j].Append(fcenters[j](0));}
    }}}}
    
    //generate submaterials
    //set centers if it's on a turned-on face that is split(only one piece is in a connected component)
    //cout << tc.turned_on << endl;
    VECTOR<bool, MaxPieces> picked;
    int NumMaterials = tc.material_ids.m;
    for (int i = 0; i < NumMaterials; i++) {
        if (!picked[tc.material_ids(i)]) {
            ARRAY<int> cc = tc.Find_CC(tc.material_ids(i), picked);
            //cout << cc << endl;
            if (cc.m != NumMaterials || overlap) {
                //cout << cc << endl;
                //find parent sim tet_cc and split it
                int parent_sim_tet_id = original_ctet2stet(tet_id);
                TET& parent_sim_tet = original_sim_elements(parent_sim_tet_id);
                
                VECTOR<bool,4> has_mat;
                VECTOR<CENTER,4> wei;//node weights in the parent sim tet
                for(int j = 0; j < cc.m; j++){
                    has_mat(material2node[cc(j)]) = 1;
                }
                //cout << has_mat << endl;
                for (int j = 0; j < NumNodesPerTet; j++){
                    sim_node_from.Append(parent_sim_tet[j]);
                    CENTER c; c[j] = 1;
                    wei[j] = Weight_In_Sim_Tet(c, tet_element, parent_sim_tet_id, original_sim_elements);
                }
                
                //sim mesh related
                int sim_size = sim_node_from.m;
                if (!sim_tet_split(parent_sim_tet_id)){//no sim tet appended
                    //if (i!=0) cout << "ridiculous!*************************************" << endl;
                    for (int j = 0; j < NumNodesPerTet; j++){
                        weights_in_sim.Append(PARENT(parent_sim_tet_id, wei[j]));
                        cutting_particle_material_space.Append(cutting_particle_material_space(tet_element(j)));
                    }
                    sim_volume->mesh.elements(parent_sim_tet_id) = VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1);
                    sim_tet_split(parent_sim_tet_id) = 1;
                }
                else {
                    ALGEBRA::MATRIX_3X3<T> am = undeformed_config_copy(parent_sim_tet_id);
                    undeformed_config_copy.Append(am);
                    for (int j = 0; j < NumNodesPerTet; j++){
                        weights_in_sim.Append(PARENT(sim_volume->mesh.elements.m, wei[j]));
                        cutting_particle_material_space.Append(cutting_particle_material_space(tet_element(j)));
                    }
                    if (i == 0) {
                        ctet2stet(tet_id) = sim_volume->mesh.elements.m;
                    }
                    else {
                        ctet2stet.Append(sim_volume->mesh.elements.m);
                    }
                    sim_volume->mesh.elements.Append(VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1));
                    sim_tet_from.Append(sim_tet_from(parent_sim_tet_id));
                }
                
                //cutting mesh related
                int size = weights_in_sim.m;
                if (i == 0) { 
                    tet_cuttings(tet_id) = tc.Generate_Sub_Tet(cc); 
                    volume->mesh.elements(tet_id) = VECTOR<int, NumNodesPerTet>(size-4, size-3, size-2, size-1);
                }
                else { 
                    tet_cuttings.Append(tc.Generate_Sub_Tet(cc)); 
                    volume->mesh.elements.Append(VECTOR<int, NumNodesPerTet>(size-4, size-3, size-2, size-1));
                    is_blue.Append(is_blue(tet_id));
                }
            }
            else {
                tet_cuttings(tet_id) = tc;
                break;
            }
        }
    }
}

template<class T> 
void MESH_CUTTING<T>::Union(const TV_TV_INT& cface, TV_TV_INT& cface_found, UNION_FIND<int>& node_classes, UNION_FIND<int>& sim_node_classes, const int shift, const bool inverse)
{
    //cout << "UNION\n";
    if (inverse) {
        for (int i = 0; i < NumMaterialsPerFace; i++) {
            int ii = ((9-i)%NumMaterialsPerFace+2*shift)%NumMaterialsPerFace;
            if (cface(ii).face(0) != cface(ii).face(1)) {
                if(cface_found(i).face(0) != cface_found(i).face(1)) {
                    int st1 = cface(ii).sim_tet_id;
                    int st2 = cface_found(i).sim_tet_id;
                    int nst1 = cface(ii).new_sim_tet_id;
                    int nst2 = cface_found(i).new_sim_tet_id;
                    for (int j = 0; j < NumNodesPerTet; j++) {
                        for (int k = 0; k < NumNodesPerTet; k++) {
                            if (original_sim_elements(st1)[j] == original_sim_elements(st2)[k]){
                                //cout << sim_volume->mesh.elements(nst1)[j] << " " << sim_volume->mesh.elements(nst2)[k] << endl;
                                sim_node_classes.Union(sim_volume->mesh.elements(nst1)[j], sim_volume->mesh.elements(nst2)[k]);
                                break;
                            }
                        }
                    }
                    
                    for (int j = 0; j < NumNodesPerTriangle; j++) {
                        node_classes.Union(cface(ii).face((NumNodesPerTriangle-1-j+shift)%NumNodesPerTriangle), cface_found(i).face(j));
                        //cout << "unioned " << cface(ii)((NumNodesPerTriangle-1-j+shift)%NumNodesPerTriangle) << " and " << cface_found(i)(j) << endl;
                    }
                }
                else{
                    for (int j = 0; j < NumNodesPerTriangle; j++){
                        cface_found(i).face(j) = cface(ii).face((NumNodesPerTriangle-1-j+shift)%NumNodesPerTriangle);
                    }
                    cface_found(i).sim_tet_id = cface(ii).sim_tet_id;
                    cface_found(i).new_sim_tet_id = cface(ii).new_sim_tet_id;
    }}}}
    else {
        for (int i = 0; i < NumMaterialsPerFace; i++) {
            int ii = (i+2*shift)%NumMaterialsPerFace;
            if (cface(ii).face(0) != cface(ii).face(1)){
                if(cface_found(i).face(0) != cface_found(i).face(1)) {
                    int st1 = cface(ii).sim_tet_id;
                    int st2 = cface_found(i).sim_tet_id;
                    int nst1 = cface(ii).new_sim_tet_id;
                    int nst2 = cface_found(i).new_sim_tet_id;
                    for (int j = 0; j < NumNodesPerTet; j++) {
                        for (int k = 0; k < NumNodesPerTet; k++) {
                            if (original_sim_elements(st1)[j] == original_sim_elements(st2)[k]){
                                //cout << sim_volume->mesh.elements(nst1)[j] << " " << sim_volume->mesh.elements(nst2)[k] << endl;
                                sim_node_classes.Union(sim_volume->mesh.elements(nst1)[j], sim_volume->mesh.elements(nst2)[k]);
                                break;
                            }
                        }
                    }
                    for (int j = 0; j < NumNodesPerTriangle; j++) {
                        node_classes.Union(cface(ii).face((j+shift)%NumNodesPerTriangle), cface_found(i).face(j));
                        //cout << "unioned " << (j+shift)%NumNodesPerTriangle << " and " << cface_found(i)(j) << endl;
                }}
                else{
                    for (int j = 0; j < NumNodesPerTriangle; j++){
                        cface_found(i).face(j) = cface(ii).face((j+shift)%NumNodesPerTriangle);
                    }
                    cface_found(i).sim_tet_id = cface(ii).sim_tet_id;
                    cface_found(i).new_sim_tet_id = cface(ii).new_sim_tet_id;
    }}}}
}

template<class T> 
void MESH_CUTTING<T>::Cut(TRIANGULATED_SURFACE<T>& cutting_surface, bool refine)
{
    cout << "cutting*************************************" << endl;
    cutFaces.Clean_Memory();
    
    //copy material space positions into meshes
//    for(int i=0;i<volume->particles.X.m;++i)
//        volume->particles.X(i)=cutting_particle_material_space(i);
    
    start_timer();
    //test bounding box hierarchy intersections
//    volume->Initialize_Hierarchy();
    cutting_surface.Initialize_Hierarchy();
//    ARRAY<ARRAY<int> > intersection_list(volume->mesh.elements.m);
//    BOX_VISITOR_TRIVIAL visitor(intersection_list);
//    volume->hierarchy->Intersection_List(*cutting_surface.hierarchy,visitor,0);
//    if(!volume->tetrahedron_list) volume->Update_Tetrahedron_List();
    
    //split
    ARRAY<TET>& elements = volume->mesh.elements;
    ARRAY<TET>& sim_elements = sim_volume->mesh.elements;
    
    //copy of the original sim elements, stay unchanged
    cout << "total particles: " << volume->particles.X.m << endl;
    //cout << elements << endl;
    
    //keep a copy of original element-wise data
    original_elements = volume->mesh.elements;
    original_sim_elements = sim_volume->mesh.elements;
    original_ctet2stet = ctet2stet;
    
    int num_elements = elements.m; cout << num_elements << " elements" << endl;
    int num_nodes = volume->particles.number;
    ARRAY<bool> need_dup(num_nodes);//whether the node needs to be duplicated
    ARRAY<bool> need_merge(num_elements);//whether the element needs to be merged
    ARRAY<int> pieces(num_elements);

    ARRAY<int> sim_node_from(sim_volume->particles.X.m);
    ARRAY<bool> sim_tet_split(sim_volume->mesh.elements.m);
    for (int i = 0; i < sim_volume->particles.X.m; i++){
        sim_node_from(i) = i;
    }
    for (int i = 0; i < sim_volume->mesh.elements.m; i++){
        sim_tet_split(i) = 0;
        sim_tet_from.Append(i);
    }
    
    //sim_node_from, weights_in_sim
    // sim_volume->mesh.elements, always completely newly created, because each child tet(split or not) gets a new parent tet
    //volume->mesh.elements, ctet2stet
    
    int old_mesh_size = num_elements;
    int split_count = 0;
    for (int i = 0; i < num_elements; i++) {
        ARRAY<int> intersection_list;
        cutting_surface.hierarchy->Intersection_List(RANGE<TV>::Bounding_Box(volume->particles.X.Subset(volume->mesh.elements(i))), intersection_list);
//        if (intersection_list(i).m != 0) {
//            Split(i, intersection_list(i), cutting_surface);
        if (intersection_list.m != 0) {
            Split(i, intersection_list, cutting_surface, sim_node_from, sim_tet_split);
            pieces(i) = elements.m - old_mesh_size;
            if (pieces(i) > 0) {// in the box
                //cout << "pices:" << pieces(i) << endl;
                split_count++;
                need_merge(i) = 1;
                need_dup.Subset(original_elements(i)).Fill(1);
                old_mesh_size = elements.m;
            }
        }
    }
    cout << "min: " << coefmin << ", max: " << coefmax << endl;
    cout << split_count << " tets split\n";
    
    //duplicated all touched nodes and split all sim tets that are cut
    for (int i = 0; i < num_elements; i++) {
        TET tet_element = original_elements(i);
        if (!need_merge(i)) {//only nodes on need_merge(i.e. split) elements are duplicated
            VECTOR<bool,4> has_mat;
            int parent_sim_tet_id = original_ctet2stet(i);
            for(int j = 0; j < tet_cuttings(i).material_ids.m; j++){
                has_mat(material2node[tet_cuttings(i).material_ids(j)]) = 1;
            }
            for (int j = 0; j < NumNodesPerTet; j++) {
                if (need_dup(elements(i)(j))) {
                    need_merge(i) = 1;
                    break;
                }
            }    
            if (1){//need_merge(i)){
                //split its parent sim tet
                for (int j = 0; j < NumNodesPerTet; j++){
                    sim_node_from.Append(original_sim_elements(parent_sim_tet_id)(j));
                    elements(i)(j) = weights_in_sim.m;
                    CENTER c; c[j] = 1;
                    CENTER w = Weight_In_Sim_Tet(c, tet_element, parent_sim_tet_id, original_sim_elements);
                    if (!sim_tet_split(parent_sim_tet_id)){   
                        weights_in_sim.Append(PARENT(parent_sim_tet_id, w));
                        cutting_particle_material_space.Append(cutting_particle_material_space(tet_element(j)));
                    }
                    else {
                        weights_in_sim.Append(PARENT(sim_volume->mesh.elements.m, w));
                        cutting_particle_material_space.Append(cutting_particle_material_space(tet_element(j)));
                    }
                }
                int sim_size = sim_node_from.m;
                if (!sim_tet_split(parent_sim_tet_id)){//no sim tet to be appended
                    sim_volume->mesh.elements(parent_sim_tet_id) = VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1);
                    sim_tet_split(parent_sim_tet_id) = 1;
                }
                else {
                    ctet2stet(i) = sim_volume->mesh.elements.m;
                    sim_volume->mesh.elements.Append(VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1));
                    sim_tet_from.Append(i);
                    undeformed_config_copy.Append(undeformed_config_copy(parent_sim_tet_id));
                }
            }
        }
    }
    
//    for (int i = 0; i < num_elements; i++) {
//        if (!need_merge(i)) {
//            int parent_sim_tet_id = original_ctet2stet(i);
//            if (sim_tet_split(parent_sim_tet_id)){
//                for (int j = 0; j < NumNodesPerTet; j++){
//                    sim_node_from.Append(original_sim_elements(parent_sim_tet_id)(j));
//                }
//                ctet2stet(i) = sim_volume->mesh.elements.m;
//                sim_volume->mesh.elements.Append(TET(sim_node_from.m-4, sim_node_from.m-3, sim_node_from.m-2, sim_node_from.m-1));
//                undeformed_config_copy.Append(undeformed_config_copy(parent_sim_tet_id));
//            }
//            need_merge(i) = 1;
//        }
//    }
    
    //cout << sim_volume->mesh.elements.m << " " << volume->mesh.elements.m << " " << ctet2stet.m << " " << undeformed_config_copy.m << endl;
//    for(int i = 0; i < weights_in_sim.m; i++){
//        cout << "(" << weights_in_sim(i).id << ", " << weights_in_sim(i).weight << ") ";
//    } cout << endl;
    
    stop_timer();
    printf("split time:    %f\n",get_time());
    cout << "total particles: " << weights_in_sim.m << endl;
    //cout << elements << endl;
    
    
    //merge nodes
    //identify which un-plit tet need merging
    start_timer();
    
    HASHTABLE<TV_INT, int> face_visited;
    HASHTABLE<TV_INT, TV_TV_INT> oface2cface;//original face mapped to cut face
    UNION_FIND<int> node_classes(weights_in_sim.m);//cluster merging nodes
    UNION_FIND<int> sim_node_classes(sim_node_from.m);
    TET_CUTTING tc;
    VECTOR<int, NumNodesPerTet> element, original_element;
    int piece_id = num_elements;
    int p, np;
    for (int i = 0; i < num_elements; i++) {
        original_element = original_elements(i);
        p = original_ctet2stet(i);
        if (1){//need_merge(i)) {
            for (int k = 0; k < NumFacesPerTet; k++) {
                TV_TV_INT cface, cface1;
                TV_INT original_face(original_element(face2nodes[k][0]), original_element(face2nodes[k][1]), original_element(face2nodes[k][2]));
                if (cutFaces.Contains(original_face.Sorted())) {
                    continue;
                }
                {
                    tc = tet_cuttings(i);
                    element = elements(i);
                    np = ctet2stet(i);
                    for (int l = 0; l < NumMaterialsPerFace; l++) {
                        if (tc.has_material(NumMaterialsPerFace*k+l)) {
                            cface(l) = CUT_FACE(TV_INT(element(face2nodes[k][0]), element(face2nodes[k][1]), element(face2nodes[k][2])), p, np);
                }}}
                for (int j = 0; j < pieces(i); j++) {
                    np = ctet2stet(piece_id + j);
                    tc = tet_cuttings(piece_id + j);
                    element = elements(piece_id + j);
                    for (int l = 0; l < NumMaterialsPerFace; l++) {
                        if (tc.has_material(NumMaterialsPerFace*k+l)) {
                            cface(l) = CUT_FACE(TV_INT(element(face2nodes[k][0]), element(face2nodes[k][1]), element(face2nodes[k][2])), p, np);
                    }}
                }
                bool found = 0;
                for (int l = 0; l < NumNodesPerTriangle; l++) {
                    TV_INT oface(original_face((2+l)%NumNodesPerTriangle), original_face((1+l)%NumNodesPerTriangle), original_face(l));
                    if (oface2cface.Get(oface, cface1)){
                        Union(cface, cface1, node_classes, sim_node_classes, l, 1);
                        oface2cface.Set(oface, cface1);
                        found = 1;
                        break;
                    }
                    oface = TV_INT(original_face(l), original_face((1+l)%NumNodesPerTriangle), original_face((2+l)%NumNodesPerTriangle));
                    if (oface2cface.Get(oface, cface1)){
                        Union(cface, cface1, node_classes, sim_node_classes, l, 0);
                        oface2cface.Set(oface, cface1);
                        found = 1;
                        break;
                }}
                if (!found) { oface2cface.Set(original_face, cface); }
            }
            piece_id+=pieces(i);
        }
    }
    
    
    //reset the particles and mesh indices for sim mesh
    ARRAY<int> sim_element_becomes(sim_elements.m);//element becomes which element in the new mesh
    HASHTABLE<TET,int> sim_tet2id;//whether the tet has been processed or not
    ARRAY<int> sim_new_node_id(sim_node_from.m);//node becomes which node in the new particles
    for (int i = 0; i < sim_new_node_id.m; i++) {
        sim_new_node_id(i) = -1;
    }
    ARRAY<TET> new_sim_elements;
    ARRAY<TV> new_sim_particles;
    ARRAY<int> new_dragging_elements;
    ARRAY<TV> new_dragging_weights;
    ARRAY<TV> new_dragging_targets;
    ARRAY<ALGEBRA::MATRIX_3X3<T> > new_undeformed_config_copy;
    for (int i = 0; i < sim_elements.m; i++){
        TET tet;
        for (int j = 0; j < NumNodesPerTet; j++){
            tet[j] = sim_node_classes.Find(sim_elements(i)(j));
        }
        int tet_id = 0;
        if (!sim_tet2id.Get(tet.Sorted(), tet_id)){//if the tet has not appeared
            sim_tet2id.Set(tet.Sorted(), new_sim_elements.m);
            sim_element_becomes(i) = new_sim_elements.m;
            TET new_tet;
            sim_tet2id.Set(tet,new_sim_elements.m);
            for (int j = 0; j < NumNodesPerTet; j++){
                int new_id = sim_new_node_id(tet[j]);
                if (new_id == -1){
                    sim_new_node_id(tet[j]) = new_sim_particles.m;
                    new_tet[j] = new_sim_particles.m;
                    new_sim_particles.Append(sim_volume->particles.X(sim_node_from(tet[j])));
                }
                else {
                    new_tet[j] = new_id;
                }
            }
            new_sim_elements.Append(new_tet);
            int pst;
            if (dragging_elements.Find(sim_tet_from(i), pst)) {
                new_dragging_elements.Append(new_sim_elements.m-1);
                new_dragging_weights.Append(dragging_weights(pst));
                new_dragging_targets.Append(dragging_targets(pst));
            }
            new_undeformed_config_copy.Append(undeformed_config_copy(i));
        }
        else {
            sim_element_becomes(i) = tet_id;
        }
    }
        
    //reset the particles and mesh indices for cutting mesh
    ARRAY<int> new_node_id(weights_in_sim.m);//node becomes which node in the new particles
    for (int i = 0; i < new_node_id.m; i++) {
        new_node_id(i) = -1;
    }
    ARRAY<PARENT> new_weights_in_sim;
    HASHTABLE<int> new_diri_nodes;
    ARRAY<TV> new_cutting_particles_material_space;
    for (int i = 0; i < elements.m; i++) {
        VECTOR<bool,4> has_mat;
        for(int j = 0; j < tet_cuttings(i).material_ids.m; j++){
            has_mat(material2node[tet_cuttings(i).material_ids(j)]) = 1;
        }
        for (int j = 0; j < NumNodesPerTet; j++) {
            if (has_mat[j]) {
                for (int k = 0; k < NumNodesPerTet; k++) {
                    if (weights_in_sim(elements(i)(j)).weight(k) == 1 && diri_nodes.Contains(sim_node_from(sim_elements(ctet2stet(i))(k)))) {
                        new_diri_nodes.Set(new_sim_elements(sim_element_becomes(ctet2stet(i)))(k));
                    }
                }
            }
            int node_class = node_classes.Find(elements(i)(j));
            int new_id = new_node_id(node_class);
            if (new_id == -1) {//the node hasn't been added
                new_node_id(node_class) = new_weights_in_sim.m;
                elements(i)(j) = new_weights_in_sim.m;
                PARENT p = weights_in_sim(node_class);
                p.id = sim_element_becomes(p.id);
                new_weights_in_sim.Append(p);
                new_cutting_particles_material_space.Append(cutting_particle_material_space(node_class));
            }
            else {
                elements(i)(j) = new_id;
            }
        }
        ctet2stet(i) = sim_element_becomes(ctet2stet(i));
    }
    
    //replace old data with new ones
    diri_nodes = new_diri_nodes;
    cutting_particle_material_space=new_cutting_particles_material_space;
    weights_in_sim = new_weights_in_sim;
    cout<<weights_in_sim.m<<" "<<cutting_particle_material_space.m<<endl;
    dragging_elements = new_dragging_elements;
    dragging_weights = new_dragging_weights;
    dragging_targets = new_dragging_targets;
    sim_volume->mesh.elements = new_sim_elements;
    undeformed_config_copy = new_undeformed_config_copy;
    sim_volume->particles.Resize(new_sim_particles.m);
    for (int i = 0; i < new_sim_particles.m; i++){
        sim_volume->particles.X(i) = new_sim_particles(i);
    }

    //cout << "sim mesh: " << sim_volume->mesh.elements << endl;
    //cout << "cutting mesh: " << volume->mesh.elements << endl;
//    for(int i = 0; i < weights_in_sim.m; i++){
//        cout << "(" << weights_in_sim(i).id << ", " << weights_in_sim(i).weight << ") ";
//    } cout << endl;
    
    stop_timer();
    printf("merge time:    %f\n",get_time());  
    
    original_ctet2stet.Remove_All();
    original_elements.Remove_All();
    original_sim_elements.Remove_All();
    sim_tet_from.Remove_All();
    Update_Cutting_Particles();
    
    //CC of the cutting mesh
    tet_cc.Remove_All();
    node_cc.Remove_All();
    int num_elements1 = volume->mesh.elements.m;
    //find tet connected components
    ARRAY<bool> picked(num_elements1);
    ARRAY<ARRAY<int> > itot(volume->particles.X.m); //node index to tet index
    for (int i = 0; i < num_elements1; i++) {
            VECTOR<int, NumNodesPerTet> element = volume->mesh.elements(i);
            for (int j = 0; j < NumNodesPerTet; j++) {
                itot(element(j)).Append(i); }}
            
    for (int i = 0; i < num_elements1; i++) {
        if (!picked(i)) {
            tet_cc.Append(Find_Cutting_Tet_CC(i, picked, itot, node_cc)); }}
    cout << tet_cc.m << " cutting mesh CCs" << endl;   
    
    if (refine) {
        Partial_Refine();
    }
    //cout << elements << endl;
    //cout << "tet_cc " << tet_cc << endl;
    Update_For_Draw();
    //cout << "tet_cc " << tet_cc << endl;
    
    start_timer();
    Reinitialize_Elasticity();
    stop_timer();
    printf("elasticity reinitialization time:    %f\n",get_time());
}

template<class T>
void MESH_CUTTING<T>::Fix_Orientation()
{
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        if (TETRAHEDRON<T>(volume->particles.X(volume->mesh.elements(i)(0)),
                        volume->particles.X(volume->mesh.elements(i)(1)),
                        volume->particles.X(volume->mesh.elements(i)(2)),
                        volume->particles.X(volume->mesh.elements(i)(3))).Signed_Volume() < 0) {
            int temp = volume->mesh.elements(i)(0);
            volume->mesh.elements(i)(0) = volume->mesh.elements(i)(1);
            volume->mesh.elements(i)(1) = temp;
        }
    }
}

template<class T> 
void MESH_CUTTING<T>::Partial_Refine()
{
    HASHTABLE<VECTOR<int,2>,int> new_nodes;//edge to edge intersection point's index
    HASHTABLE<VECTOR<int,2>, std::set<int> > edge2tets;
    HASHTABLE<VECTOR<int,3>,int> new_nodes2;
    ARRAY<TET> new_elements;
    ARRAY<int> new_ctet2stet;
    ARRAY<bool> new_is_blue;
    cuttingFaces.Clean_Memory();
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        if (tet_cuttings(i).material_ids.m < MaxPieces){
            const int parent_sim_tet = ctet2stet(i);
            VECTOR<int,4> &element = volume->mesh.elements(i);
            TET_CUTTING& tc = tet_cuttings(i);
            
            //tet center
            int cid = weights_in_sim.m;
            {
                CENTER c = Weight_In_Sim_Tet(tc.tet_center.value, element, parent_sim_tet);
                weights_in_sim.Append(PARENT(parent_sim_tet, c));
            }
            
            //edge centers
            VECTOR<int, NumEdgesPerTet> ec;
            for (int j = 0; j < NumEdgesPerTet; j++){
                CENTER c = Weight_In_Sim_Tet(tc.edge_centers[j].value, element, parent_sim_tet);
                //cout << "e weight: " << tc.edge_centers[j] << endl;
                VECTOR<int,2> edge = VECTOR<int,2>(element(edge2node[j][0]),element(edge2node[j][1])).Sorted();
                //store the weight
                if (!(new_nodes.Get(edge,ec[j]))){
                    int eid = weights_in_sim.m;
                    ec[j] = eid;
                    new_nodes.Set(edge,eid);
                    weights_in_sim.Append(PARENT(parent_sim_tet, c));
                }
            }
            
            //face centers
            VECTOR<int, NumFacesPerTet> fc;
            for (int j = 0; j < NumFacesPerTet; j++){
                CENTER c = Weight_In_Sim_Tet(tc.face_centers[j].value, element, parent_sim_tet);
                TV_INT face = TV_INT(element(face2node[j][0]),element(face2node[j][1]),element(face2node[j][2])).Sorted();
                if (!(new_nodes2.Get(face,fc[j]))){
                    int fid = weights_in_sim.m;
                    fc[j] = fid;
                    new_nodes2.Set(face,fid);
                    weights_in_sim.Append(PARENT(parent_sim_tet, c));
                }
            }
            
            for (int j = 0; j < tc.material_ids.m; j++) {
                int mid = tc.material_ids(j);
                int fid = mid / 6;
                new_elements.Append(VECTOR<int,4>(element(face2node[fid][((mid%6+1)%6)/2]), cid, fc[fid], ec[face2edge[fid][(mid%6)/2]]));
                new_ctet2stet.Append(parent_sim_tet);
                new_is_blue.Append(is_blue(i));
            }
            
//            VECTOR<bool,4> has_node;
//            for (int j = 0; j < tc.material_ids.m; j++) {
//                has_node(material2node[tc.material_ids(j)]) = 1;
//            }
//            for (int j = 0; j < 4; j++){
//                if(has_node(j)) {
//                    for (int k = 0; k < 3; k++){
//                        new_elements.Append(VECTOR<int,4>(element(j), cid, fc[node2face[j][k]], ec[node2edge[j][k][0]]));
//                        new_elements.Append(VECTOR<int,4>(element(j), cid, fc[node2face[j][k]], ec[node2edge[j][k][1]]));
//                        cuttingFaces.Set(VECTOR<int,3>(cid, fc[node2face[j][k]], ec[node2edge[j][k][0]]).Sorted());
//                        cuttingFaces.Set(VECTOR<int,3>(cid, fc[node2face[j][k]], ec[node2edge[j][k][1]]).Sorted());
//                        new_ctet2stet.Append(parent_sim_tet);
//                        new_ctet2stet.Append(parent_sim_tet);
//                        new_is_blue.Append(is_blue(i));
//                        new_is_blue.Append(is_blue(i));
//                    }
//                }
//            }  
        }
    }
    
    //subdivide by face intersections
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        if (tet_cuttings(i).material_ids.m == MaxPieces){
            const int parent_sim_tet = ctet2stet(i);
            TET element = volume->mesh.elements(i);
            ARRAY<TET> new_e;
            new_e.Append(element);
            
            int nn = -1;
            HASHTABLE<TV_INT, int> ht;
            for (int j = 0; j < NumFacesPerTet; j++){
                TV_INT face = TV_INT(element(face2node[j][0]),element(face2node[j][1]),element(face2node[j][2])).Sorted();
                int new_node = -1;//subdivide the tet by this node
                if (new_nodes2.Get(face, new_node)){
                    if(new_e.m == 1) {//subdivide the original tet
                        new_e.Resize(NumNodesPerTriangle);
                        for (int k = 0; k < NumNodesPerTriangle; k++){
                            new_e(k) = TET(new_node, element(face2node[j][k]), element(face2node[j][(k+1)%NumNodesPerTriangle]), element(face2opposite_node[j]));
                            ht.Set(TV_INT(element(face2node[j][k]), element(face2node[j][(k+1)%NumNodesPerTriangle]), element(face2opposite_node[j])).Sorted(), k);
                        }
                        nn = new_node;
                    }
                    else {
                        int k;
                        if (!ht.Get(face, k)) 
                            cout << "should have been set!!!\n";
                        else {
                            new_e(k) = TET(new_node, element(face2node[j][0]), element(face2node[j][1]), nn);
                            new_e.Append(TET(new_node, element(face2node[j][1]), element(face2node[j][2]), nn));
                            new_e.Append(TET(new_node, element(face2node[j][2]), element(face2node[j][0]), nn));
                        }
                    }
                }
            }
            
            for (int j = 0; j < new_e.m; j++){
                new_ctet2stet.Append(parent_sim_tet);
                new_elements.Append(new_e(j));
                new_is_blue.Append(is_blue(i));
            }
        } 
    }
    
    //subdivide by edge intersections
    for (int i = 0; i < new_elements.m; i++) {
        VECTOR<int,4> element = new_elements(i);
        for (int j = 0; j < 6; ++j) {
            VECTOR<int,2> edge = VECTOR<int,2>(element(edge2node[j][0]),element(edge2node[j][1])).Sorted();
            if (new_nodes.Contains(edge)) {
                edge2tets.Get_Or_Insert(edge).insert(i);
            }
        }
    }
    
    for(typename HASHTABLE<VECTOR<int,2>, std::set<int> >::ITERATOR iterator(edge2tets);iterator.Valid();iterator.Next()){
        VECTOR<int,2> edge = iterator.Key();
        int e1 = edge(0);
        int e2 = edge(1);
        int nn = new_nodes.Get(edge);
        set<int>& tets = iterator.Data();
        for (set<int>::iterator it = tets.begin(); it != tets.end(); ++it) {
            int tid = *it;
            VECTOR<int,4> element = new_elements(tid);
            int n1 = -1, n2 = -1;
            for (int i = 0; i < 4; ++i) {
                if (element(i) != e1 && element(i) != e2) {
                    if (n1 < 0) {
                        n1 = element(i);
                    }
                    else {
                        n2 = element(i);
                    }
                }
            }

            new_elements.Append(TET(n1, n2, nn, e2));
            new_ctet2stet.Append(new_ctet2stet(tid));
            new_is_blue.Append(new_is_blue(tid));
            new_elements(tid) = TET(n1, n2, nn, e1);
            
            VECTOR<int, 2> ee = VECTOR<int, 2>(n1, n2).Sorted();
            if (edge2tets.Contains(ee)) {
                edge2tets.Get(ee).insert(new_elements.m-1);
            }
            ee = VECTOR<int,2>(n1, e2).Sorted();
            if (edge2tets.Contains(ee)) {
                edge2tets.Get(ee).erase(tid);
                edge2tets.Get(ee).insert(new_elements.m-1);
            }
            ee = VECTOR<int,2>(n2, e2).Sorted();
            if (edge2tets.Contains(ee)) {
                edge2tets.Get(ee).erase(tid);
                edge2tets.Get(ee).insert(new_elements.m-1);
            }
        }
    }

    volume->mesh.elements = new_elements;
    tet_cuttings.Remove_All();
    tet_cuttings.Resize(volume->mesh.elements.m);
    ctet2stet = new_ctet2stet;
    is_blue = new_is_blue;
    volume->Update_Number_Nodes();
    Update_Cutting_Particles();
    cutting_particle_material_space.Resize(volume->particles.X.m);
    for(int i=0;i<volume->particles.X.m;++i)
        cutting_particle_material_space(i)=volume->particles.X(i);
}

template<class T>
void MESH_CUTTING<T>::Refine_And_Save_To(TETRAHEDRALIZED_VOLUME<T>* refined_volume)
{
    HASHTABLE<VECTOR<int,2>,int> new_nodes;//edge to edge intersection point's index
    HASHTABLE<VECTOR<int,2>, std::set<int> > edge2tets;
    HASHTABLE<VECTOR<int,3>,int> new_nodes2;
    ARRAY<TET> new_elements;
    ARRAY<PARENT> new_weights_in_sim = weights_in_sim;
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        if (tet_cuttings(i).material_ids.m < MaxPieces){
            const int parent_sim_tet = ctet2stet(i);
            VECTOR<int,4> &element = volume->mesh.elements(i);
            TET_CUTTING& tc = tet_cuttings(i);
            
            //tet center
            int cid = new_weights_in_sim.m;
            {
                CENTER c = Weight_In_Sim_Tet(tc.tet_center.value, element, parent_sim_tet);
                new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
            }
            
            //edge centers
            VECTOR<int, NumEdgesPerTet> ec;
            for (int j = 0; j < NumEdgesPerTet; j++){
                CENTER c = Weight_In_Sim_Tet(tc.edge_centers[j].value, element, parent_sim_tet);
                VECTOR<int,2> edge = VECTOR<int,2>(element(edge2node[j][0]),element(edge2node[j][1])).Sorted();
                //store the weight
                if (!(new_nodes.Get(edge,ec[j]))){
                    int eid = new_weights_in_sim.m;
                    ec[j] = eid;
                    new_nodes.Set(edge,eid);
                    new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
                }
            }
            
            //face centers
            VECTOR<int, NumFacesPerTet> fc;
            for (int j = 0; j < NumFacesPerTet; j++){
                CENTER c = Weight_In_Sim_Tet(tc.face_centers[j].value, element, parent_sim_tet);
                TV_INT face = TV_INT(element(face2node[j][0]),element(face2node[j][1]),element(face2node[j][2])).Sorted();
                if (!(new_nodes2.Get(face,fc[j]))){
                    int fid = new_weights_in_sim.m;
                    fc[j] = fid;
                    new_nodes2.Set(face,fid);
                    new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
                }
            }
            
            VECTOR<bool,4> has_node;
            for (int j = 0; j < tc.material_ids.m; j++) {
                has_node(material2node[tc.material_ids(j)]) = 1;
            }
            for (int j = 0; j < 4; j++){
                if(has_node(j)) {
                    for (int k = 0; k < 3; k++){
                        new_elements.Append(VECTOR<int,4>(element(j), cid, fc[node2face[j][k]], ec[node2edge[j][k][0]]));
                        new_elements.Append(VECTOR<int,4>(element(j), cid, fc[node2face[j][k]], ec[node2edge[j][k][1]]));
                    }
                }
            }
        }
    }
    
    //subdivide by face intersections
    for (int i = 0; i < volume->mesh.elements.m; i++) {
        if (tet_cuttings(i).material_ids.m == MaxPieces){
            TET element = volume->mesh.elements(i);
            ARRAY<TET> new_e;
            new_e.Append(element);
            
            int nn = -1;
            HASHTABLE<TV_INT, int> ht;
            for (int j = 0; j < NumFacesPerTet; j++){
                TV_INT face = TV_INT(element(face2node[j][0]),element(face2node[j][1]),element(face2node[j][2])).Sorted();
                int new_node = -1;//subdivide the tet by this node
                if (new_nodes2.Get(face, new_node)){
                    if(new_e.m == 1) {//subdivide the original tet
                        new_e.Resize(NumNodesPerTriangle);
                        for (int k = 0; k < NumNodesPerTriangle; k++){
                            new_e(k) = TET(new_node, element(face2node[j][k]), element(face2node[j][(k+1)%NumNodesPerTriangle]), element(face2opposite_node[j]));
                            ht.Set(TV_INT(element(face2node[j][k]), element(face2node[j][(k+1)%NumNodesPerTriangle]), element(face2opposite_node[j])).Sorted(), k);
                        }
                        nn = new_node;
                    }
                    else {
                        int k;
                        if (!ht.Get(face, k))
                            cout << "should have been set!!!\n";
                        else {
                            new_e(k) = TET(new_node, element(face2node[j][0]), element(face2node[j][1]), nn);
                            new_e.Append(TET(new_node, element(face2node[j][1]), element(face2node[j][2]), nn));
                            new_e.Append(TET(new_node, element(face2node[j][2]), element(face2node[j][0]), nn));
                        }
                    }
                }
            }
            
            for (int j = 0; j < new_e.m; j++){
                new_elements.Append(new_e(j));
            }
        }
    }
    
    //subdivide by edge intersections
    for (int i = 0; i < new_elements.m; i++) {
        VECTOR<int,4> element = new_elements(i);
        for (int j = 0; j < 6; ++j) {
            VECTOR<int,2> edge = VECTOR<int,2>(element(edge2node[j][0]),element(edge2node[j][1])).Sorted();
            if (new_nodes.Contains(edge)) {
                edge2tets.Get_Or_Insert(edge).insert(i);
            }
        }
    }
    
    for(typename HASHTABLE<VECTOR<int,2>, std::set<int> >::ITERATOR iterator(edge2tets);iterator.Valid();iterator.Next()){
        VECTOR<int,2> edge = iterator.Key();
        int e1 = edge(0);
        int e2 = edge(1);
        int nn = new_nodes.Get(edge);
        set<int>& tets = iterator.Data();
        for (set<int>::iterator it = tets.begin(); it != tets.end(); ++it) {
            int tid = *it;
            VECTOR<int,4> element = new_elements(tid);
            int n1 = -1, n2 = -1;
            for (int i = 0; i < 4; ++i) {
                if (element(i) != e1 && element(i) != e2) {
                    if (n1 < 0) {
                        n1 = element(i);
                    }
                    else {
                        n2 = element(i);
                    }
                }
            }
            
            new_elements.Append(TET(n1, n2, nn, e2));
            new_elements(tid) = TET(n1, n2, nn, e1);
            
            VECTOR<int, 2> ee = VECTOR<int, 2>(n1, n2).Sorted();
            if (edge2tets.Contains(ee)) {
                edge2tets.Get(ee).insert(new_elements.m-1);
            }
            ee = VECTOR<int,2>(n1, e2).Sorted();
            if (edge2tets.Contains(ee)) {
                edge2tets.Get(ee).erase(tid);
                edge2tets.Get(ee).insert(new_elements.m-1);
            }
            ee = VECTOR<int,2>(n2, e2).Sorted();
            if (edge2tets.Contains(ee)) {
                edge2tets.Get(ee).erase(tid);
                edge2tets.Get(ee).insert(new_elements.m-1);
            }
        }
    }
    
    refined_volume->particles.Resize(new_weights_in_sim.m);
    for (int i = 0; i < new_weights_in_sim.m; i++){
        refined_volume->particles.X(i) = weight2vec_sim(new_weights_in_sim(i).id, new_weights_in_sim(i).weight);
    }
    refined_volume->Update_Number_Nodes();
    refined_volume->mesh.elements = new_elements;
}

template<class T> 
ARRAY<int> MESH_CUTTING<T>::Find_Cutting_Tet_CC(int m, ARRAY<bool>& picked, ARRAY<ARRAY<int> >& itot, ARRAY<ARRAY<int> >& node_cc)
{
    ARRAY<bool> visited(volume->particles.X.m);
    ARRAY<int> cc, n_cc;
    cc.Append(m);
    picked(m) = 1;
    for (int i = 0; i < cc.m; i++){
        int r = cc(i);
        for (int j = 0; j < NumNodesPerTet; j++) {
            int node_id = volume->mesh.elements(r)(j);
            if (!visited(node_id)){
                visited(node_id) = 1;
                n_cc.Append(node_id);
                ARRAY<int>& v = itot(volume->mesh.elements(r)(j));
                int n = v.m;
                for (int k = 0; k < n; k++){
                    if (!picked(v(k))) {
                        cc.Append(v(k));
                        picked(v(k)) = 1; }}}}}
    node_cc.Append(n_cc);
    return cc;
}

template<class T> 
ARRAY<int> MESH_CUTTING<T>::Find_Tet_CC(int m, ARRAY<bool>& picked, ARRAY<ARRAY<int> >& itot, ARRAY<ARRAY<int> >& node_cc)
{
    ARRAY<bool> visited(sim_volume->particles.X.m);
    ARRAY<int> cc, n_cc;
    cc.Append(m);
    picked(m) = 1;
    for (int i = 0; i < cc.m; i++){
        int r = cc(i);
        for (int j = 0; j < NumNodesPerTet; j++) {
            int node_id = sim_volume->mesh.elements(r)(j);
            if (!visited(node_id)){
                visited(node_id) = 1;
                n_cc.Append(node_id);
                ARRAY<int>& v = itot(sim_volume->mesh.elements(r)(j));
                int n = v.m;
                for (int k = 0; k < n; k++){
                    if (!picked(v(k))) {
                        cc.Append(v(k));
                        picked(v(k)) = 1; }}}}}
    node_cc.Append(n_cc);
    return cc;
}

//for drawing
template<class T> 
void MESH_CUTTING<T>::Update_For_Draw()
{
    cout << "Updating for draw...\n";
    sim_volume->Update_Number_Nodes();
    sim_volume->mesh.Initialize_Boundary_Mesh(); //cout << "sim boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
    sim_volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
    
    Fix_Orientation();
    volume->Update_Number_Nodes();
    volume->mesh.Initialize_Boundary_Mesh(); //cout << "cutting boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
    volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
    
    tet_cc.Remove_All();
    node_cc.Remove_All();
    int num_elements = sim_volume->mesh.elements.m;
    cout << "total sim elements: " << num_elements << endl;
    cout << "total sim particles: " << sim_volume->particles.X.m << endl;
    //find tet connected components
    ARRAY<bool> picked(num_elements);
    ARRAY<ARRAY<int> > itot(sim_volume->particles.X.m); //node index to tet index
    for (int i = 0; i < num_elements; i++) {
            VECTOR<int, NumNodesPerTet> element = sim_volume->mesh.elements(i);
            for (int j = 0; j < NumNodesPerTet; j++) {
                itot(element(j)).Append(i); }}
            
    for (int i = 0; i < num_elements; i++) {
        if (!picked(i)) {
            tet_cc.Append(Find_Tet_CC(i, picked, itot, node_cc)); }}
    cout << tet_cc.m << " CCs" << endl;      
      
    element2cc.Resize(sim_volume->mesh.elements.m);
    for (int i = 0; i < tet_cc.m; i++){
        for (int j = 0; j < tet_cc(i).m; j++) {
            element2cc(tet_cc(i)(j)) = i;
        }
    }
    
    node2cc.Resize(sim_volume->particles.X.m);
    for (int i = 0; i < node_cc.m; i++){
        for (int j = 0; j < node_cc(i).m; j++) {
            node2cc(node_cc(i)(j)) = i;
        }
    }

    //faces and interfaces to be drawn
    start_timer();
    faces_to_draw.Remove_All();
    interfaces_to_draw.Remove_All();
    HASHTABLE<TV_INT,bool> htf;
    bool b = 0;
    for (int i = 0; i < volume->mesh.elements.m; i++){
        //faces
        for (int j = 0; j < 4; j++){
            TV_INT f(volume->mesh.elements(i)[face2node[j][0]], volume->mesh.elements(i)[face2node[j][1]], volume->mesh.elements(i)[face2node[j][2]]);
            TV_INT ff = f.Sorted();
            if (htf.Get(ff,b)){
                htf.Set(ff,0);
            }
            else {
                htf.Set(ff,1);
            }
        }
        
        //interfaces
        for (int j = 0; j < NumSplitTurnOns; j++){
            if (tet_cuttings(i).has_material[face2material[j][0]] ^ tet_cuttings(i).has_material[face2material[j][1]]){
                interfaces_to_draw.Append(VECTOR<int,2>(i,j));
            }
        }
    }
    
    for (int i = 0; i < volume->mesh.elements.m; i++){
        //faces
        for (int j = 0; j < 4; j++){
            htf.Get(TV_INT(volume->mesh.elements(i)[face2node[j][0]], volume->mesh.elements(i)[face2node[j][1]], volume->mesh.elements(i)[face2node[j][2]]).Sorted(),b);
            if (b){
                faces_to_draw.Append(VECTOR<int,2>(i,j));
            }
        }
    }
    stop_timer();
    printf("update for draw time:    %f\n",get_time());
}

template<typename T> 
inline typename MESH_CUTTING<T>::TV& MESH_CUTTING<T>::Position(const int element_id, const int node_id) const
{
    return volume->particles.X(volume->mesh.elements(element_id)(node_id));
}

template<typename T> 
inline int MESH_CUTTING<T>::Index(const int element_id, const int node_id) const
{
    return volume->mesh.elements(element_id)(node_id);
}
//need modification
template<class T> 
int MESH_CUTTING<T>::Compute_Intersection(const T& x, const T& y)
{
    T z = 100;
    T zz;
    VECTOR<TV,3> v;
    int element_id = -1;
    for (int i = 0; i < faces_to_draw.m; i++){
        //find intersection point and its weight in the tet
        ARRAY<TV> mesh_vec, color_vec, refine_edges, normal_vec;
        Faces_For_OpenGL(faces_to_draw(i), mesh_vec, color_vec, normal_vec, refine_edges, picking);
        for (int j = 0; j < mesh_vec.m/3; j++) {
            if (intersects(x, y, mesh_vec(3*j), mesh_vec(3*j+1), mesh_vec(3*j+2), zz)){
                if (zz < z) {
                    z = zz;
                    element_id = faces_to_draw(i)(0);
                    v = VECTOR<TV,3>(mesh_vec(3*j), mesh_vec(3*j+1), mesh_vec(3*j+2));}}}} 

    for (int i = 0; i < interfaces_to_draw.m; i++){
        //find intersection point
        ARRAY<TV> mesh_vec, refine_edges;
        Interfaces_OpenGL(interfaces_to_draw(i), mesh_vec, refine_edges);
        if (intersects(x, y, mesh_vec(0), mesh_vec(1), mesh_vec(2), zz)){
            if (zz < z) {
                z = zz;
                element_id = interfaces_to_draw(i)(0);
                v = VECTOR<TV,3>(mesh_vec(0), mesh_vec(1), mesh_vec(2));}}}
    
    if (element_id >= 0) {
        int sim_tet = ctet2stet(element_id);
        int id;
        if (dragging_elements.Find(sim_tet, id)) {
            return id;
        }
        else {
            dragging_elements.Append(sim_tet);
            VECTOR<T,3> cc = weight_in_tet(TV(x, y, z), Position(element_id,0), Position(element_id,1), Position(element_id,2), Position(element_id,3));
            CENTER c(cc[0],cc[1],cc[2],1-cc[0]-cc[1]-cc[2]); 
            CENTER w = Weight_In_Sim_Tet(c, volume->mesh.elements(element_id), sim_tet);
            dragging_weights.Append(TV(w[0],w[1],w[2]));
            dragging_targets.Append(TV(x,y,z));
            cout << "dragging sim_tet " << sim_tet << endl;
            return dragging_elements.m-1;
        }
    }
    return -1;
}

template<class T> 
typename MESH_CUTTING<T>::TV MESH_CUTTING<T>::weight2vec_sim(int tet_id, CENTER c) const
{
    TV v;
    for (int i = 0; i < NumNodesPerTet; i++) {
        v += sim_volume->particles.X(sim_volume->mesh.elements(tet_id)(i))*c[i];
    }
    return v;
}

template<class T> 
typename MESH_CUTTING<T>::TV MESH_CUTTING<T>::weight2vec(int tet_id, CENTER c) const
{
    TV v;
    for (int i = 0; i < NumNodesPerTet; i++) {
        v += volume->particles.X(volume->mesh.elements(tet_id)(i))*c[i];
    }
    return v;
}

double MAG = 1.001;
template<class T> 
void MESH_CUTTING<T>::Colors_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& color_vec, const DRAWING_MODE& mode) const
{
    if (mode == normal){
//        if (dragging_elements.Contains(ctet2stet(f[0])))
//            color_vec.Append(TV(0,1,1));
        if (is_blue(f[0]))
            color_vec.Append(TV(0,0,1));
        else
            color_vec.Append(TV(1,1,1));
    }
    else if (mode == picking){
 cout <<       (T)element2cc(ctet2stet(f[0]))/256 << endl;
        color_vec.Append(TV((T)element2cc(ctet2stet(f[0]))/256,0,0));
    } 
}

template<class T> 
void MESH_CUTTING<T>::Faces_For_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& mesh_vec, ARRAY<TV>& color_vec, ARRAY<TV>& normal_vec, ARRAY<TV>& refine_edges, const DRAWING_MODE& mode) const
{   
    bool b = 1;
    for (int i = 0; i < NumMaterialsPerFace; i++){
        if (!tet_cuttings(f[0]).has_material(f[1]*6+i)){
            b = 0;
            break;
        }
    }
    if (b) {
        //adjust the ordering of nodes so face normal will point outwards
        TV_INT a;
        if (f[1] == 0) {
            a[0] = face2node[f[1]][1];
            a[1] = face2node[f[1]][0];
            a[2] = face2node[f[1]][2];
        }
        else {
            a[0] = face2node[f[1]][0];
            a[1] = face2node[f[1]][1];
            a[2] = face2node[f[1]][2];
        }
        
        for (int i = 0; i < NumNodesPerTriangle; i++) {
            mesh_vec.Append(volume->particles.X(volume->mesh.elements(f[0])(a[i])));
            Colors_OpenGL(f, color_vec, mode);
        }
        
        TV normal = TRIANGLE_3D<T>(volume->particles.X(volume->mesh.elements(f[0])(a[0])),
                                   volume->particles.X(volume->mesh.elements(f[0])(a[1])),
                                   volume->particles.X(volume->mesh.elements(f[0])(a[2]))).Normal();
        normal_vec.Append(normal);
        normal_vec.Append(normal);
        normal_vec.Append(normal);
    }
    else {
        for (int j = 0; j < NumMaterialsPerFace; j++){
            if (tet_cuttings(f[0]).has_material(f[1]*6+j)){
                int node_id = face2node[f[1]][((j+1)/2)%3];
                mesh_vec.Append(volume->particles.X(volume->mesh.elements(f[0])(node_id)));
                
                int edge_id = face2edge[f[1]][j/2];
                mesh_vec.Append(weight2vec(f[0], tet_cuttings(f[0]).edge_centers[edge_id].value));
                
                mesh_vec.Append(weight2vec(f[0], tet_cuttings(f[0]).face_centers[f[1]].value));
                
                refine_edges.Append(volume->particles.X(volume->mesh.elements(f[0])(node_id))*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).edge_centers[edge_id].value)*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).edge_centers[edge_id].value)*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).face_centers[f[1]].value)*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).face_centers[f[1]].value)*MAG);
                refine_edges.Append(volume->particles.X(volume->mesh.elements(f[0])(node_id))*MAG);
                
                for (int i = 0; i < NumNodesPerTriangle; i++) {
                    Colors_OpenGL(f, color_vec, mode);
                }
            }
        }
    }
}

template<class T> 
void MESH_CUTTING<T>::Interfaces_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& mesh_vec, ARRAY<TV>& refine_edges) const
{
    TV v1 = weight2vec(f[0], tet_cuttings(f[0]).tet_center.value);
    mesh_vec.Append(v1);
    refine_edges.Append(v1);
    
    if (f[1] < 24) {
        int face_id = f[1]/6;
        TV v2 = weight2vec(f[0], tet_cuttings(f[0]).face_centers[face_id].value);
        mesh_vec.Append(v2);
        refine_edges.Append(v2);
        refine_edges.Append(v2);
        
        int id_within_face = (f[1]%6)/2;
        if ((f[1]%6) & 1){
            TV v3 = weight2vec(f[0], tet_cuttings(f[0]).edge_centers[face2edge[face_id][id_within_face]].value);
            mesh_vec.Append(v3);
            refine_edges.Append(v3);
            refine_edges.Append(v3);
        }
        else {
            TV v3 = volume->particles.X(volume->mesh.elements(f[0])(face2node[face_id][id_within_face]));
            mesh_vec.Append(v3);
            refine_edges.Append(v3);
            refine_edges.Append(v3);
        }
    }
    
    refine_edges.Append(v1);
}

template<class T> 
void MESH_CUTTING<T>::Interfaces_For_OpenGL(const VECTOR<int,2>& f, ARRAY<TV>& mesh_vec, ARRAY<TV>& refine_edges, ARRAY<TV>& color_vec, const DRAWING_MODE& mode) const
{
    Interfaces_OpenGL(f, mesh_vec, refine_edges);
    
    for (int i = 0; i < NumNodesPerTriangle; i++) {
        Colors_OpenGL(f, color_vec, mode);
    }
}

template<class T> 
void MESH_CUTTING<T>::Draw(bool draw_cutting_mesh, const DRAWING_MODE& mode) const
{
    ARRAY<typename OPENGL_POLICY<T>::T_GL> vertices;
    start_timer();
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    
    //material
    glEnableClientState(GL_VERTEX_ARRAY);
//    for(int t=0;t<volume->mesh.boundary_mesh->elements.m;t++){
//        int i,j,k;volume->mesh.boundary_mesh->elements(t).Get(i,j,k);
//        VECTOR<T,3> xi=volume->particles.X(i),xj=volume->particles.X(j),xk=volume->particles.X(k);
//        OpenGL_Triangle(xi,xj,xk,vertices);}  
//    glVertexPointer(3, GL_DOUBLE,0,vertices.base_pointer);
//    glColor4f(1.0, 0.0, 0.0, 1.0);
//    cout << vertices.m << endl;
//    glDrawArrays(GL_TRIANGLES,0,vertices.m/3);
    
    ARRAY<TV> mesh_vec, color_vec, refine_edges, normal_vec;
    for (int i = 0; i < faces_to_draw.m; i++) {
        Faces_For_OpenGL(faces_to_draw(i), mesh_vec, color_vec, normal_vec, refine_edges, mode);
    }
    
    for (int i = 0; i < interfaces_to_draw.m; i++) {
        Interfaces_For_OpenGL(interfaces_to_draw(i), mesh_vec, refine_edges, color_vec, mode);
    }
    
    
//    cout << mesh_vec << endl;
    glVertexPointer(3, GL_FLOAT,0,mesh_vec.base_pointer);
    glColorPointer(3, GL_FLOAT, 0, color_vec.base_pointer);
    glNormalPointer(GL_FLOAT, 0, normal_vec.base_pointer);
    cout << mesh_vec.m << " color: " << color_vec.m << endl;
    glDrawArrays(GL_TRIANGLES,0,mesh_vec.m);
    
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    
    if (mode == normal){
        //sim edges
        if(0) {
            vertices.Remove_All();
            for(int t=0;t<sim_volume->mesh.boundary_mesh->segment_mesh->elements.m;t++){
                int i,j;sim_volume->mesh.boundary_mesh->segment_mesh->elements(t).Get(i,j);
                VECTOR<T,3> xi=sim_volume->particles.X(i),xj=sim_volume->particles.X(j);
                OpenGL_Line(xi*MAG,xj*MAG,vertices);}  
            glVertexPointer(3, GL_FLOAT,0,vertices.base_pointer);
            glColor4f(0.0, 0.0, 1.0, 1.0);
            glDrawArrays(GL_LINES,0,vertices.m/3);
        }
        
        //cutting mesh edges
        if (1){
            vertices.Remove_All();
            for(int t=0;t<volume->mesh.boundary_mesh->segment_mesh->elements.m;t++){
                int i,j;volume->mesh.boundary_mesh->segment_mesh->elements(t).Get(i,j);
                VECTOR<T,3> xi=volume->particles.X(i),xj=volume->particles.X(j);
                OpenGL_Line(xi*MAG,xj*MAG,vertices);}  
            glVertexPointer(3, GL_FLOAT,0,vertices.base_pointer);
            glColor4f(0.0, 1.0, 1.0, 1.0);
            glDrawArrays(GL_LINES,0,vertices.m/3);    
        }
        
        //refine_edges
        if(1){
            glVertexPointer(3, GL_FLOAT,0,refine_edges.base_pointer);
            glColor4f(0.0, 1.0, 1.0, 1.0);
            glDrawArrays(GL_LINES,0,refine_edges.m);
        }
        
        //cutting mesh
        if(draw_cutting_mesh){
            glVertexPointer(3, GL_FLOAT,0,cutting_vertices.base_pointer);
            glColor4f(0.0, 1.0, 0.0, 1.0);
            glDrawArrays(GL_LINES,0,cutting_vertices.m);
        }
        glutSwapBuffers();  
    }
    
    glDisableClientState(GL_VERTEX_ARRAY);    
    stop_timer();
    //printf("drawing time:    %f\n",get_time());
}

template<class T>
void MESH_CUTTING<T>::Write_To_File(const string& writing_directory, int frame) const
{
    start_timer();
    stringstream ss;
    ss << frame;
    
    ARRAY<TV> mesh_vec, color_vec, normal_vec, refine_edges;
    for (int i = 0; i < faces_to_draw.m; i++) {
        Faces_For_OpenGL(faces_to_draw(i), mesh_vec, color_vec, normal_vec, refine_edges, normal);
    }

    for (int i = 0; i < interfaces_to_draw.m; i++) {
        Interfaces_For_OpenGL(interfaces_to_draw(i), mesh_vec, refine_edges, color_vec, normal);
    }
    
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/mesh"+ss.str()+string(".tet.gz"), mesh_vec);
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/color"+ss.str()+string(".tet.gz"), color_vec);
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/normal"+ss.str()+string(".tet.gz"), normal_vec);
    
    stop_timer();
    printf("Writing mcut to file time:    %f\n",get_time());
}

template<class T>
void MESH_CUTTING<T>::Write_Boundary_Mesh_To_File(const string& writing_directory, int frame) const
{
    start_timer();
    stringstream ss;
    ss << frame;

    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/cutting_volume"+ss.str()+string(".tet.gz"), volume->particles, volume->mesh.boundary_mesh);
    
    stop_timer();
    printf("Writing boundary mesh to file time:    %f\n",get_time());
}

template<class T> 
void MESH_CUTTING<T>::Draw_For_Picking()
{
    start_timer();
    Draw(0, picking);
    stop_timer();
    printf("draw for picking time:    %f\n",get_time());
}

template<class T> 
void MESH_CUTTING<T>::Translate_CC(int picked_cc_id, TV translation)
{
    for (int i = 0; i < my_constrained->Size(); i++) {
        if (node2cc(my_constrained->operator()(i)/3) == picked_cc_id) {
            my_constrained_locations->operator()(i) += translation(i%3);
        }
    }
    be->Set_Boundary_Conditions(*my_constrained, *my_constrained_locations);
    
    for (int i = 0; i < node_cc(picked_cc_id).m; i++) {
        sim_volume->particles.X(node_cc(picked_cc_id)(i)) += translation;
    }
    //update cutting nodes
    Update_Cutting_Particles();
}

template<class T> 
void MESH_CUTTING<T>::Split_By_Levelset(const int& tet_id, ARRAY<int>& lsi)
{
//    TET_CUTTING tc = tet_cuttings(tet_id);
//    
////    for(int i = 0; i < 4; i++){ 
////        if(fabs(levelset(original_elements(tet_id)(i))) <= 10e-3){
////            levelset(original_elements(tet_id)(i)) = -0.01;
////////            if(levelset(original_elements(tet_id)(i)) == 0)
////////                cout << "00000000000000000\n" << endl;
////        }
////    }

////    if(levelsets(levelsets.m-1)(original_elements(tet_id)(0))<0 && levelsets(levelsets.m-1)(original_elements(tet_id)(1))<0 && levelsets(levelsets.m-1)(original_elements(tet_id)(2))<0 && levelsets(levelsets.m-1)(original_elements(tet_id)(3))<0) {
////        for (int i = 0; i < NumEdgesPerTet; i++) {
////            if(seed_colors(original_elements(tet_id)(edge2node[i][0])) != seed_colors(original_elements(tet_id)(edge2node[i][1]))){
////                tc.edge_centers[i].Append(Edge_Intersection(i,original_elements(tet_id)(edge2node[i][0]), original_elements(tet_id)(edge2node[i][1])));   
////                tc.edge_centers[i].Set();
////                tc.turned_on[edge2turnon[i][0]] = 1;
////                tc.turned_on[edge2turnon[i][1]] = 1;
////            }
////        }
////    }

//    
//    for (int ii = 0; ii < lsi.m; ii++){
//        ARRAY<T>& levelset = levelsets(lsi(ii));
//        int s = 0;
//        for (int i = 0; i < NumEdgesPerTet; i++) {
//            int n1 = edge2node[i][0];
//            int n2 = edge2node[i][1];
//            
//            int i1 = original_elements(tet_id)(n1);
//            int i2 = original_elements(tet_id)(n2);
//            
//            T& lv1 = levelset(i1);
//            T& lv2 = levelset(i2);
//            if (lv1!=lv1 || lv2!=lv2) cout << tet_id << "*********************************\n";
//            T lvm = lv1*lv2;
//            if(lvm<0){
//                s++;
//                CENTER c;
//                T sum = fabs(lv1) + fabs(lv2);
//                c[edge2node[i][0]] = fabs(lv2)/sum;
//                c[edge2node[i][1]] = fabs(lv1)/sum;
//                tc.edge_centers[i].Append(c);   
//                tc.edge_centers[i].Set();
//                tc.turned_on[edge2turnon[i][0]] = 1;
//                tc.turned_on[edge2turnon[i][1]] = 1;
//            }
//            else if (lv1 == 0) {
//                for(int j = 0; j < 3; j++){
//                    tc.turned_on[node2turnon[n1][j]] = 1;
//                }
//            }
//            else if (lv2 == 0) {
//                for(int j = 0; j < 3; j++){
//                    tc.turned_on[node2turnon[n2][j]] = 1;
//                }
//            }
//        }
//        
//        for (int i = 0; i < NumFacesPerTet; i++) {
//            for (int j = 0; j < NumEdgesPerTriangle; j++) {
//                if (tc.edge_centers[face2edge[i][j]].n == -1){
//                    tc.face_centers[i].Append(tc.edge_centers[face2edge[i][j]].value);
//                }
//            }
//        }
//        for (int i = 0; i < NumEdgesPerTet; i++) {
//            if (tc.edge_centers[i].n == -1){
//                tc.tet_center.Append(tc.edge_centers[i].value);
//            }
//        }
//    }
//    
//    //generate submaterials
//    //set centers if it's on a turned-on face that is split(only one piece is in a connected component)
//    //cout << tc.turned_on << endl;
//    VECTOR<bool, MaxPieces> picked;
//    int NumMaterials = tc.material_ids.m;
//    for (int i = 0; i < NumMaterials; i++) {
//        if (!picked[tc.material_ids(i)]) {
//            ARRAY<int> cc = tc.Find_CC(tc.material_ids(i), picked);
//            if (cc.m != NumMaterials) {
//                VECTOR<bool,4> has_mat;
//                for(int j = 0; j < cc.m; j++){
//                    has_mat(material2node[cc(j)]) = 1;
//                }
//                //cout << has_mat << endl;
//                for (int j = 0; j < NumNodesPerTet; j++){
//                    volume->particles.Append(volume->particles,volume->mesh.elements(tet_id)(j));
//                    if (has_mat(j)) {
//                        //cout << original_elements(tet_id)(j) << " has material\n" ;
//                        if (diri_nodes.find(DIRI(original_elements(tet_id)(j),0))!=diri_nodes.end()) {
//                            //cout << "find " << DIRI(original_elements(tet_id)(j),0) << " inserting " << DIRI(volume->particles.X.m-1,0) << endl;
//                            diri_nodes.insert(DIRI(volume->particles.X.m-1,0));
//                        }
//                        else if (diri_nodes.find(DIRI(original_elements(tet_id)(j),1))!=diri_nodes.end()) {
//                            //cout << "find " << DIRI(original_elements(tet_id)(j),1) << " inserting " << DIRI(volume->particles.X.m-1,1) << endl;
//                            diri_nodes.insert(DIRI(volume->particles.X.m-1,1));
//                        }
//                    }
//                }
//                int size = volume->particles.X.m;
//                if (i == 0) { 
//                    tet_cuttings(tet_id) = tc.Generate_Sub_Tet(cc); 
//                    volume->mesh.elements(tet_id) = VECTOR<int, NumNodesPerTet>(size-4, size-3, size-2, size-1);
//                }
//                else { 
//                    tet_cuttings.Append(tc.Generate_Sub_Tet(cc)); 
//                    volume->mesh.elements.Append(VECTOR<int, NumNodesPerTet>(size-4, size-3, size-2, size-1)); 
//                    //ALGEBRA::MATRIX_3X3<T> am = undeformed_config_copy(tet_id);
//                    //undeformed_config_copy.Append(am);
//                }
//            }
//            else {
//                tet_cuttings(tet_id) = tc; 
//            }
//        }
//    }
}

template<class T> 
void MESH_CUTTING<T>::Cut_By_Levelset(ARRAY<int>& lsi)
{
//    cout << "levelset cutting*************************************" << endl;
//    
//    start_timer();
//    //split
//    ARRAY<VECTOR<int, NumNodesPerTet> >& elements = volume->mesh.elements;

//    original_elements = volume->mesh.elements;
//    int num_elements = elements.m; cout << num_elements << " elements" << endl;
//    int num_nodes = volume->particles.number;
//    ARRAY<bool> need_dup(num_nodes);//whether the node needs to be duplicated
//    ARRAY<bool> need_merge(num_elements);//whether the element needs to be merged
//    ARRAY<int> pieces(num_elements);
//    int old_mesh_size = num_elements;
//    int split_count = 0;
//    for (int i = 0; i < num_elements; i++) {
//        Split_By_Levelset(i, lsi);
//        pieces(i) = elements.m - old_mesh_size;
//        if (pieces(i) > 0) {// in the box
//            //cout << "pices:" << pieces(i) << endl;
//            split_count++;
//            need_merge(i) = 1;
//            need_dup.Subset(original_elements(i)).Fill(1);
//            old_mesh_size = elements.m;
//        }
//    }
//    cout << "min: " << coefmin << ", max: " << coefmax << endl;
//    cout << split_count << " tets split\n";
//    for (int i = 0; i < num_elements; i++) {
//        if (!need_merge(i)) {//only nodes on need_merge(i.e. split) elements are duplicated
//            VECTOR<bool,4> has_mat;
//            for(int j = 0; j < tet_cuttings(i).material_ids.m; j++){
//                has_mat(material2node[tet_cuttings(i).material_ids(j)]) = 1;
//            }
//            for (int j = 0; j < NumNodesPerTet; j++) {
//                int& node_index = elements(i)(j);
//                if (need_dup(node_index)) {
//                    if (has_mat(j)) {
//                        if (diri_nodes.find(DIRI(node_index,0))!=diri_nodes.end()) {
//                            diri_nodes.insert(DIRI(volume->particles.X.m,0));
//                        }
//                        else if (diri_nodes.find(DIRI(node_index,1))!=diri_nodes.end()) {
//                            diri_nodes.insert(DIRI(volume->particles.X.m,1));
//                        }
//                    }
//                    node_index = volume->particles.Append(volume->particles,node_index);
//                    need_merge(i) = 1;
//                }
//            }
//        }
//    }
//                
//    stop_timer();
//    printf("split time:    %f\n",get_time());
//    cout << "total particles: " << volume->particles.X.m << endl;
//    //cout << elements << endl;

//    //merge nodes
//    //identify which un-plit tet need merging
//    start_timer();
//    
//    HASHTABLE<TV_INT, int> face_visited;
//    HASHTABLE<TV_INT, TV_TV_INT> oface2cface;//original face mapped to cut face
//    UNION_FIND<int> node_classes(volume->particles.number);//cluster merging nodes
//    TET_CUTTING tc;
//    VECTOR<int, NumNodesPerTet> element, original_element;
//    int piece_id = num_elements;
//    for (int i = 0; i < num_elements; i++) {
//        original_element = original_elements(i);
//        if (need_merge(i)) {
//            for (int k = 0; k < NumFacesPerTet; k++) {
//                TV_TV_INT cface, cface1;
//                tc = tet_cuttings(i);
//                element = elements(i);
//                TV_INT original_face(original_element(face2nodes[k][0]), original_element(face2nodes[k][1]), original_element(face2nodes[k][2]));
//                {
//                    for (int l = 0; l < NumMaterialsPerFace; l++) {
//                        if (tc.has_material(NumMaterialsPerFace*k+l)) {
//                            cface(l) = TV_INT(element(face2nodes[k][0]), element(face2nodes[k][1]), element(face2nodes[k][2]));
//                }}}
//                for (int j = 0; j < pieces(i); j++) {
//                    tc = tet_cuttings(piece_id + j);
//                    element = elements(piece_id + j);
//                    for (int l = 0; l < NumMaterialsPerFace; l++) {
//                        if (tc.has_material(NumMaterialsPerFace*k+l)) {
//                            cface(l) = TV_INT(element(face2nodes[k][0]), element(face2nodes[k][1]), element(face2nodes[k][2]));
//                    }}
//                }
//                bool found = 0;
//                for (int l = 0; l < NumNodesPerTriangle; l++) {
//                    TV_INT oface(original_face((2+l)%NumNodesPerTriangle), original_face((1+l)%NumNodesPerTriangle), original_face(l));
//                    if (oface2cface.Get(oface, cface1)){
//                        Union(cface, cface1, node_classes, l, 1);
//                        oface2cface.Set(oface, cface1);
//                        found = 1;
//                        break;
//                    }
//                    oface = TV_INT(original_face(l), original_face((1+l)%NumNodesPerTriangle), original_face((2+l)%NumNodesPerTriangle));
//                    if (oface2cface.Get(oface, cface1)){
//                        Union(cface, cface1, node_classes, l, 0);
//                        oface2cface.Set(oface, cface1);
//                        found = 1;
//                        break;
//                }}
//                if (!found) { oface2cface.Set(original_face, cface); }
//            }
//            piece_id+=pieces(i);
//        }
//    }
//    
//    int element_id = num_elements;
//    for (int i = 0; i < num_elements; i++) {
//        if (need_merge(i)) {
//            for (int k = 0; k < NumNodesPerTet; k++) {
//                elements(i)(k) = node_classes.Find(elements(i)(k));
//            }
//            for (int j = 0; j < pieces(i); j++) {
//                for (int k = 0; k < NumNodesPerTet; k++) {
//                    elements(element_id)(k) = node_classes.Find(elements(element_id)(k)); 
//                }
//                element_id++;
//    }}}
//    stop_timer();
//    printf("merge time:    %f\n",get_time());
//    cout << "total particles: " << volume->particles.X.m << endl;
//    //cout << elements << endl;
//    //cout << "tet_cc " << tet_cc << endl;
//    
//    //Refine();
//    
//    //Update_For_Draw();
//    //cout << "tet_cc " << tet_cc << endl;
//    
////    start_timer();
////    Reinitialize_Elasticity();
////    stop_timer();
////    printf("elasticity reinitialization time:    %f\n",get_time());
}

template<class T> 
void MESH_CUTTING<T>::Keep_CC(int cc_id)
{
    ARRAY<VECTOR<int,4> > new_elements;
    ARRAY<TET_CUTTING> new_cuttings;
    for (int i = 0; i < tet_cc(cc_id).m; i++){
        new_elements.Append(volume->mesh.elements(tet_cc(cc_id)(i)));
        new_cuttings.Append(tet_cuttings(tet_cc(cc_id)(i)));
    }
    volume->mesh.elements = new_elements;
    tet_cuttings = new_cuttings;
    Update_For_Draw();
}

template<class T> 
void MESH_CUTTING<T>::Become_Sphere()
{
    ARRAY<int> lsi; lsi.Append(0);
    Cut_By_Levelset(lsi);
    Update_For_Draw();
    Keep_CC(1);
}

template class MESH_CUTTING<float>;
template class MESH_CUTTING<double>;
//template class MESH_CUTTING<double>;
