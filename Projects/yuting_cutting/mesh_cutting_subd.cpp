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
#include <OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include "CONSISTENT_INTERSECTIONS.h"
#include <iostream>
#include "mesh_cutting_subd.h"
#include "DEFORMABLE_OBJECTS.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <set>
    
void read_tsc(){__asm__("rdtsc");}
struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);read_tsc();}
void stop_timer(){gettimeofday(&stoptime,NULL);read_tsc();}
double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}
#define DE cout<<"file "<<__FILE__<<"   line "<<__LINE__<<"  "<<&volume->particles.X<<"   "<<volume->particles.X<<endl;

using namespace PhysBAM;
using namespace std;

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

int edge2node[NumEdgesPerTet][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
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

int face2edge_center[24] = {-1,0,-1,1,-1,2,
                            -1,0,-1,4,-1,3,
                            -1,1,-1,5,-1,4,
                            -1,2,-1,3,-1,5};

int face2face_center[24] = {0,0,0,0,0,0,
                            1,1,1,1,1,1,
                            2,2,2,2,2,2,
                            3,3,3,3,3,3};

int tri_node2material[3][2] = {{0, 5}, {2, 1}, {4, 3}};

int face2nodes[NumFacesPerTet][NumNodesPerTriangle] = {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
int face2opposite_node[NumFacesPerTet] = {3, 2, 0, 1};

int node2face[NumNodesPerTet][3] = {{0,1,3},{0,1,2},{0,3,2},{1,2,3}};
int node2edge[NumNodesPerTet][3][2] = {{{0,2},{0,3},{2,3}},{{0,1},{0,4},{1,4}},{{1,2},{2,5},{1,5}},{{3,4},{4,5},{3,5}}};

int edge2turnon[NumEdgesPerTet][2] = {{1,7},{3,13},{5,19},{11,21},{9,17},{15,23}};
int node2turnon[NumNodesPerTet][3] = {{0,6,20},{2,8,12},{4,14,18},{10,16,22}};

template<class T>
void MESH_CUTTING<T>::Connected_Components(TETRAHEDRALIZED_VOLUME<T>* v, ARRAY<int>& labels) {
    int n = v->mesh.elements.m;
    labels.Resize(n);
    labels.Fill(0);
    UNION_FIND<int> uf(n);
    HASHTABLE<I3, int> ht;
    for (int i = 0; i < n; ++i) {
        I4 e = v->mesh.elements(i);
        for (int j = 0; j < 4; ++j) {
            I3 f(e(face2node[j][0]), e(face2node[j][1]), e(face2node[j][2]));
            f.Sort();
            int id = i;
            if (ht.Get(f, id)) {
                uf.Union(i, id);
            }
            else {
                ht.Set(f, id);
            }
        }
    }
    int c = 1;
    for (int i = 0; i < n; ++i) {
        int j = uf.Find(i);
        int l = labels(j);
        if (l > 0) {
            labels(i) = l;
        }
        else {
            labels(j) = c++;
            labels(i) = labels(j);
        }
    }
}

template<class T> 
MESH_CUTTING<T>::TET_CUTTING::TET_CUTTING()
{
    material_ids = IDENTITY_ARRAY<>(MaxPieces);
    for (int l = 0; l < material_ids.m; l++) { has_material(material_ids(l)) = 1; }
    tet_center.sum = CENTER::Constant_Vector(1.0/NumNodesPerTet);
    for (int i = 0; i < NumFacesPerTet; i++){
        for (int j = 0; j < NumNodesPerTriangle; j++){
            face_centers(i).sum(face2node[i][j]) = 1.0/3;
    }}
    for (int i = 0; i < NumEdgesPerTet; i++){
        edge_centers(i).sum(edge2node[i][0]) = 0.5;
        edge_centers(i).sum(edge2node[i][1]) = 0.5;
    }
}

template<class T> 
typename MESH_CUTTING<T>::TET_CUTTING MESH_CUTTING<T>::TET_CUTTING::Generate_Sub_Tet(const ARRAY<int>& material_ids_input)
{
    //set centers
    for (int j = 0; j < 24; j++){
        //has split along this face
        if(material_ids_input.Contains(face2material[j][0]) ^ material_ids_input.Contains(face2material[j][1])){
            tet_center.set = 1;
            face_centers[face2face_center[j]].set = 1;
            if (face2edge_center[j] >= 0) {
                edge_centers[face2edge_center[j]].set = 1;
            }
        }
    }
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
MESH_CUTTING<T>::MESH_CUTTING(TETRAHEDRALIZED_VOLUME<T>* sim_volume_input, T timestep_input, int ratio_input, bool interactive_input): sim_volume(sim_volume_input), interactive(interactive_input), timestep(timestep_input), ratio(ratio_input)
{
    Initialize_Cutting_Volume();
    undeformed_config_copy.Resize(sim_volume->mesh.elements.m);
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
                int c2s = ctet2stet(tets(i));
                ctet2stet.Append(c2s);
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
void MESH_CUTTING<T>::Initialize_Cutting_Volume()
{    
    //create cutting volume, initially copied from sim_volume
    volume = TETRAHEDRALIZED_VOLUME<T>::Create();    
    int initial_num_nodes = sim_volume->particles.X.m;
    weights_in_sim.Resize(initial_num_nodes);
    for (int i = 0; i < sim_volume->mesh.elements.m; i++){
        TET e = sim_volume->mesh.elements(i);
        ctet2stet.Append(i);
        for (int j = 0; j < NumNodesPerTet; j++){
            T4 c;
            c(j) = 1;
            weights_in_sim(e(j)) = PARENT(i, c);
        }
        volume->mesh.elements.Append(e);
    }
    tet_cuttings.Resize(sim_volume->mesh.elements.m);
    Update_Cutting_Particles();
    
    //refine cutting volume
    int num_refinements = 0;
    for (int i = 0; i < num_refinements; i++){
        Refine_Cutting_Volume();
    }
    Fix_Orientation();
    
    for(int i=0;i<volume->particles.X.m;++i)
        cutting_particle_material_space.Append(volume->particles.X(i));
        
    //subdevide cutting volume so it can be an eyeball...
    is_blue.Resize(volume->mesh.elements.m);
//    Subdivide_Cutting_Mesh_Into_Eyeball();
    
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
            CENTER c = Weight_In_Sim_Tet(tc.tet_center.Value(), element, parent_sim_tet);
            weights_in_sim.Append(PARENT(parent_sim_tet, c));
        }
        
        //edge centers
        VECTOR<int, NumEdgesPerTet> ec;
        for (int j = 0; j < NumEdgesPerTet; j++){
            CENTER c = Weight_In_Sim_Tet(tc.edge_centers[j].Value(), element, parent_sim_tet);
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
            CENTER c = Weight_In_Sim_Tet(tc.face_centers[j].Value(), element, parent_sim_tet);
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
    volume->Update_Number_Nodes();
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

    deformable_object = new DEFORMABLE_OBJECT_3D<ST>(SIZE_MESH, NUM_NODES, 0);
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
    fem = new FEM_HYPERELASTICITY_3D<ST>(deformable_object->Tetrahedron_Mesh(),deformable_object->Positions());
    le = new FIXED_COROTATED_ELASTICITY_3D<ST>((ST)100000,(ST).3,fem->F());
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
    be = new BACKWARD_EULER_TIME_STEPPING_3D<ST>(dt,end_time,start_time,*deformable_object);
    be->Set_Elastic_Forces(*fem);
    nodal_volumes = new ALGEBRA::VECTOR<ST>(deformable_object->Positions().Size());
    fem->Nodal_Volume_Fractions(*nodal_volumes);
    
    be->Initialize_BE_Matrix(*nodal_volumes);
    //be->Initialize_CG();
    be->Initialize_MINRES();
    
    if (1) {
        my_constrained = new ALGEBRA::VECTOR<int>(3*diri_nodes.Size());
        my_constrained_locations = new ALGEBRA::VECTOR<ST>(3*diri_nodes.Size());
        int i = 0;
        for (HASHTABLE_ITERATOR<int> it(diri_nodes); it.Valid(); it.Next()) {
            int fixed_node = it.Key();
            my_constrained->Set_Value(3*i,3*fixed_node); my_constrained->Set_Value(3*i+1,3*fixed_node+1); my_constrained->Set_Value(3*i+2,3*fixed_node+2); // first node constrained
            my_constrained_locations->Set_Value(3*i,sim_volume->particles.X(fixed_node)(0));
            my_constrained_locations->Set_Value(3*i+1,sim_volume->particles.X(fixed_node)(1));
            my_constrained_locations->Set_Value(3*i+2,sim_volume->particles.X(fixed_node)(2));
            i++;
        }
    }
    else {
        //z component of every node fixed, x,y component of all diri nodes fixed.
        my_constrained = new ALGEBRA::VECTOR<int>(2*diri_nodes.Size()+sim_volume->particles.X.m);
        my_constrained_locations = new ALGEBRA::VECTOR<ST>(2*diri_nodes.Size()+sim_volume->particles.X.m);
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
    }
    
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

    deformable_object = new DEFORMABLE_OBJECT_3D<ST>(SIZE_MESH, NUM_NODES, 0);
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
    fem = new FEM_HYPERELASTICITY_3D<ST>(deformable_object->Tetrahedron_Mesh(),deformable_object->Positions());
    le = new FIXED_COROTATED_ELASTICITY_3D<ST>((T)100000,(T).3,fem->F());
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
    be = new BACKWARD_EULER_TIME_STEPPING_3D<ST>(dt,end_time,start_time,*deformable_object);
    be->Set_Elastic_Forces(*fem);
    nodal_volumes = new ALGEBRA::VECTOR<ST>(deformable_object->Positions().Size());
    fem->Nodal_Volume_Fractions(*nodal_volumes);
    be->Initialize_BE_Matrix(*nodal_volumes);
    //be->Initialize_CG();
    be->Initialize_MINRES();
    
    delete my_constrained;
    delete my_constrained_locations;
    
    if (1) {
        my_constrained = new ALGEBRA::VECTOR<int>(3*diri_nodes.Size());
        my_constrained_locations = new ALGEBRA::VECTOR<ST>(3*diri_nodes.Size());
        int i = 0;
        for (HASHTABLE_ITERATOR<int> it(diri_nodes); it.Valid(); it.Next()) {
            int fixed_node = it.Key();
            my_constrained->Set_Value(3*i,3*fixed_node); my_constrained->Set_Value(3*i+1,3*fixed_node+1); my_constrained->Set_Value(3*i+2,3*fixed_node+2); // first node constrained
            my_constrained_locations->Set_Value(3*i,sim_volume->particles.X(fixed_node)(0));
            my_constrained_locations->Set_Value(3*i+1,sim_volume->particles.X(fixed_node)(1));
            my_constrained_locations->Set_Value(3*i+2,sim_volume->particles.X(fixed_node)(2));
            i++;
        }
    }
    else {
        //z component of every node fixed, x,y component of all diri nodes fixed.
        my_constrained = new ALGEBRA::VECTOR<int>(2*diri_nodes.Size()+sim_volume->particles.X.m);
        my_constrained_locations = new ALGEBRA::VECTOR<ST>(2*diri_nodes.Size()+sim_volume->particles.X.m);
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
    }
    be->Set_Boundary_Conditions(*my_constrained, *my_constrained_locations);
}

template<class T>
void MESH_CUTTING<T>::Split(const int& tet_id, HASHTABLE<int,H>& tri2inter, ARRAY<int>& sim_node_from, ARRAY<bool>& sim_tet_split, ARRAY<int>& material_node_from)
{
    //cout << "******" << tri2inter << endl;
    TET_CUTTING tc = tet_cuttings(tet_id);
    VECTOR<bool, NumSplitTurnOns> can_split;
    for (int i = 0; i < NumSplitTurnOns; i++){
        if(tc.has_material(face2material[i][0]) && tc.has_material(face2material[i][1])) {
            can_split(i) = 1;
    }}

    const VECTOR<int, NumNodesPerTet> tet_element = volume->mesh.elements(tet_id);
    for(typename HASHTABLE<int,H>::ITERATOR it(tri2inter);it.Valid();it.Next()) {
        H& intersects=it.Data();

        //centers will eventually be set as average of the related intersections, which works for cutting surface being approximately a plane inside the tet
        ARRAY<CENTER> vcenters[NumNodesPerTet];
        ARRAY<CENTER> ecenters[NumEdgesPerTet];
        ARRAY<CENTER> fcenters[NumFacesPerTet];
        ARRAY<CENTER> tcenters;
        
        for (typename H::ITERATOR it2(intersects); it2.Valid(); it2.Next()) {
            I4 t = it2.Key();
            T4 w = it2.Data();
            if (t(1) == -1) {
                T4 c;
                int id = tet_element.Find(t(0));
                c[id] = 1;
                vcenters[id].Append(c);
            }
            else if (t(2) == -1) {
                T4 c;
                c[tet_element.Find(t(0))] = w(0);
                c[tet_element.Find(t(1))] = 1 - w(0);
                for (int i = 0; i < 6; ++i) {
                    if (I2(tet_element(edge2node[i][0]),tet_element(edge2node[i][1])).Sorted() == I2(t(0),t(1))) {
                        ecenters[i].Append(c);
                        break;
                    }
                }
            }
            else if (t(3) == -1) {
                T4 c;
                for (int i = 0; i < 3; ++i) {
                    c(tet_element.Find(t(i))) = w(i);
                }
                for (int i = 0; i < 4; ++i) {
                    if (I3(tet_element(face2node[i][0]),tet_element(face2node[i][1]), tet_element(face2node[i][2])).Sorted() == t.Remove_Index(3)) {
                        fcenters[i].Append(c);
                        break;
                    }
                }
            }
            else {
                T4 c;
                for (int i = 0; i < 4; ++i) {
                    c(tet_element.Find(t(i))) = w(i);
                }
                tcenters.Append(c);
            }
        }
        
        //whether all centers are on one face. if so, flag them as no_merge and return
        bool cut_through = 1;
        for (int j = 0; j < 4; ++j) {
            int e1 = face2edge[j][0];
            int e2 = face2edge[j][1];
            int e3 = face2edge[j][2];
            int n1 = face2node[j][0];
            int n2 = face2node[j][1];
            int n3 = face2node[j][2];
            if (vcenters[n1].m+vcenters[n2].m+vcenters[n3].m+ecenters[e1].m+ecenters[e2].m+ecenters[e3].m+fcenters[j].m == intersects.Size()) {
//                cout << "no merge: " << tet_id << " " << tet_element << " " << j << " " << I3(tet_element(n1), tet_element(n2), tet_element(n3)) << endl;
//                cout << "ffff: " << intersects << endl;
//                cout <<vcenters[n1].m<<vcenters[n2].m<<vcenters[n3].m<<ecenters[e1].m<<ecenters[e2].m<<ecenters[e3].m<<fcenters[j].m << endl;
                bool cf = 0;
                if (vcenters[n1].m && (vcenters[n2].m || ecenters[e1].m) && (vcenters[n3].m || ecenters[e2].m || ecenters[e3].m || fcenters[j].m)) {
                    tc.no_merge(j*6) = 1;
                    cf = 1;
                }
                if (vcenters[n2].m && (vcenters[n3].m || ecenters[e2].m) && (vcenters[n1].m || ecenters[e1].m || ecenters[e3].m || fcenters[j].m)) {
                    tc.no_merge(j*6+2) = 1;
                    cf = 1;
                }
                if (vcenters[n3].m && (vcenters[n1].m || ecenters[e3].m) && (vcenters[n2].m || ecenters[e1].m || ecenters[e2].m || fcenters[j].m)) {
                    tc.no_merge(j*6+4) = 1;
                    cf = 1;
                }
                if (vcenters[n2].m && (vcenters[n1].m || ecenters[e1].m) && (vcenters[n3].m || ecenters[e2].m || ecenters[e3].m || fcenters[j].m)) {
                    tc.no_merge(j*6+1) = 1;
                    cf = 1;
                }
                if (vcenters[n3].m && (vcenters[n2].m || ecenters[e2].m) && (vcenters[n1].m || ecenters[e1].m || ecenters[e3].m || fcenters[j].m)) {
                    tc.no_merge(j*6+3) = 1;
                    cf = 1;
                }
                if (vcenters[n1].m && (vcenters[n3].m || ecenters[e3].m) && (vcenters[n2].m || ecenters[e1].m || ecenters[e2].m || fcenters[j].m)) {
                    tc.no_merge(j*6+5) = 1;
                    cf = 1;
                }
                cut_through = 0;
                if (cf) {
                    cutFaces.Set(I3(tet_element(n1), tet_element(n2), tet_element(n3)).Sorted());
                }
                break;
            }
        }
        if (cut_through) { //turn on flags and add the related centers
            //turn on edge segments
            for (int j = 0; j < NumEdgesPerTet; j++){
                if (vcenters[edge2node[j][0]].m && (ecenters[j].m || vcenters[edge2node[j][1]].m)) {
                    int id = 2 * j + 24;
                    tc.turned_on[id] = 1;
                }
                if (vcenters[edge2node[j][1]].m && (ecenters[j].m || vcenters[edge2node[j][0]].m)) {
                    int id = 2 * j + 1 + 24;
                    tc.turned_on[id] = 1;
                }
                for (int k = 0; k < ecenters[j].m; ++k) {
                    tc.edge_centers(j).Add(ecenters[j](k));
                }
            }

            //turn on segments on the faces
            for (int j = 0; j < 4; j++) {
                VECTOR<bool, 7> hit;
                for (int k = 0; k < 3; ++k) {
                    if (ecenters[face2edge[j][k]].m) {
                        hit(2*k+1) = 1;
                    }
                    if (vcenters[face2node[j][k]].m) {
                        hit(2*k) = 1;
                    }
                }
                if (fcenters[j].m) {
                    hit(6) = 1;
                    for (int k = 0; k < fcenters[j].m; ++k) {
                        tc.face_centers(j).Add(fcenters[j](k));
                    }
                }
                bool add = 0;
                for(int k = 0; k < 6; ++k){
                    if (hit(k) && (hit((k+3)%6) || hit(6) || ((k%2) && (hit((k+4)%6) || hit((k+2)%6))))){
                        tc.turned_on(j*6+k)=1;
                        add = 1;
                    }
                }
                if (add) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < ecenters[face2edge[j][k]].m; ++l) {
                            tc.face_centers(j).Add(ecenters[face2edge[j][k]](l));
                        }
                        for (int l = 0; l < vcenters[face2node[j][k]].m; ++l) {
                            tc.face_centers(j).Add(vcenters[face2node[j][k]](l));
                        }
                    }
                }
            }
            
            //add tet center
            if (tcenters.m) {
                for (int j = 0; j < tcenters.m; ++j) {
                    tc.tet_center.Add(tcenters(j));
                }
            }
            else {
                for (int j = 0; j < 4; ++j) {
                    for (int k = 0; k < vcenters[j].m; ++k) {
                        tc.tet_center.Add(vcenters[j](k));
                    }
                    for (int k = 0; k < fcenters[j].m; ++k) {
                        tc.tet_center.Add(fcenters[j](k));
                    }
                }
                for (int j = 0; j < 6; ++j) {
                    for (int k = 0; k < ecenters[j].m; ++k) {
                        tc.tet_center.Add(ecenters[j](k));
                    }
                }
            }
        }
    }
    
    //generate submaterials
    //set centers if it's on a turned-on face that is split(only one piece is in a connected component)
    VECTOR<bool, MaxPieces> picked;
    int NumMaterials = tc.material_ids.m;
    for (int i = 0; i < NumMaterials; i++) {
        if (!picked[tc.material_ids(i)]) {
            ARRAY<int> cc = tc.Find_CC(tc.material_ids(i), picked);
            //cout << cc << endl;
//            if (cc.m == NumMaterials) {
//                cout << tri2inter << endl;
//                cout << tc.turned_on << " not split " << NumMaterials << endl;
//            }
            if (1){//cc.m != NumMaterials || face_cut) {
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
                    if (i!=0) cout << "ridiculous!*************************************" << endl;
                    for (int j = 0; j < NumNodesPerTet; j++){
                        weights_in_sim.Append(PARENT(parent_sim_tet_id, wei[j]));
                        material_node_from.Append(tet_element(j));
                        cutting_particle_material_space.Append((TV)cutting_particle_material_space(tet_element(j)));
                    }
                    sim_volume->mesh.elements(parent_sim_tet_id) = VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1);
                    sim_tet_split(parent_sim_tet_id) = 1;
                }
                else {
                    ALGEBRA::MATRIX_3X3<ST> am = undeformed_config_copy(parent_sim_tet_id);
                    undeformed_config_copy.Append(am);
                    for (int j = 0; j < NumNodesPerTet; j++){
                        weights_in_sim.Append(PARENT(sim_volume->mesh.elements.m, wei[j]));
                        material_node_from.Append(tet_element(j));
                        cutting_particle_material_space.Append((TV)cutting_particle_material_space(tet_element(j)));
                    }
                    if (i == 0) {
                        ctet2stet(tet_id) = sim_volume->mesh.elements.m;
                    }
                    else {
                        ctet2stet.Append(sim_volume->mesh.elements.m);
                    }
                    sim_volume->mesh.elements.Append(VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1));
                    int stf = sim_tet_from(parent_sim_tet_id);
                    sim_tet_from.Append(stf);
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
                    bool ib = is_blue(tet_id);
                    is_blue.Append(ib);
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
void MESH_CUTTING<T>::Split2(const int& tet_id, HASHTABLE<int,H>& tri2inter, ARRAY<int>& sim_node_from, ARRAY<bool>& sim_tet_split, ARRAY<int>& material_node_from)
{
    //cout << "******" << tri2inter << endl;
    TET_CUTTING tc = tet_cuttings(tet_id);
    VECTOR<bool, NumSplitTurnOns> can_split;
    for (int i = 0; i < NumSplitTurnOns; i++){
        if(tc.has_material(face2material[i][0]) && tc.has_material(face2material[i][1])) {
            can_split(i) = 1;
        }}
    
    const VECTOR<int, NumNodesPerTet> tet_element = volume->mesh.elements(tet_id);
    for(typename HASHTABLE<int,H>::ITERATOR it(tri2inter);it.Valid();it.Next()) {
        H& intersects=it.Data();
        
        //centers will eventually be set as average of the related intersections, which works for cutting surface being approximately a plane inside the tet
        ARRAY<CENTER> vcenters[NumNodesPerTet];
        ARRAY<CENTER> ecenters[NumEdgesPerTet];
        ARRAY<CENTER> fcenters[NumFacesPerTet];
        ARRAY<CENTER> tcenters;
        int hits = 0;
        
        for (typename H::ITERATOR it2(intersects); it2.Valid(); it2.Next()) {
            I4 t = it2.Key();
            T4 w = it2.Data();
            if (t(1) == -1) {
                T4 c;
                int id = tet_element.Find(t(0));
                c[id] = 1;
                vcenters[id].Append(c);
                hits |= (1<<(1<<id));
            }
            else if (t(2) == -1) {
                T4 c;
                c[tet_element.Find(t(0))] = w(0);
                c[tet_element.Find(t(1))] = 1 - w(0);
                for (int i = 0; i < 6; ++i) {
                    if (I2(tet_element(edge2node[i][0]),tet_element(edge2node[i][1])).Sorted() == I2(t(0),t(1))) {
                        ecenters[i].Append(c);
                        int i1 = (1 << edge2node[i][0]);
                        int i2 = (1 << edge2node[i][1]);
                        hits |= (1 << (i1 | i2));
                        break;
                    }
                }
            }
            else if (t(3) == -1) {
                T4 c;
                for (int i = 0; i < 3; ++i) {
                    c(tet_element.Find(t(i))) = w(i);
                }
                for (int i = 0; i < 4; ++i) {
                    if (I3(tet_element(face2node[i][0]),tet_element(face2node[i][1]), tet_element(face2node[i][2])).Sorted() == t.Remove_Index(3)) {
                        fcenters[i].Append(c);
                        int i1 = (1 << face2node[i][0]);
                        int i2 = (1 << face2node[i][1]);
                        int i3 = (1 << face2node[i][2]);
                        hits |= (1<<(i1+i2+i3));
                        break;
                    }
                }
            }
            else {
                T4 c;
                for (int i = 0; i < 4; ++i) {
                    c(tet_element.Find(t(i))) = w(i);
                }
                tcenters.Append(c);
                hits |= (1<<15);
            }
        }
        
        for (int j = 0; j < 4; ++j) {
            int e1 = face2edge[j][0];
            int e2 = face2edge[j][1];
            int e3 = face2edge[j][2];
            int n1 = face2node[j][0];
            int n2 = face2node[j][1];
            int n3 = face2node[j][2];
            bool cf = 0;
            if (vcenters[n1].m && (vcenters[n2].m || ecenters[e1].m) && (vcenters[n3].m || ecenters[e2].m || ecenters[e3].m || fcenters[j].m)) {
                tc.no_merge(j*6) = 1;
                cf = 1;
            }
            if (vcenters[n2].m && (vcenters[n3].m || ecenters[e2].m) && (vcenters[n1].m || ecenters[e1].m || ecenters[e3].m || fcenters[j].m)) {
                tc.no_merge(j*6+2) = 1;
                cf = 1;
            }
            if (vcenters[n3].m && (vcenters[n1].m || ecenters[e3].m) && (vcenters[n2].m || ecenters[e1].m || ecenters[e2].m || fcenters[j].m)) {
                tc.no_merge(j*6+4) = 1;
                cf = 1;
            }
            if (vcenters[n2].m && (vcenters[n1].m || ecenters[e1].m) && (vcenters[n3].m || ecenters[e2].m || ecenters[e3].m || fcenters[j].m)) {
                tc.no_merge(j*6+1) = 1;
                cf = 1;
            }
            if (vcenters[n3].m && (vcenters[n2].m || ecenters[e2].m) && (vcenters[n1].m || ecenters[e1].m || ecenters[e3].m || fcenters[j].m)) {
                tc.no_merge(j*6+3) = 1;
                cf = 1;
            }
            if (vcenters[n1].m && (vcenters[n3].m || ecenters[e3].m) && (vcenters[n2].m || ecenters[e1].m || ecenters[e2].m || fcenters[j].m)) {
                tc.no_merge(j*6+5) = 1;
                cf = 1;
            }
            if (cf) {
                cutFaces.Set(I3(tet_element(n1), tet_element(n2), tet_element(n3)).Sorted());
            }
        }

        for (int j = 0; j < NumEdgesPerTet; j++){
            int n1 = edge2node[j][0];
            int n2 = edge2node[j][1];
            int n3 = 15 ^ ((1<<n1) | (1<<n2));
            int n4 = 15 ^ (1<<n1);
            int n5 = 15 ^ (1<<n2);
            int m = (1<<n3) | (1<<n4) | (1<<n5) | (1<<15);
            if (vcenters[n1].m && (ecenters[j].m || vcenters[n2].m) && (hits & m)) {
                int id = 2 * j + 24;
                tc.turned_on[id] = 1;
            }
            if (vcenters[n2].m && (ecenters[j].m || vcenters[n1].m) && (hits & m)) {
                int id = 2 * j + 1 + 24;
                tc.turned_on[id] = 1;
            }
            for (int k = 0; k < ecenters[j].m; ++k) {
                tc.edge_centers(j).Add(ecenters[j](k));
            }
        }
        
        //turn on segments on the faces
        for (int j = 0; j < 4; j++) {
            int e1 = face2edge[j][0];
            int e2 = face2edge[j][1];
            int e3 = face2edge[j][2];
            int n1 = face2node[j][0];
            int n2 = face2node[j][1];
            int n3 = face2node[j][2];
            if ((vcenters[n1].m+vcenters[n2].m+vcenters[n3].m+ecenters[e1].m+ecenters[e2].m+ecenters[e3].m+fcenters[j].m) != intersects.Size()) {
                VECTOR<bool, 7> hit;
                for (int k = 0; k < 3; ++k) {
                    if (ecenters[face2edge[j][k]].m) {
                        hit(2*k+1) = 1;
                    }
                    if (vcenters[face2node[j][k]].m) {
                        hit(2*k) = 1;
                    }
                }
                if (fcenters[j].m) {
                    hit(6) = 1;
                    for (int k = 0; k < fcenters[j].m; ++k) {
                        tc.face_centers(j).Add(fcenters[j](k));
                    }
                }
                bool add = 0;
                for(int k = 0; k < 6; ++k){
                    if (hit(k) && (hit((k+3)%6) || hit(6) || ((k%2) && (hit((k+4)%6) || hit((k+2)%6))))){
                        tc.turned_on(j*6+k)=1;
                        add = 1;
                    }
                }
                if (add) {
                    for (int k = 0; k < 3; ++k) {
                        for (int l = 0; l < ecenters[face2edge[j][k]].m; ++l) {
                            tc.face_centers(j).Add(ecenters[face2edge[j][k]](l));
                        }
                        for (int l = 0; l < vcenters[face2node[j][k]].m; ++l) {
                            tc.face_centers(j).Add(vcenters[face2node[j][k]](l));
                        }
                    }
                }
            }
        }
        
        //add tet center
        if (tcenters.m) {
            for (int j = 0; j < tcenters.m; ++j) {
                tc.tet_center.Add(tcenters(j));
            }
        }
        else {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < vcenters[j].m; ++k) {
                    tc.tet_center.Add(vcenters[j](k));
                }
                for (int k = 0; k < fcenters[j].m; ++k) {
                    tc.tet_center.Add(fcenters[j](k));
                }
            }
            for (int j = 0; j < 6; ++j) {
                for (int k = 0; k < ecenters[j].m; ++k) {
                    tc.tet_center.Add(ecenters[j](k));
                }
            }
        }
    }
    
    //generate submaterials
    //set centers if it's on a turned-on face that is split(only one piece is in a connected component)
    VECTOR<bool, MaxPieces> picked;
    int NumMaterials = tc.material_ids.m;
    for (int i = 0; i < NumMaterials; i++) {
        if (!picked[tc.material_ids(i)]) {
            ARRAY<int> cc = tc.Find_CC(tc.material_ids(i), picked);
            //cout << cc << endl;
//            if (cc.m == NumMaterials) {
//                cout << endl << "hit but not split tet " << tet_id << ": " << volume->mesh.elements(tet_id) << "    " << volume->particles.X.Subset(volume->mesh.elements(tet_id)) << endl;
//                cout <<  tri2inter << endl;
//                cout << "turned on: " << tc.turned_on << endl;
//            }
            if (1){//cc.m != NumMaterials || face_cut) {
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
                    if (i!=0) cout << "ridiculous!*************************************" << endl;
                    for (int j = 0; j < NumNodesPerTet; j++){
                        weights_in_sim.Append(PARENT(parent_sim_tet_id, wei[j]));
                        material_node_from.Append(tet_element(j));
                        cutting_particle_material_space.Append((TV)cutting_particle_material_space(tet_element(j)));
                    }
                    sim_volume->mesh.elements(parent_sim_tet_id) = VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1);
                    sim_tet_split(parent_sim_tet_id) = 1;
                }
                else {
                    ALGEBRA::MATRIX_3X3<ST> am = undeformed_config_copy(parent_sim_tet_id);
                    undeformed_config_copy.Append(am);
                    for (int j = 0; j < NumNodesPerTet; j++){
                        weights_in_sim.Append(PARENT(sim_volume->mesh.elements.m, wei[j]));
                        material_node_from.Append(tet_element(j));
                        cutting_particle_material_space.Append((TV)cutting_particle_material_space(tet_element(j)));
                    }
                    if (i == 0) {
                        ctet2stet(tet_id) = sim_volume->mesh.elements.m;
                    }
                    else {
                        ctet2stet.Append(sim_volume->mesh.elements.m);
                    }
                    sim_volume->mesh.elements.Append(VECTOR<int, NumNodesPerTet>(sim_size-4, sim_size-3, sim_size-2, sim_size-1));
                    int stf = sim_tet_from(parent_sim_tet_id);
                    sim_tet_from.Append(stf);
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
                    bool ib = is_blue(tet_id);
                    is_blue.Append(ib);
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
int MESH_CUTTING<T>::Sorted_Id(const I3& sorted_tri, const I3& tri, int material_id) {
    I3 in(sorted_tri.Find(tri(0)), sorted_tri.Find(tri(1)), sorted_tri.Find(tri(2)));
    bool inv = (in == I3(0, 2, 1) || in == I3(1, 0, 2) || in == I3(2, 1, 0));
    return tri_node2material[sorted_tri.Find(tri(((material_id+1)%6)/2))][inv^(material_id&1)];
}

template<class T>
void MESH_CUTTING<T>::Cut(TRIANGULATED_SURFACE<T>& cutting_surface, bool refine, bool material_space)
{
    cout << "cutting*************************************" << endl;
    
    start_timer();
    
    //preprocessing
    cutFaces.Clean_Memory();
    
    //copy the volume to get rid of old hierarchy...
    TETRAHEDRALIZED_VOLUME<T> *tv = TETRAHEDRALIZED_VOLUME<T>::Create();
    tv->mesh.elements = volume->mesh.elements;
    tv->particles.Add_Elements(volume->particles.X.m);
    tv->Update_Number_Nodes();
    for (int i = 0; i < volume->particles.X.m; ++i) {
        tv->particles.X(i) = volume->particles.X(i);
    }
    delete volume;
    volume = tv;
    
    //do this if doing incremental cutting in material space
    if (material_space) {
        for(int i=0;i<volume->particles.X.m;++i)
            volume->particles.X(i)=cutting_particle_material_space(i);
    }

    //keep a copy of original element-wise data
    original_elements = volume->mesh.elements;
    original_sim_elements = sim_volume->mesh.elements;
    original_ctet2stet = ctet2stet;
    
    //intersections
    volume->Update_Number_Nodes();
    CONSISTENT_INTERSECTIONS<TV> intersections(*volume,cutting_surface);
    intersections.Compute();
    
    HASHTABLE<I4,int> tet_from_tet;
    HASHTABLE<I3,ARRAY<int>> tet_from_face;
    HASHTABLE<I2,ARRAY<int> > tet_from_edge;
    ARRAY<ARRAY<int> > tet_from_vertex(volume->particles.number);
    
    for(int i=0;i<volume->mesh.elements.m;i++){
        I4 e=volume->mesh.elements(i).Sorted();
        tet_from_tet.Set(e,i);
        for(int j=0;j<4;j++) tet_from_face.Get_Or_Insert(e.Remove_Index(j)).Append(i);
        for(int j=0;j<3;j++) {
            for (int k = j+1; k < 4; ++k) {
                tet_from_edge.Get_Or_Insert(I2(e(j),e(k))).Append(i);
            }
        }
        for(int j=0;j<4;j++) tet_from_vertex(e(j)).Append(i);
    }
    
    HASHTABLE<I3,int> tri_from_face;
    HASHTABLE<I2,ARRAY<int>> tri_from_edge;
    ARRAY<ARRAY<int> > tri_from_vertex(cutting_surface.particles.number);
    
    for(int i=0;i<cutting_surface.mesh.elements.m;i++){
        I3 e=cutting_surface.mesh.elements(i).Sorted();
        tri_from_face.Set(e,i);
        for (int j = 0; j < 3; ++j) {
            tri_from_edge.Get_Or_Insert(e.Remove_Index(j)).Append(i);
        }
        for(int j=0;j<3;j++) tri_from_vertex(e(j)).Append(i);
    }

    HASHTABLE<int,HASHTABLE<int,H> > components;
    
    for(typename HASHTABLE<I2>::ITERATOR it(intersections.hash_vv);it.Valid();it.Next()){
        I2 key=it.Key();
        const ARRAY<int>& a=tet_from_vertex(key.x);
        const ARRAY<int>& b=tri_from_vertex(key.y);
        for(int i=0;i<a.m;i++)
            for(int j=0;j<b.m;j++)
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Get_Or_Insert(I4(key(0),-1,-1,-1))=T4();
    }
    
    for(typename HASHTABLE<I3,T>::ITERATOR it(intersections.hash_ve);it.Valid();it.Next()){
        I3 key=it.Key();
        const ARRAY<int>& a=tet_from_vertex(key.x);
        const ARRAY<int>& b=tri_from_edge.Get(key.Remove_Index(0));
        for(int i=0;i<a.m;i++){
            for(int j=0;j<b.m;++j) {
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Get_Or_Insert(I4(key(0),-1,-1,-1))=T4();
            }
        }
    }
    
    for(typename HASHTABLE<I3,T>::ITERATOR it(intersections.hash_ev);it.Valid();it.Next()){
        I3 key=it.Key();
        const ARRAY<int>& a=tet_from_edge.Get(key.Remove_Index(2));
        const ARRAY<int>& b=tri_from_vertex(key.z);
        for(int i=0;i<a.m;i++)
            for(int j=0;j<b.m;j++)
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Get_Or_Insert(I4(key(0),key(1),-1,-1))=T4(1-it.Data(),0,0,0);
    }
    
    for(typename HASHTABLE<I4,T2>::ITERATOR it(intersections.hash_ee);it.Valid();it.Next()){
        I4 key=it.Key();
        const ARRAY<int>& a=tet_from_edge.Get(I2(key(0),key(1)));
        const ARRAY<int>& b=tri_from_edge.Get(I2(key(2),key(3)));
        for(int i=0;i<a.m;i++)
            for(int j=0;j<b.m;j++)
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Get_Or_Insert(I4(key(0),key(1),-1,-1))=T4(1-it.Data()(0),0,0,0);
    }
    
    for(typename HASHTABLE<I4,TV>::ITERATOR it(intersections.hash_fv);it.Valid();it.Next()){
        I4 key=it.Key();
        const ARRAY<int>& a=tet_from_face.Get(key.Remove_Index(3));
        const ARRAY<int>& b=tri_from_vertex(key.Last());
        key(3)=-1;
        for(int i=0;i<a.m;i++)
            for(int j=0;j<b.m;j++)
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Get_Or_Insert(key)=it.Data().Append(0);
    }
    
    for(typename HASHTABLE<I4,TV>::ITERATOR it(intersections.hash_vf);it.Valid();it.Next()){
        I4 key=it.Key();
        const ARRAY<int>& a=tet_from_vertex(key(0));
        int b=tri_from_face.Get(key.Remove_Index(0));
        for(int i=0;i<a.m;i++)
            components.Get_Or_Insert(a(i)).Get_Or_Insert(b).Get_Or_Insert(I4(key(0),-1,-1,-1))=T4();
    }
    
    for(typename HASHTABLE<I5,T4>::ITERATOR it(intersections.hash_ef);it.Valid();it.Next()){
        I5 key=it.Key();
        const ARRAY<int>& a=tet_from_edge.Get(I2(key(0),key(1)));
        int b=tri_from_face.Get(I3(key(2),key(3),key(4)));
        for(int i=0;i<a.m;i++){
            components.Get_Or_Insert(a(i)).Get_Or_Insert(b).Get_Or_Insert(I4(key(0),key(1),-1,-1))=T4(1-it.Data()(0),0,0,0);
        }
    }
    
    for(typename HASHTABLE<I5,T4>::ITERATOR it(intersections.hash_fe);it.Valid();it.Next()){
        I5 key=it.Key();
        const ARRAY<int>& a=tet_from_face.Get(I3(key(0),key(1),key(2)));
        const ARRAY<int>& b=tri_from_edge.Get(I2(key(3),key(4)));
        it.Data()(3)=0;
        for(int i=0;i<a.m;i++)
            for(int j=0;j<b.m;j++)
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Get_Or_Insert(I4(key(0),key(1),key(2),-1))=it.Data();
    }
    
    for(typename HASHTABLE<I5,T4>::ITERATOR it(intersections.hash_tv);it.Valid();it.Next()){
        I5 key=it.Key();
        int a=tet_from_tet.Get(key.Remove_Index(4));
        const ARRAY<int>& b=tri_from_vertex(key.Last());
        for(int j=0;j<b.m;j++)
            components.Get_Or_Insert(a).Get_Or_Insert(b(j)).Get_Or_Insert(key.Remove_Index(4))=it.Data();
    }
    //cout << components << endl;
    
    //split
    ARRAY<TET>& elements = volume->mesh.elements;
    ARRAY<TET>& sim_elements = sim_volume->mesh.elements;
    
    int num_elements = elements.m; cout << num_elements << " elements" << endl;
    int num_nodes = volume->particles.number;
    ARRAY<int> material_node_from(num_nodes);
    for (int i = 0; i < num_nodes; ++i) {
        material_node_from(i) = i;
    }
    ARRAY<bool> need_dup(num_nodes);//whether the node needs to be duplicated
    ARRAY<bool> need_merge(num_elements);//whether the element needs to be merged

    ARRAY<int> sim_node_from(sim_volume->particles.X.m);
    ARRAY<bool> sim_tet_split(sim_volume->mesh.elements.m);
    
    for (int i = 0; i < sim_volume->particles.X.m; i++){
        sim_node_from(i) = i;
    }
    for (int i = 0; i < sim_volume->mesh.elements.m; i++){
        sim_tet_split(i) = 0;
        sim_tet_from.Append(i);
    }
    
    //split: data to take care of: sim_node_from, weights_in_sim, sim_volume->mesh.elements(always completely newly created, because each child tet(split or not) gets a new parent tet), volume->mesh.elements, ctet2stet
    for(typename HASHTABLE<int,HASHTABLE<int,H> >::ITERATOR it(components);it.Valid();it.Next()){
        int i=it.Key();
        Split(i, it.Data(), sim_node_from, sim_tet_split, material_node_from);
        need_merge(i) = 1;
        need_dup.Subset(original_elements(i)).Fill(1);
    }
    cout << components.Size() << " tets touched\n";

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
                    PHYSBAM_ASSERT(cutting_particle_material_space.m == weights_in_sim.m);
                    if (!sim_tet_split(parent_sim_tet_id)){   
                        weights_in_sim.Append(PARENT(parent_sim_tet_id, w));
                        material_node_from.Append(original_elements(i)(j));
                        cutting_particle_material_space.Append((TV)cutting_particle_material_space(tet_element(j)));
                    }
                    else {
                        weights_in_sim.Append(PARENT(sim_volume->mesh.elements.m, w));
                        material_node_from.Append(original_elements(i)(j));
                        cutting_particle_material_space.Append((TV)cutting_particle_material_space(tet_element(j)));
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
                    ALGEBRA::MATRIX_3X3<ST> ucc = undeformed_config_copy(parent_sim_tet_id);
                    undeformed_config_copy.Append(ucc);
                }
            }
        }
    }
    
    stop_timer();
    printf("split time:    %f\n",get_time());
    cout << "total particles: " << weights_in_sim.m << " " << material_node_from.m << endl;
//    cout << elements << endl;
//    cout << material_node_from << endl;
    
    //merge nodes
    //identify which un-plit tet need merging
    start_timer();
    
    HASHTABLE<I4,int> ht;
    UNION_FIND<int> node_classes(weights_in_sim.m);//cluster merging nodes
    UNION_FIND<int> sim_node_classes(sim_node_from.m);
    cout << "cut faces: " << cutFaces << endl;

    for (int i = 0; i < elements.m; i++) {
        //cout << tet_cuttings(i).material_ids << endl;
        if (1){//need_merge(i)) {
            I4 pt = sim_elements(ctet2stet(i));
            for (int j = 0; j < 4; ++j) {
                I3 f(elements(i)(face2node[j][0]), elements(i)(face2node[j][1]), elements(i)(face2node[j][2]));
                I3 of(material_node_from(f(0)), material_node_from(f(1)), material_node_from(f(2)));
                I3 sof = of.Sorted();
                for (int k = 0; k < 6; ++k) {
                    int material_id = j * 6 + k;
                    if (!tet_cuttings(i).no_merge(material_id)) {
                        int sorted_k = Sorted_Id(sof, of, k);
                        if (tet_cuttings(i).has_material(material_id)) {
                            int tid = i;
                            //cout << sof.Append(sorted_k) << endl;
                            if (ht.Get(sof.Append(sorted_k),tid)) {
                                for (int l = 0; l < 3; ++l) {
                                    for (int m = 0; m < 4; ++m) {
                                        if (of(l) == material_node_from(elements(tid)(m))) {
                                            node_classes.Union(f(l), elements(tid)(m));
                                            //cout << "union " << f(l) << " " << elements(tid)(m) << endl;
                                            break;
                                        }
                                    }
                                }
                                I4 pt2 = sim_elements(ctet2stet(tid));
                                for (int l = 0; l < 4; ++l) {
                                    for (int m = 0; m < 4; ++m) {
                                        if (sim_node_from(pt(l)) == sim_node_from(pt2(m))) {
                                            sim_node_classes.Union(pt(l), pt2(m));
                                            break;
                                        }
                                    }
                                }
                            }
                            else {
                                ht.Set(sof.Append(sorted_k),i);
                            }
                        }
                    }
                    else {
                        cout << "no merge: " << i << " " << material_id << endl;
                    }
                }
            }
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
    ARRAY<ALGEBRA::MATRIX_3X3<ST> > new_undeformed_config_copy;
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
    dragging_elements = new_dragging_elements;
    dragging_weights = new_dragging_weights;
    dragging_targets = new_dragging_targets;
    sim_volume->mesh.elements = new_sim_elements;
    undeformed_config_copy = new_undeformed_config_copy;
    sim_volume->particles.Resize(new_sim_particles.m);
    for (int i = 0; i < new_sim_particles.m; i++){
        sim_volume->particles.X(i) = new_sim_particles(i);
    }
//    cout << volume->mesh.elements << endl;
//    cout << sim_volume->mesh.elements << endl;
    
    stop_timer();
    printf("merge time:    %f\n",get_time());  
    
    original_ctet2stet.Remove_All();
    original_elements.Remove_All();
    original_sim_elements.Remove_All();
    sim_tet_from.Remove_All();
    Update_Cutting_Particles();
    
    if (refine) {
        Partial_Refine();
    }
    
    //reinitialize boundary mesh
    sim_volume->Update_Number_Nodes();
    sim_volume->mesh.Initialize_Boundary_Mesh(); //cout << "sim boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
    sim_volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
    
    Fix_Orientation();
    volume->Update_Number_Nodes();
    volume->mesh.Initialize_Boundary_Mesh(); //cout << "cutting boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
    volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
    
    //reinitialize elasticity
    if (!interactive) {
        start_timer();
        Reinitialize_Elasticity();
        stop_timer();
        printf("elasticity reinitialization time:    %f\n",get_time());
    }
    
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
    int n = volume->mesh.elements.m;
    
    for (int i = 0; i < n; i++) {
        if (tet_cuttings(i).material_ids.m < MaxPieces){
            const int parent_sim_tet = ctet2stet(i);
            VECTOR<int,4> &element = volume->mesh.elements(i);
            TET_CUTTING& tc = tet_cuttings(i);
            ARRAY<I4> ne;
            
            //tet center
            int cid = weights_in_sim.m;
            {
                CENTER c = Weight_In_Sim_Tet(tc.tet_center.Value(), element, parent_sim_tet);
                weights_in_sim.Append(PARENT(parent_sim_tet, c));
                cutting_particle_material_space.Append(weight2vec_material_space(i, c));
            }
            
            for (int j = 0; j < 4; ++j) {
                int fc = weights_in_sim.m;
                int n1 = element(face2node[j][0]);
                int n2 = element(face2node[j][1]);
                int n3 = element(face2node[j][2]);
                bool add = 0;
                for (int k = 0; k < 6; ++k) {
                    if (!tet_cuttings(i).has_material(j*6+k)) {
                        CENTER c = Weight_In_Sim_Tet(tc.face_centers[j].Value(), element, parent_sim_tet);
                        add = 1;
                        I3 face = I3(n1, n2, n3).Sorted();
                        if (!new_nodes2.Get(face, fc)) {
                            new_nodes2.Set(face, fc);
                            weights_in_sim.Append(PARENT(parent_sim_tet, c));
                            cutting_particle_material_space.Append(weight2vec_material_space(i, c));
                        }
                        break;
                    }
                }
                if (add) {
                    for (int k = 0; k < 3; ++k) {
                        int fn1 = element(face2node[j][k]);
                        int fn2 = element(face2node[j][(k+1)%3]);
                        int m1 = j * 6 + k * 2;
                        int m2 = j * 6 + k * 2 + 1;
                        if (tet_cuttings(i).has_material(m1) && tet_cuttings(i).has_material(m2)) {
                            ne.Append(I4(cid, fc, fn1, fn2));
                        }
                        else if (tet_cuttings(i).has_material(m1) || tet_cuttings(i).has_material(m2)){
                            int eid = weights_in_sim.m;
                            I2 edge = I2(fn1, fn2).Sorted();
                            if (!new_nodes.Get(edge, eid)) {
                                CENTER c = Weight_In_Sim_Tet(tc.edge_centers[face2edge[j][k]].Value(), element, parent_sim_tet);
                                new_nodes.Set(edge,eid);
                                weights_in_sim.Append(PARENT(parent_sim_tet, c));
                                cutting_particle_material_space.Append(weight2vec_material_space(i, c));
                            }
                            if (tet_cuttings(i).has_material(m1)) {
                                ne.Append(I4(cid, fc, fn1, eid));
                            }
                            else {
                                ne.Append(I4(cid, fc, fn2, eid));
                            }
                        }
                    };
                }
                else {
                    ne.Append(I4(cid, n1, n2, n3));
                }
            }
            element = ne(0);
            for (int j = 1; j < ne.m; ++j) {
                volume->mesh.elements.Append(ne(j));
                bool ib = is_blue(i);
                is_blue.Append(ib);
                int c2s = ctet2stet(i);
                ctet2stet.Append(c2s);
            }
        }
    }
    
    //subdivide by face intersections
    for (int i = 0; i < volume->mesh.elements.m; i++) {
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
                    if (!ht.Get(face, k)) {
                        cout << "should have been set!!!\n";
                    }
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
            int nc2s = new_ctet2stet(tid);
            new_ctet2stet.Append(nc2s);
            int nib = new_is_blue(tid);
            new_is_blue.Append(nib);
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

    //delete unused nodes in volume
    HASHTABLE<int,int> new_pids;
    int new_pid=0;
    ARRAY<PARENT> new_weights_in_sim;
    ARRAY<TV> new_cutting_particles_material_space;
    for(int i=0;i<new_elements.m;++i)
        for(int j=0;j<4;++j){
            int& id=new_elements(i)(j);
            if(!new_pids.Get(id,id)){
                new_weights_in_sim.Append(weights_in_sim(id));
                new_cutting_particles_material_space.Append(cutting_particle_material_space(id));
                new_pids.Set(id,new_pid);
                id=new_pid;
                ++new_pid;}}
    weights_in_sim=new_weights_in_sim;

    //reset data
    volume->mesh.elements = new_elements;
    tet_cuttings.Remove_All();
    tet_cuttings.Resize(volume->mesh.elements.m);
    ctet2stet = new_ctet2stet;
    is_blue = new_is_blue;
    volume->Update_Number_Nodes();
    Update_Cutting_Particles();
    
    cutting_particle_material_space = new_cutting_particles_material_space;
//    cutting_particle_material_space.Resize(volume->particles.X.m);
//    for(int i=0;i<volume->particles.X.m;++i)
//        cutting_particle_material_space(i)=volume->particles.X(i);
}

template<class T>
void MESH_CUTTING<T>::Refine_And_Save_To(TETRAHEDRALIZED_VOLUME<T>* refined_volume)
{
    refined_volume->mesh.elements = volume->mesh.elements;
    
    HASHTABLE<VECTOR<int,2>,int> new_nodes;//edge to edge intersection point's index
    HASHTABLE<VECTOR<int,2>, std::set<int> > edge2tets;
    HASHTABLE<VECTOR<int,3>,int> new_nodes2;
    ARRAY<TET> new_elements;
    ARRAY<PARENT> new_weights_in_sim = weights_in_sim;
    
    int n = refined_volume->mesh.elements.m;
    for (int i = 0; i < n; i++) {
        if (tet_cuttings(i).material_ids.m < MaxPieces){
            const int parent_sim_tet = ctet2stet(i);
            VECTOR<int,4> &element = refined_volume->mesh.elements(i);
            TET_CUTTING& tc = tet_cuttings(i);
            ARRAY<I4> ne;
            
            //tet center
            int cid = new_weights_in_sim.m;
            {
                CENTER c = Weight_In_Sim_Tet(tc.tet_center.Value(), element, parent_sim_tet);
                new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
            }
            
            for (int j = 0; j < 4; ++j) {
                int fc = new_weights_in_sim.m;
                int n1 = element(face2node[j][0]);
                int n2 = element(face2node[j][1]);
                int n3 = element(face2node[j][2]);
                bool add = 0;
                for (int k = 0; k < 6; ++k) {
                    if (!tet_cuttings(i).has_material(j*6+k)) {
                        CENTER c = Weight_In_Sim_Tet(tc.face_centers[j].Value(), element, parent_sim_tet);
                        add = 1;
                        I3 face = I3(n1, n2, n3).Sorted();
                        if (!new_nodes2.Get(face, fc)) {
                            new_nodes2.Set(face, fc);
                            new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
                        }
                        break;
                    }
                }
                if (add) {
                    for (int k = 0; k < 3; ++k) {
                        int fn1 = element(face2node[j][k]);
                        int fn2 = element(face2node[j][(k+1)%3]);
                        int m1 = j * 6 + k * 2;
                        int m2 = j * 6 + k * 2 + 1;
                        if (tet_cuttings(i).has_material(m1) && tet_cuttings(i).has_material(m2)) {
                            ne.Append(I4(cid, fc, fn1, fn2));
                        }
                        else if (tet_cuttings(i).has_material(m1) || tet_cuttings(i).has_material(m2)){
                            int eid = new_weights_in_sim.m;
                            I2 edge = I2(fn1, fn2).Sorted();
                            if (!new_nodes.Get(edge, eid)) {
                                CENTER c = Weight_In_Sim_Tet(tc.edge_centers[face2edge[j][k]].Value(), element, parent_sim_tet);
                                new_nodes.Set(edge,eid);
                                new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
                            }
                            if (tet_cuttings(i).has_material(m1)) {
                                ne.Append(I4(cid, fc, fn1, eid));
                            }
                            else {
                                ne.Append(I4(cid, fc, fn2, eid));
                            }
                        }
                    };
                }
                else {
                    ne.Append(I4(cid, n1, n2, n3));
                }
            }
            element = ne(0);
            for (int j = 1; j < ne.m; ++j) {
                refined_volume->mesh.elements.Append(ne(j));
            }
        }
    }
    
    //subdivide by face intersections
    for (int i = 0; i < refined_volume->mesh.elements.m; i++) {
        TET element = refined_volume->mesh.elements(i);
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
                    if (!ht.Get(face, k)) {
                        cout << "should have been set!!!\n";
                    }
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
    
    refined_volume->mesh.elements = new_elements;
    refined_volume->particles.Resize(new_weights_in_sim.m);
    for (int i = 0; i < new_weights_in_sim.m; i++){
        refined_volume->particles.X(i) = weight2vec_sim(new_weights_in_sim(i).id, new_weights_in_sim(i).weight);
    }
    refined_volume->Update_Number_Nodes();
}

template<class T>
void MESH_CUTTING<T>::Refine_And_Save_To(TETRAHEDRALIZED_VOLUME<T>* refined_volume, HASHTABLE<I3>& cutting_faces, HASHTABLE<I3>& new_cutting_faces)
{
    refined_volume->mesh.elements = volume->mesh.elements;
    
    HASHTABLE<VECTOR<int,2>,int> new_nodes;//edge to edge intersection point's index
    HASHTABLE<VECTOR<int,2>, std::set<int> > edge2tets;
    HASHTABLE<VECTOR<int,3>,int> new_nodes2;
    ARRAY<TET> new_elements;
    ARRAY<PARENT> new_weights_in_sim = weights_in_sim;
    new_cutting_faces.Remove_All();
    
    int n = refined_volume->mesh.elements.m;
    for (int i = 0; i < n; i++) {
        if (tet_cuttings(i).material_ids.m < MaxPieces){
            const int parent_sim_tet = ctet2stet(i);
            VECTOR<int,4> &element = refined_volume->mesh.elements(i);
            TET_CUTTING& tc = tet_cuttings(i);
            ARRAY<I4> ne;
            
            //tet center
            int cid = new_weights_in_sim.m;
            {
                CENTER c = Weight_In_Sim_Tet(tc.tet_center.Value(), element, parent_sim_tet);
                new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
            }
            
            for (int j = 0; j < 4; ++j) {
                int fc = new_weights_in_sim.m;
                int n1 = element(face2node[j][0]);
                int n2 = element(face2node[j][1]);
                int n3 = element(face2node[j][2]);
                bool add = 0;
                for (int k = 0; k < 6; ++k) {
                    if (!tet_cuttings(i).has_material(j*6+k)) {
                        CENTER c = Weight_In_Sim_Tet(tc.face_centers[j].Value(), element, parent_sim_tet);
                        add = 1;
                        I3 face = I3(n1, n2, n3).Sorted();
                        if (!new_nodes2.Get(face, fc)) {
                            new_nodes2.Set(face, fc);
                            new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
                        }
                        break;
                    }
                }
                if (add) {
                    for (int k = 0; k < 3; ++k) {
                        int fn1 = element(face2node[j][k]);
                        int fn2 = element(face2node[j][(k+1)%3]);
                        int m1 = j * 6 + k * 2;
                        int m2 = j * 6 + k * 2 + 1;
                        if (tet_cuttings(i).has_material(m1) && tet_cuttings(i).has_material(m2)) {
                            ne.Append(I4(cid, fc, fn1, fn2));
                        }
                        else if (tet_cuttings(i).has_material(m1) || tet_cuttings(i).has_material(m2)){
                            int eid = new_weights_in_sim.m;
                            I2 edge = I2(fn1, fn2).Sorted();
                            if (!new_nodes.Get(edge, eid)) {
                                CENTER c = Weight_In_Sim_Tet(tc.edge_centers[face2edge[j][k]].Value(), element, parent_sim_tet);
                                new_nodes.Set(edge,eid);
                                new_weights_in_sim.Append(PARENT(parent_sim_tet, c));
                            }
                            if (tet_cuttings(i).has_material(m1)) {
                                ne.Append(I4(cid, fc, fn1, eid));
                            }
                            else {
                                ne.Append(I4(cid, fc, fn2, eid));
                            }
                        }
                    };
                }
                else {
                    ne.Append(I4(cid, n1, n2, n3));
                }
            }
            element = ne(0);
            for (int j = 1; j < ne.m; ++j) {
                refined_volume->mesh.elements.Append(ne(j));
            }
        }
    }
    
    //subdivide by face intersections
    for (int i = 0; i < refined_volume->mesh.elements.m; i++) {
        TET element = refined_volume->mesh.elements(i);
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
                    if (!ht.Get(face, k)) {
                        cout << "should have been set!!!\n";
                    }
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
    
    //delete unused nodes in volume
    HASHTABLE<int,int> new_pids;
    int new_pid=0;
    ARRAY<PARENT> new_weights_in_sim2;
    for(int i=0;i<new_elements.m;++i)
        for(int j=0;j<4;++j){
            int& id=new_elements(i)(j);
            if(!new_pids.Get(id,id)){
                new_weights_in_sim2.Append(new_weights_in_sim(id));
                new_pids.Set(id,new_pid);
                id=new_pid;
                ++new_pid;}}
    
    //reset indices in cutting_faces
//    for (HASHTABLE_ITERATOR<I3> it(new_cutting_faces); it.Valid(); it.Next()) {
//        I3& f = it.Key();
//        for (int i = 0; i < 3; ++i) {
//            new_pids.Get(f(i), f(i));
//        }
//    }
    
    refined_volume->mesh.elements = new_elements;
    refined_volume->particles.Resize(new_weights_in_sim2.m);
    for (int i = 0; i < new_weights_in_sim2.m; i++){
        refined_volume->particles.X(i) = weight2vec_sim(new_weights_in_sim2(i).id, new_weights_in_sim2(i).weight);
    }
    refined_volume->Update_Number_Nodes();
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

//for virtual surgery dragging
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
        color_vec.Append(TV((T)element2cc(ctet2stet(f[0]))/tet_cc.m,0,0));
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
                mesh_vec.Append(weight2vec(f[0], tet_cuttings(f[0]).edge_centers[edge_id].Value()));
                
                mesh_vec.Append(weight2vec(f[0], tet_cuttings(f[0]).face_centers[f[1]].Value()));
                
                refine_edges.Append(volume->particles.X(volume->mesh.elements(f[0])(node_id))*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).edge_centers[edge_id].Value())*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).edge_centers[edge_id].Value())*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).face_centers[f[1]].Value())*MAG);
                refine_edges.Append(weight2vec(f[0], tet_cuttings(f[0]).face_centers[f[1]].Value())*MAG);
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
    TV v1 = weight2vec(f[0], tet_cuttings(f[0]).tet_center.Value());
    mesh_vec.Append(v1);
    refine_edges.Append(v1);
    
    if (f[1] < 24) {
        int face_id = f[1]/6;
        TV v2 = weight2vec(f[0], tet_cuttings(f[0]).face_centers[face_id].Value());
        mesh_vec.Append(v2);
        refine_edges.Append(v2);
        refine_edges.Append(v2);
        
        int id_within_face = (f[1]%6)/2;
        if ((f[1]%6) & 1){
            TV v3 = weight2vec(f[0], tet_cuttings(f[0]).edge_centers[face2edge[face_id][id_within_face]].Value());
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

template<class T>
typename MESH_CUTTING<T>::TV MESH_CUTTING<T>::weight2vec_material_space(int tet_id, CENTER c) const
{
    TV v;
    for (int i = 0; i < NumNodesPerTet; i++) {
        v += cutting_particle_material_space(volume->mesh.elements(tet_id)(i))*c[i];
    }
    return v;
}

template class MESH_CUTTING<double>;
template class MESH_CUTTING<float>;
