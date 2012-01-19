//#####################################################################
// Copyright 2002, Robert Bridson.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
//
//#####################################################################
// Bridson - November 20, 2002
//#####################################################################
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <fstream>
#include <iostream>
#include "float.h"

typedef float GEOMETRY_TYPE;

using namespace std;
using namespace PhysBAM;

class SURFACE_CHECKER
{
public:
    const char *filename;
    TRIANGLE_MESH mesh;
    PARTICLE_3D<GEOMETRY_TYPE> particles;
    TRIANGULATED_SURFACE<GEOMETRY_TYPE> surface;

    explicit SURFACE_CHECKER(const char* filename_input)
        :filename(filename_input),surface(mesh,particles)
    {}

    ~SURFACE_CHECKER()
    {}

//#######################################################################
// Function Run
//#######################################################################
void Run()
{
    cout<<"Input filename: "<<filename<<endl;
    ifstream input(filename,ios::binary);
    if(!input) {cout<<"#### Could not open file"<<endl;return;}
    surface.Read(input);
    bool surface_problems=false;

    // check particles and mesh match up
    if(particles.number != mesh.number_nodes) {cout<<"#### "<<endl;surface_problems=true;}
    cout<<"particles.number="<<particles.number<<" and mesh.number_nodes="<<mesh.number_nodes<<endl;

    cout<<"number of triangles is "<<mesh.triangles.m<<endl;

    int i,j,k;
    // check particle positions for NaN
    for(i=0;i<particles.number;i++) if(particles.active(i)){
        if(!Number_Is_Finite(particles.X(i).x) || !Number_Is_Finite(particles.X(i).y) || !Number_Is_Finite(particles.X(i).z)){
            cout<<"#### particle "<<i<<" has suspicious coordinates "<<particles.X(i)<<endl;surface_problems=true;}}

    // check for boundary
    mesh.Initialize_Boundary_Mesh();
    if(mesh.boundary_mesh->segments.m) cout<<"???? mesh has "<<mesh.boundary_mesh->segments.m<<" boundary segments."<<endl;

    // check for non-manifold edges (incident to 3 or more triangles) and orientation problems
    mesh.Initialize_Segment_Mesh();
    mesh.Initialize_Edge_Triangles();
    for(i=0;i<mesh.segment_mesh->segments.m;i++){
        int p=mesh.segment_mesh->segments(1,i),q=mesh.segment_mesh->segments(2,i);
        if((*mesh.edge_triangles)(i).m>=3){
            cout<<"#### edge "<<i<<" ("<<p<<","<<q<<") has too many incident triangles:"<<endl;
            for(j=0;j<(*mesh.edge_triangles)(i).m;j++){
                cout<<"  triangle "<<(*mesh.edge_triangles)(i)(j)<<" (";
                for(k=0;k<3;k++) cout<<" "<<mesh.triangles(k,(*mesh.edge_triangles)(i)(j));
                cout<<")"<<endl;}
            surface_problems=true;}
        else if((*mesh.edge_triangles)(i).m==2){
            int t1=(*mesh.edge_triangles)(i)(1),t2=(*mesh.edge_triangles)(i)(2);
            bool ascending1,ascending2;
            if(mesh.triangles(1,t1)==p){
                if(mesh.triangles(2,t1)==q) ascending1=true;
                else if(mesh.triangles(3,t1)==q) ascending1=false;
                else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}}
            else if(mesh.triangles(2,t1)==p){
                if(mesh.triangles(3,t1)==q) ascending1=true;
                else if(mesh.triangles(1,t1)==q) ascending1=false;
                else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}}
            else if(mesh.triangles(3,t1)==p){
                if(mesh.triangles(1,t1)==q) ascending1=true;
                else if(mesh.triangles(2,t1)==q) ascending1=false;
                else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}}
            else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}
            if(mesh.triangles(1,t2)==p){
                if(mesh.triangles(2,t2)==q) ascending2=true;
                else if(mesh.triangles(3,t2)==q) ascending2=false;
                else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}}
            else if(mesh.triangles(2,t2)==p){
                if(mesh.triangles(3,t2)==q) ascending2=true;
                else if(mesh.triangles(1,t2)==q) ascending2=false;
                else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}}
            else if(mesh.triangles(3,t2)==p){
                if(mesh.triangles(1,t2)==q) ascending2=true;
                else if(mesh.triangles(2,t2)==q) ascending2=false;
                else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}}
            else {cout<<"BUG IN EDGE_TRIANGLES CODE!"<<endl;return;}
            if(ascending1 == ascending2) {
                cout<<"#### triangles "<<t1<<" and "<<t2<<" disagree on orientation for common edge ("<<p<<","<<q<<")"<<endl;
                surface_problems=true;}}}

    // double check neighbor_nodes vs. toplogically sorted neighbor_nodes
    mesh.Initialize_Neighbor_Nodes();
    assert(mesh.neighbor_nodes->m==mesh.number_nodes);
    ARRAY<ARRAY<int> > &nbrs=*mesh.neighbor_nodes;
    mesh.neighbor_nodes=0;
    mesh.Initialize_Topologically_Sorted_Neighbor_Nodes();
    assert(mesh.neighbor_nodes->m==mesh.number_nodes);
    ARRAY<ARRAY<int> > &top_nbrs=*mesh.neighbor_nodes;
    for(i=0;i<mesh.number_nodes;i++){
        if(nbrs(i).m!=top_nbrs(i).m){
            cout<<"Topologically sorted neighbors has "<<top_nbrs(i).m<<" entries for node "<<i<<" versus ";
            cout<<nbrs(i).m<<" for regular neighbors"<<endl;
            surface_problems=true;}
        else{
            for(j=1; j<=nbrs(i).m; ++j) if(!top_nbrs(i).Find(nbrs(i)(j),k)) {
                cout<<"Topologically sorted neighbors doesn't have neighbor "<<nbrs(i)(j)<<" for node "<<i<<endl;
                surface_problems=true;}
            for(j=1; j<=top_nbrs(i).m; ++j) if(!nbrs(i).Find(top_nbrs(i)(j),k)) {
                cout<<"Topologically sorted neighbors has neighbor "<<nbrs(i)(j)<<" for node "<<i<<" which is false"<<endl;
                surface_problems=true;}}}

    // check for self-intersections
    surface.Initialize_Triangle_Hierarchy();
    ARRAYS<VECTOR<int,1> > triangle_pairs(2,1,0);
    surface.Update_Bounding_Box();
    double thickness=1e-8*surface.bounding_box->Size().Magnitude();
    if(surface.Check_For_Self_Intersection(thickness,&triangle_pairs)){
        surface_problems=true;
        cout<<"Surface self-intersects for the following pairs of triangles"<<endl;
        for(j=0;j<triangle_pairs.m;j++) cout<<"  "<<triangle_pairs(1,j)<<" "<<triangle_pairs(2,j)<<endl;}

    if(!surface_problems) cout<<filename<<" appears to be ok"<<endl;
    else cout<<"#### "<<filename<<" has problems!"<<endl;
}
//#######################################################################
// Function Number_Is_Finite
//#######################################################################
bool Number_Is_Finite(double number)
{
    return number<=FLT_MAX && number>=-DBL_MAX && number==number; // check for inf,-inf, and NaN
}
//#######################################################################
};
//#######################################################################
// Function main
//#######################################################################
int main(int argc,char *argv[])
{
    if(argc<=1){
        SURFACE_CHECKER checker("../Public_Library/Data/Rigid_Bodies/Bones/Cranium.tri");
        checker.Run();}
    else for(int i=1; i<argc; i++){
        if(i>1) cout<<"\n----------------------------------------------------------------------------"<<endl;
        SURFACE_CHECKER checker(argv[i]);
        checker.Run();}
    return 0;
}
