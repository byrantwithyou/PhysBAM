//#####################################################################
// Copyright 2014, Yuting Wang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_2D
//#####################################################################
#include <Tools/Math_Tools/cube.h>
#include <Tools/Matrices/MATRIX.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Tools/Data_Structures/UNION_FIND.h>
#include "CUTTING_2D.h"
#include "CONSISTENT_INTERSECTIONS.h"
using namespace PhysBAM;
using namespace std;
//#####################################################################
// Constructor
//#####################################################################
template<class T> CUTTING<VECTOR<T,2> >::
CUTTING(TRIANGULATED_AREA<T>& ta,SEGMENTED_CURVE<TV>& sc)
    :ta(ta),sc(sc)
{
}
//#####################################################################
// Function Run
//#####################################################################
template<class T> void CUTTING<VECTOR<T,2> >::
Run(T tol)
{
    //preprocessing
    CONSISTENT_INTERSECTIONS<TV> intersections(ta,sc);
    intersections.Set_Tol();
    intersections.Compute();

    cout << "here" << endl;
    ARRAY<ARRAY<int> > inter_list;
    BOX_VISITOR_TRIVIAL v(inter_list);
    int num_old_tris=ta.mesh.elements.m;
    inter_list.Resize(num_old_tris);
    ARRAY<int> parent_particles(num_old_tris);
    for(int i=0;i<num_old_tris;++i)
        parent_particles(i)=i;
    ta.hierarchy->Intersection_List(*sc.hierarchy,v,tol);
    
    //turn on
    HASHTABLE<int> split_tris;
    for(int i=0;i<num_old_tris;i++){
        TRI_CUTTING tc=tri_cuttings(i);
        I3 tri=ta.mesh.elements(i);
        //turn on by each cutting segment based on itersections
        for(int j=0;j<inter_list(i).m;j++){
            I2 e=sc.mesh.elements(inter_list(i)(j)).Sorted();
            //hit nodes
            VECTOR<bool,7> hit;
            for(int k=0;k<3;++k){
                int hit_id=2*k;
                for(int l=0;l<2;++l)
                    if (intersections.hash_vv.Contains(I2(tri(k),e(l))))
                        hit(hit_id)=1;
                if(intersections.hash_ve.Contains(I3(tri(k),e(0),e(1))))
                    hit(hit_id)=1;}
            //hit edges
            for(int k=0;k<3;++k){
                int hit_id=2*k+1;
                I2 te(tri(k),tri((k+1)%3));
                te.Sort();
                for(int l=0;l<2;++l)
                    if(intersections.hash_ev.Contains(te.Append(e(l))))
                        hit(hit_id)=1;
                if(intersections.hash_ee.Contains(te.Append_Elements(e)))
                    hit(hit_id)=1;}
            //hit interior
            if (intersections.hash_fv.Contains(tri.Sorted().Append(e(0)))||intersections.hash_fv.Contains(tri.Sorted().Append(e(1))))
                hit(6)=1;
            //turn on
            for(int k=0;k<6;++k)
                if(hit(k)&&(hit((k+3)%6)||hit(6)))
                    tc.turned_on(k)=1;
            for(int k=0;k<3;++k){
                int k1=2*k;
                int k2=2*k+1;
                int k3=(2*k+2)%6;
                if(hit(k1)&&(hit(k2)||hit(k3)))
                    tc.turned_on(k1+6)=1;
                if(hit(k3)&&(hit(k1)||hit(k2)))
                    tc.turned_on(k2+6)=1;}}
        
        //split based on turn-on and intersections
        ARRAY<int> a;
        for(int j=0;j<6;++j)
            if (tc.turned_on(j))
                a.Append(j);
        if(a.m>1){
            for(int j=0;j<a.m;++j){
                TRI_CUTTING tc_new=tc;
                //material of new element
                tc_new.material.Remove_All();
                if (j==a.m-1){
                    for(int k=a(j);k<6;++k)
                        tc_new.material.Append(k);
                    for(int k=0;k<a(0);++k)
                        tc_new.material.Append(k);}
                else
                    for(int k=a(j);k<a(j+1);++k)
                        tc_new.material.Append(k);
                //new particles and element
                for(int k=0;k<3;++k){
                    int pid=ta.mesh.elements(i)(k);
                    ta.particles.Append(ta.particles,pid);
                    parent_particles.Append(pid);}
                I3 new_tri(ta.particles.X.m-3,ta.particles.X.m-2,ta.particles.X.m-1);
                if(j==0){
                    tri_cuttings(i)=tc_new;
                    ta.mesh.elements(i)=new_tri;
                    split_tris.Set(i);}
                else{
                    tri_cuttings.Append(tc_new);
                    ta.mesh.elements.Append(new_tri);
                    split_tris.Set(tri_cuttings.m-1);}}
            //neighbors also need to duplicate
            
        }
    }
    
    //union
    HASHTABLE<I3,I2> ht;
    UNION_FIND<int> uf(parent_particles.m);
    for(HASHTABLE_ITERATOR<int> it(split_tris);it.Valid();it.Next()){
        int tri_id=it.Key();
        I3 tri=ta.mesh.elements(tri_id);
        for(int i=0;i<6;++i){
            if(tri_cuttings(tri_id).material.Contains(i) && !tri_cuttings(tri_id).turned_on(i+6)){
                int j1=tri(i/2);
                int j2=tri((i/2+1)%3);
                int i1=parent_particles(j1);
                int i2=parent_particles(j2);
                int i3=i%2;
                if(i1>i2){
                    int temp=i1;
                    i1=i2;
                    i2=temp;
                    temp=j1;
                    j1=j2;
                    j2=j1;
                    i3=1-i3;}
                I2 saved;
                if(ht.Get(I3(i1,i2,i3),saved)){
                    uf.Union(i1,saved(0));
                    uf.Union(i2,saved(1));}
                else
                    ht.Set(I3(i1,i2,i3),I2(j1,j2));}}}
    //merge
    HASHTABLE<int,int> new_pids;
    int new_pid=0;
    ARRAY<TV> new_par;
    for(int i=0;i<tri_cuttings.m;++i){
        for(int j=0;j<3;++j){
            int a=uf.Find(ta.mesh.elements(i)(j));
            if(new_pids.Contains(a))
                ta.mesh.elements(i)(j)=new_pids.Get(a);
            else{
                new_par.Append(ta.particles.X(ta.mesh.elements(i)(j)));
                ta.mesh.elements(i)(j)=new_pid;
                new_pids.Set(a,new_pid);
                ++new_pid;}}}
    ta.particles.X=new_par;
    
    //subdivide
    HASHTABLE<I2,int> new_edge_particles;
    HASHTABLE<I3,int> new_tri_particles;
    ARRAY<I3> new_elements;
    for(int i=0;i<tri_cuttings.m;++i){
        TRI_CUTTING tc=tri_cuttings(i);
        if(tc.material.m<6){
            I3 tri=ta.mesh.elements(i);
            int p;
            if(!new_tri_particles.Get(tri.Sorted(),p)){
                p=ta.particles.X.m;
                new_tri_particles.Set(tri.Sorted(),p);
                TV par;
                for(int j=0;j<3;++j)
                    par+=ta.particles.X(tri(i))*tc.face_center.Value()(j);
                int pp=ta.particles.Add_Element();
                ta.particles.X(pp)=par;}
            for(int j=0;j<3;++j){
                int q;
                I2 e(tri(j),tri((j+1)%3));
                if(!new_edge_particles.Get(e.Sorted(),q)){
                    q=ta.particles.X.m;
                    new_edge_particles.Get(e.Sorted(),p);
                    TV par;
                    for(int k=0;k<2;++k)
                        par+=ta.particles.X(e(k))*tc.edge_centers(j).Value()(k);
                    int pp=ta.particles.Add_Element();
                    ta.particles.X(pp)=par;}
                if(tc.material.Contains(2*j))
                    new_elements.Append(I3(e(0),q,p));
                if(tc.material.Contains(2*j+1))
                    new_elements.Append(I3(p,e(1),q));}}}
    for(int i=0;i<tri_cuttings.m;++i){
        I3 tri=ta.mesh.elements(i);
        if(tri_cuttings(i).material.m==6){
            int p;
            for(int j=0;j<3;++j){
                I2 e(tri(j),tri((j+1)%3));
                if(new_edge_particles.Get(e.Sorted(),p)){
                    int n=tri((j+2)%3);
                    int q;
                    if(new_edge_particles.Get(I2(e(0),n).Sorted(),q)){
                        new_elements.Append(I3(n,q,p));
                        new_elements.Append(I3(q,e(0),p));}
                    else
                        new_elements.Append(I3(e(0),p,n));
                    if(new_edge_particles.Get(I2(e(1),n).Sorted(),q)){
                        new_elements.Append(I3(n,p,q));
                        new_elements.Append(I3(q,p,e(1)));}
                    else
                        new_elements.Append(I3(p,e(1),n));
                    break;}}}
        else
            new_elements.Append(tri);}
    ta.mesh.elements=new_elements;
    ta.Update_Number_Nodes();
}
//#####################################################################
template class CUTTING<VECTOR<double,2> >;
