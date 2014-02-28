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
    cout<<"preprocessing"<<endl;
    CONSISTENT_INTERSECTIONS<TV> intersections(ta,sc);
    intersections.Set_Tol();
    intersections.Compute();
    
    HASHTABLE<I3,int> tri_from_face;
    HASHTABLE<I2,ARRAY<int> > tri_from_edge;
    ARRAY<ARRAY<int> > tri_from_vertex(ta.particles.number);

    for(int i=0;i<ta.mesh.elements.m;i++){
        I3 e=ta.mesh.elements(i).Sorted();
        tri_from_face.Set(e,i);
        for(int j=0;j<3;j++) tri_from_edge.Get_Or_Insert(e.Remove_Index(j)).Append(i);
        for(int j=0;j<3;j++) tri_from_vertex(e(j)).Append(i);}

    HASHTABLE<I2,int> seg_from_edge;
    ARRAY<ARRAY<int> > seg_from_vertex(sc.particles.number);

    for(int i=0;i<sc.mesh.elements.m;i++){
        I2 e=sc.mesh.elements(i).Sorted();
        seg_from_edge.Set(e,i);
        for(int j=0;j<2;j++) seg_from_vertex(e(j)).Append(i);}

    typedef PAIR<I5,T3> P;
    HASHTABLE<int,HASHTABLE<int,ARRAY<P> > > components;

    for(typename HASHTABLE<I2>::ITERATOR it(intersections.hash_vv);it.Valid();it.Next()){
        I2 key=it.Key();
        const ARRAY<int>& a=tri_from_vertex(key.x);
        const ARRAY<int>& b=seg_from_vertex(key.y);
        for(int i=0;i<a.m;i++)
            for(int j=0;j<b.m;j++)
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Append(P(I5(key.x,-1,-1,key.y,-1),T3()));}

    for(typename HASHTABLE<I3,T>::ITERATOR it(intersections.hash_ve);it.Valid();it.Next()){
        I3 key=it.Key();
        const ARRAY<int>& a=tri_from_vertex(key.x);
        int b=seg_from_edge.Get(key.Remove_Index(0));
        for(int i=0;i<a.m;i++)
            components.Get_Or_Insert(a(i)).Get_Or_Insert(b).Append(P(I5(key.x,-1,-1,key.y,key.z),T3(it.Data(),0,0)));}

    for(typename HASHTABLE<I3,T>::ITERATOR it(intersections.hash_ev);it.Valid();it.Next()){
        I3 key=it.Key();
        const ARRAY<int>& a=tri_from_edge.Get(key.Remove_Index(2));
        const ARRAY<int>& b=seg_from_vertex(key.z);
        for(int i=0;i<a.m;i++)
            for(int j=0;j<b.m;j++)
                components.Get_Or_Insert(a(i)).Get_Or_Insert(b(j)).Append(P(I5(key.x,key.y,-1,key.z,-1),T3(it.Data(),0,0)));}

    for(typename HASHTABLE<I4,TV>::ITERATOR it(intersections.hash_ee);it.Valid();it.Next()){
        I4 key=it.Key();
        const ARRAY<int>& a=tri_from_edge.Get(I2(key(0),key(1)));
        int b=seg_from_edge.Get(I2(key(2),key(3)));
        for(int i=0;i<a.m;i++)
            components.Get_Or_Insert(a(i)).Get_Or_Insert(b).Append(P(key.Insert(-1,2),it.Data().Append(0)));}

    for(typename HASHTABLE<I4,T3>::ITERATOR it(intersections.hash_fv);it.Valid();it.Next()){
        I4 key=it.Key();
        int a=tri_from_face.Get(key.Remove_Index(3));
        const ARRAY<int>& b=seg_from_vertex(key.Last());
        for(int j=0;j<b.m;j++)
            components.Get_Or_Insert(a).Get_Or_Insert(b(j)).Append(P(key.Append(-1),it.Data()));}

    //split
    cout<<"splitting"<<endl;
    int num_old_tris=ta.mesh.elements.m;
    ARRAY<int> parent_particles(IDENTITY_ARRAY<>(ta.particles.number));
    HASHTABLE<int> split_tris;
    tri_cuttings.Resize(num_old_tris);
    for(typename HASHTABLE<int,HASHTABLE<int,ARRAY<P> > >::ITERATOR it(components);it.Valid();it.Next()){
        int i=it.Key();
        TRI_CUTTING tc=tri_cuttings(i);
        I3 tri=ta.mesh.elements(i);
        //turn on by each cutting segment based on itersections
        for(typename HASHTABLE<int,ARRAY<P> >::ITERATOR seg_it(it.Data());seg_it.Valid();seg_it.Next()){
            const ARRAY<P>& intersects=seg_it.Data();
            VECTOR<bool,7> hit;
            for(int k=0;k<intersects.m;++k){
                const I5& p=intersects(k).x;
                const T3& weight=intersects(k).y;
                if(p(1)==-1){
                    hit(tri.Find(p(0))*2)=1;
                    T3 c;
                    c(tri.Find(p(0)))=1;
                    tc.face_center.Add(c);
                }
                else if(p(2)==-1){
                    I2 e(tri.Find(p(0)),tri.Find(p(1)));
                    e.Sort();
                    int eid=-1;
                    if(e(0)==0){
                        if(e(1)==1){
                            hit(1)=1;
                            eid=0;}
                        else{
                            hit(5)=1;
                            eid=2;}}
                    else{
                        hit(3)=1;
                        eid=1;}
                    if(tri(eid)==p(0)){
                        tc.edge_centers(eid).Add(TV(1-weight(0),weight(0)));
                        T3 c;
                        c(eid)=1-weight(0);
                        c((eid+1)%3)=weight(0);
                        tc.face_center.Add(c);
                    }
                    else{
                        tc.edge_centers(eid).Add(TV(weight(0),1-weight(0)));
                        T3 c;
                        c(eid)=weight(0);
                        c((eid+1)%3)=1-weight(0);
                        tc.face_center.Add(c);
                    }
                }
                else{
                    hit(6)=1;
                    T3 c;
                    for(int j=0;j<3;++j)
                        c(tri.Find(p(j)))=weight(j);
                    tc.face_center.Add(c);}}
            //turn on
            for(int k=0;k<6;++k)
                if(hit(k)&&(hit((k+3)%6)||hit(6)||(k+4)%6))
                    tc.turned_on(k)=1;
            for(int k=0;k<3;++k){
                int k1=2*k;
                int k2=2*k+1;
                int k3=(2*k+2)%6;
                if(hit(k1)&&(hit(k2)||hit(k3)))
                    tc.turned_on(k1+6)=1;
                if(hit(k3)&&(hit(k1)||hit(k2)))
                    tc.turned_on(k2+6)=1;
            }
        }
        
        //split based on turn-on and intersections
        ARRAY<int> a;
        for(int j=0;j<6;++j)
            if (tc.turned_on(j))
                a.Append(j);
        if(a.m>1){
            //set centers
            tc.face_center.set=true;
            for(int j=0;j<a.m;++j){
                if(a(j)%2)
                    tc.edge_centers(a(j)/2).set=true;}
            //split
            for(int j=0;j<a.m;++j){
                TRI_CUTTING tc_new=tc;
                //material of new element
                tc_new.materials.Fill(0);
                if (j==a.m-1){
                    for(int k=a(j);k<6;++k)
                        tc_new.materials(k)=1;
                    for(int k=0;k<a(0);++k)
                        tc_new.materials(k)=1;}
                else
                    for(int k=a(j);k<a(j+1);++k)
                        tc_new.materials(k)=1;
                //new particles and element
                for(int k=0;k<3;++k){
                    int pid=tri(k);
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
            
        }
    }
    //neighbors' node also need to duplicate
    HASHTABLE<int> dup_nodes;
    for(HASHTABLE_ITERATOR<int> it(split_tris);it.Valid();it.Next()){
        int i=it.Key();
        I3 tri=ta.mesh.elements(i);
        for(int j=0;j<3;++j){
            dup_nodes.Set(parent_particles(tri(j)));
        }    
    }
    for(int i=0;i<ta.mesh.elements.m;++i){
        I3& tri=ta.mesh.elements(i);
        if(!split_tris.Contains(i)){
            for(int j=0;j<3;++j){
                if(dup_nodes.Contains(tri(j))){
                    int p=ta.particles.Add_Element();
                    ta.particles.X(p)=ta.particles.X(tri(j));
                    parent_particles.Append(tri(j));
                    tri(j)=p;
                }
            }
            split_tris.Set(i);
        }
    }
    
    //union
    cout<<"merging"<<endl;
    HASHTABLE<I3,I2> ht;
    UNION_FIND<int> uf(parent_particles.m);
    for(HASHTABLE_ITERATOR<int> it(split_tris);it.Valid();it.Next()){
        int tri_id=it.Key();
        I3 tri=ta.mesh.elements(tri_id);
        for(int i=0;i<6;++i){
            if(tri_cuttings(tri_id).materials(i) && !tri_cuttings(tri_id).turned_on(i+6)){
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
                    j2=temp;
                    i3=1-i3;}
                I2 saved(j1,j2);
                //cout << I2(i1, i2) << saved << endl;
                if(ht.Get(I3(i1,i2,i3),saved)){
                    uf.Union(j1,saved(0));
                    uf.Union(j2,saved(1));}
                else
                    ht.Set(I3(i1,i2,i3),saved);}}}

    //merge
    HASHTABLE<int,int> new_pids;
    int new_pid=0;
    ARRAY<TV> new_par;
    for(int i=0;i<tri_cuttings.m;++i){
        for(int j=0;j<3;++j){
            int a=uf.Find(ta.mesh.elements(i)(j));
            if(!new_pids.Get(a,ta.mesh.elements(i)(j))){
                new_par.Append(ta.particles.X(ta.mesh.elements(i)(j)));
                ta.mesh.elements(i)(j)=new_pid;
                new_pids.Set(a,new_pid);
                ++new_pid;}}}
    ta.particles.Resize(new_par.m);
    ta.particles.X=new_par;
    
    //subdivide
    cout<<"subdividing"<<endl;
    HASHTABLE<I2,int> new_edge_particles;
    HASHTABLE<I3,int> new_tri_particles;
    ARRAY<I3> new_elements;
    //subdivide split tris
    for(int i=0;i<tri_cuttings.m;++i){
        TRI_CUTTING tc=tri_cuttings(i);
        if(tc.materials.Find(0)!=-1){
            I3 tri=ta.mesh.elements(i);
            int p;
            if(!new_tri_particles.Get(tri.Sorted(),p)){
                new_tri_particles.Set(tri.Sorted(),p);
                TV par;
                for(int j=0;j<3;++j)
                    par+=ta.particles.X(tri(j))*tc.face_center.Value()(j);
                p=ta.particles.Add_Element();
                ta.particles.X(p)=par;}
            for(int j=0;j<3;++j){
                if(tc.materials(2*j)||tc.materials(2*j+1)){
                    int q;
                    I2 e(tri(j),tri((j+1)%3));
                    if(!new_edge_particles.Get(e.Sorted(),q)){
                        TV par;
                        for(int k=0;k<2;++k){
                            par+=ta.particles.X(e(k))*tc.edge_centers(j).Value()(k);
                        }
                        q=ta.particles.Add_Element();
                        new_edge_particles.Set(e.Sorted(),q);
                        ta.particles.X(q)=par;}
                    if(tc.materials(2*j))
                        new_elements.Append(I3(e(0),q,p));
                    if(tc.materials(2*j+1))
                        new_elements.Append(I3(q,e(1),p));}}}}
    //subdivide neighbors
    for(int i=0;i<tri_cuttings.m;++i){
        I3 tri=ta.mesh.elements(i);
        if(tri_cuttings(i).materials.Find(0)==-1){
            int p=-1;
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
                    break;}}
            if(p==-1)
                new_elements.Append(tri);}}
    //get rid of unused particles
    new_pids.Clean_Memory();
    new_pid=0;
    new_par.Remove_All();
    for(int i=0;i<new_elements.m;++i)
        for(int j=0;j<3;++j){
            int& id=new_elements(i)(j);
            if(!new_pids.Get(id,id)){
                new_par.Append(ta.particles.X(id));
                new_pids.Set(id,new_pid); 
                id=new_pid;
                ++new_pid;}}
    ta.particles.Resize(new_par.m);
    ta.particles.X=new_par;
    ta.mesh.elements=new_elements;
    ta.Update_Number_Nodes();
}
//#####################################################################
template class CUTTING<VECTOR<double,2> >;
