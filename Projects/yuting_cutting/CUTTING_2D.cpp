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
CUTTING(TRIANGULATED_AREA<T>* sim_ta_,SEGMENTED_CURVE<TV>* sc_)
:sim_ta(sim_ta_),sc(sc_)
{
    //ta
    ta=TRIANGULATED_AREA<T>::Create();
    int p=sim_ta->particles.X.m;
    ta->particles.Resize(p);
    ta->Update_Number_Nodes();
    for(int i=0;i<p;++i){
        ta->particles.X(i)=sim_ta->particles.X(i);
    }
    ta->mesh.elements=sim_ta->mesh.elements;
    
    //tri_in_mesh
    tri_in_sim.Resize(ta->mesh.elements.m);
    for(int i=0;i<ta->mesh.elements.m;++i)
        tri_in_sim(i)=i;
    
    //particle_in_sim
    particle_in_sim.Resize(p);
    for(int i=0;i<ta->mesh.elements.m;++i)
        for(int j=0;j<3;++j){
            T3 w;
            w(j)=1;
            particle_in_sim(ta->mesh.elements(i)(j))=PS(i,w);
        }
    
    //tri_cuttings
    tri_cuttings.Resize(ta->mesh.elements.m);
}

//#####################################################################
// Function Run
//#####################################################################
template<class T> void CUTTING<VECTOR<T,2> >::
Run(T tol)
{
    cout<<"*********cutting**************"<<endl;
    //preprocessing
    cout<<"preprocessing"<<endl;
    CONSISTENT_INTERSECTIONS<TV> intersections(*ta,*sc);
    intersections.Set_Tol();
    intersections.Compute();
    
    HASHTABLE<I3,int> tri_from_face;
    HASHTABLE<I2,ARRAY<int> > tri_from_edge;
    ARRAY<ARRAY<int> > tri_from_vertex(ta->particles.number);

    for(int i=0;i<ta->mesh.elements.m;i++){
        I3 e=ta->mesh.elements(i).Sorted();
        tri_from_face.Set(e,i);
        for(int j=0;j<3;j++) tri_from_edge.Get_Or_Insert(e.Remove_Index(j)).Append(i);
        for(int j=0;j<3;j++) tri_from_vertex(e(j)).Append(i);}

    HASHTABLE<I2,int> seg_from_edge;
    ARRAY<ARRAY<int> > seg_from_vertex(sc->particles.number);

    for(int i=0;i<sc->mesh.elements.m;i++){
        I2 e=sc->mesh.elements(i).Sorted();
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
    ARRAY<int> parent_particles(IDENTITY_ARRAY<>(ta->particles.number));
    ARRAY<int> sim_parent_particles(IDENTITY_ARRAY<>(sim_ta->particles.number));
    HASHTABLE<int> split_tris, duplicated_sim_tris;
    HASHTABLE<int> dup_nodes;
    ARRAY<I3> original_sim_elements=sim_ta->mesh.elements;
    ARRAY<int> original_tri_in_sim=tri_in_sim;
    for(typename HASHTABLE<int,HASHTABLE<int,ARRAY<P> > >::ITERATOR it(components);it.Valid();it.Next()){
        int i=it.Key();
        TRI_CUTTING tc=tri_cuttings(i);
        I3 tri=ta->mesh.elements(i);
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
        cout << tc.turned_on << endl;
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
            I3 parent_tri=original_sim_elements(tri_in_sim(i));
            for(int j=0;j<a.m;++j){
                //duplicate sim tri
                for(int k=0;k<3;++k)
                    sim_parent_particles.Append(parent_tri(k));
                int parent_tri_id=original_tri_in_sim(i);
                I3 new_sim_tri(sim_parent_particles.m-3,sim_parent_particles.m-2,sim_parent_particles.m-1);
                if(duplicated_sim_tris.Contains(parent_tri_id)){
                    parent_tri_id=sim_ta->mesh.elements.m;
                    sim_ta->mesh.elements.Append(new_sim_tri);
                }
                else{
                    duplicated_sim_tris.Set(parent_tri_id);
                    sim_ta->mesh.elements(parent_tri_id)=new_sim_tri;
                }

                //split cutting tri, take care of parent-child relationship
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
                    parent_particles.Append(pid);
                    dup_nodes.Set(pid);
                    int s=particle_in_sim(pid).x;
                    I3 tri1=original_sim_elements(s);
                    T3 w1=particle_in_sim(pid).y;
                    I3 tri=original_sim_elements(parent_tri_id);
                    particle_in_sim.Append(PS(parent_tri_id,Weight_In_Tri(tri,tri1,w1)));
                }
                I3 new_tri(parent_particles.m-3,parent_particles.m-2,parent_particles.m-1);
                if(j==0){
                    tri_cuttings(i)=tc_new;
                    ta->mesh.elements(i)=new_tri;
                    split_tris.Set(i);
                    tri_in_sim(i)=parent_tri_id;
                }
                else{
                    tri_cuttings.Append(tc_new);
                    ta->mesh.elements.Append(new_tri);
                    split_tris.Set(tri_cuttings.m-1);
                    tri_in_sim.Append(parent_tri_id);
                }
            }
        }
    }
    
    
    //NOT DONE YET!!!
    //split triangles' neighbors' sharing node also need to duplicate, and parent needs to be cuplicated too
    for(int i=0;i<ta->mesh.elements.m;++i){
        I3& tri=ta->mesh.elements(i);
        if(!split_tris.Contains(i)){
            for(int j=0;j<3;++j){
                if(dup_nodes.Contains(tri(j))){
//                    parent_particles.Append(tri(j));
//                    tri(j)=parent_particles.m-1;
                    split_tris.Set(i);
                }
            }
        }
    }
    
    //union
    cout<<"merging"<<endl;
    HASHTABLE<I3,I3> ht;
    UNION_FIND<int> uf(parent_particles.m);
    UNION_FIND<int> sim_uf(sim_parent_particles.m);
    for(HASHTABLE_ITERATOR<int> it(split_tris);it.Valid();it.Next()){
        int tri_id=it.Key();
        I3 tri=ta->mesh.elements(tri_id);
        int parent_tri_id=tri_in_sim(tri_id);
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
                I3 saved(j1,j2,parent_tri_id);
                //cout << I2(i1, i2) << saved << endl;
                if(ht.Get(I3(i1,i2,i3),saved)){
                    //union material nodes
                    uf.Union(j1,saved(0));
                    uf.Union(j2,saved(1));
                    //union sim nodes
                    for(int j=0;j<3;++j){
                        int p1=sim_ta->mesh.elements(parent_tri_id)(j);
                        for(int k=0;k<3;++k){
                            int p2=sim_ta->mesh.elements(saved(2))(k);
                            if(sim_parent_particles(p1)==sim_parent_particles(p2))
                                sim_uf.Union(p1,p2);
                        }
                    }
                }
                else
                    ht.Set(I3(i1,i2,i3),saved);}}}

    //merge
    //sim_ta
    ARRAY<int> sim_tri_becomes(sim_ta->mesh.elements.m);
    {
        //delete unused nodes
        HASHTABLE<int,int> new_pids;
        int new_pid=0;
        ARRAY<TV> new_par;
        for(int i=0;i<sim_ta->mesh.elements.m;++i){
            for(int j=0;j<3;++j){
                int a=sim_uf.Find(sim_ta->mesh.elements(i)(j));
                if(!new_pids.Get(a,sim_ta->mesh.elements(i)(j))){
                    new_par.Append(sim_ta->particles.X(sim_parent_particles(a)));
                    sim_ta->mesh.elements(i)(j)=new_pid;
                    new_pids.Set(a,new_pid);
                    ++new_pid;}}}
        sim_ta->particles.Resize(new_par.m);
        sim_ta->particles.X=new_par;
        
        //delete duplicated sim triangles
        HASHTABLE<I3,int> tri_becomes;
        int new_tid=0;
        for(int i=0;i<sim_ta->mesh.elements.m;++i){
            I3 tri=sim_ta->mesh.elements(i);
            if(!tri_becomes.Get(tri,sim_tri_becomes(i))){
                tri_becomes.Set(tri,new_tid);
                sim_tri_becomes(i)=new_tid;
                sim_ta->mesh.elements(new_tid)=tri;
                ++new_tid;
            }
        }
        sim_ta->mesh.elements.Resize(new_tid);
    }
    
    //ta
    {
        HASHTABLE<int,int> new_pids;
        int new_pid=0;
        ARRAY<PS> new_particle_in_sim;
        for(int i=0;i<tri_cuttings.m;++i){
            tri_in_sim(i)=sim_tri_becomes(tri_in_sim(i));
            for(int j=0;j<3;++j){
                int& n=ta->mesh.elements(i)(j);
                int a=uf.Find(n);
                if(!new_pids.Get(a,n)){
                    PS ps=particle_in_sim(n);
                    new_particle_in_sim.Append(PS(sim_tri_becomes(ps.x),ps.y));
                    n=new_pid;
                    new_pids.Set(a,new_pid);
                    ++new_pid;
                }}}
        particle_in_sim=new_particle_in_sim;
    }
    
    //subdivide and delete unused nodes
    {
        cout<<"subdividing"<<endl;
        HASHTABLE<I2,int> new_edge_particles;
        ARRAY<I3> new_elements;
        ARRAY<int> new_tri_in_sim;
        //subdivide split tris
        for(int i=0;i<tri_cuttings.m;++i){
            TRI_CUTTING tc=tri_cuttings(i);
            if(tc.materials.Find(0)!=-1){
                I3 tri=ta->mesh.elements(i);
                int p=particle_in_sim.m;
                particle_in_sim.Append(PS(tri_in_sim(i),Weight_In_Sim(i,tc.face_center.Value())));
                for(int j=0;j<3;++j){
                    if(tc.materials(2*j)||tc.materials(2*j+1)){
                        int q;
                        I2 e(tri(j),tri((j+1)%3));
                        if(!new_edge_particles.Get(e.Sorted(),q)){
                            new_edge_particles.Set(e.Sorted(),particle_in_sim.m);
                            T3 w;
                            for(int k=0;k<2;++k)
                                w[(j+k)%3]=tc.edge_centers(j).Value()(k);
                            particle_in_sim.Append(PS(tri_in_sim(i),Weight_In_Sim(i,w)));
                        }
                        if(tc.materials(2*j)){
                            new_elements.Append(I3(e(0),q,p));
                            new_tri_in_sim.Append(tri_in_sim(i));
                        }
                        if(tc.materials(2*j+1)){
                            new_elements.Append(I3(q,e(1),p));
                            new_tri_in_sim.Append(tri_in_sim(i));
                        }}}}}
        //subdivide neighbors
        for(int i=0;i<tri_cuttings.m;++i){
            I3 tri=ta->mesh.elements(i);
            if(tri_cuttings(i).materials.Find(0)==-1){
                int p=-1;
                int n=new_elements.m;
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
                        break;
                    }
                }
                if(p==-1)
                    new_elements.Append(tri);
                for(int j=n;j<new_elements.m;++j)
                    new_tri_in_sim.Append(tri_in_sim(i));
            }
        }
        //reset per-triangle data
        ta->mesh.elements=new_elements;
        tri_in_sim=new_tri_in_sim;
        tri_cuttings.Remove_All();
        tri_cuttings.Resize(tri_in_sim.m);
        //get rid of unused particles produced by subdivision
        HASHTABLE<int,int> new_pids;
        int new_pid=0;
        ARRAY<PS> new_particle_in_sim;
        for(int i=0;i<ta->mesh.elements.m;++i)
            for(int j=0;j<3;++j){
                int& id=ta->mesh.elements(i)(j);
                if(!new_pids.Get(id,id)){
                    new_particle_in_sim.Append(particle_in_sim(id));
                    new_pids.Set(id,new_pid); 
                    id=new_pid;
                    ++new_pid;}}
        particle_in_sim=new_particle_in_sim;
        Update_Material_Particles();
    }
    
    cout<<"*********cutting done**************"<<endl;
}

//#####################################################################
// Function Update_Material_Particles
//#####################################################################
template<class T> void CUTTING<VECTOR<T,2> >::
Update_Material_Particles()
{
    int p=particle_in_sim.m;
    ta->particles.Resize(p);
    ta->Update_Number_Nodes();
    for(int i=0;i<p;++i){
        T3 w=particle_in_sim(i).y;
        I3 tri=sim_ta->mesh.elements(particle_in_sim(i).x);
        TV pa;
        for(int j=0;j<3;++j){
            pa+=(w(j)*sim_ta->particles.X(tri(j)));
        }
        ta->particles.X(i)=pa;
    }
}
//#####################################################################
template class CUTTING<VECTOR<double,2> >;
