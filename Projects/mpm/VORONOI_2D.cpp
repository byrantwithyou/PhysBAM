//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "TIMING.h"
#include "VORONOI_2D.h"
namespace PhysBAM{
//#####################################################################
// Function Initialize_With_A_Regular_Grid_Of_Particles
//#####################################################################
template<class T> void VORONOI_2D<T>::
Initialize_With_A_Regular_Grid_Of_Particles(const GRID<TV>& particle_grid)
{
    TIMING_START;
    GRID<TV> node_grid(TV_INT(particle_grid.counts.x*2+1,particle_grid.counts.y*2+1),RANGE<TV>(particle_grid.domain.min_corner-particle_grid.dX/2.0,TV(particle_grid.domain.max_corner+particle_grid.dX/2.0)));
    HASHTABLE<TV_INT,int> particle_index;
    int ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+particle_grid.counts));it.Valid();it.Next()) particle_index.Get_Or_Insert(it.index)=ID++;
    HASHTABLE<TV_INT,int> node_index;
    ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+node_grid.counts));it.Valid();it.Next()){
        int n_old_component=(int)(it.index.x%2==1)+(int)(it.index.y%2==1);
        if(n_old_component==0){
            node_index.Get_Or_Insert(it.index)=ID++;
            Xm.Append(node_grid.Node(it.index));
            if(it.index.x==0||it.index.x==node_grid.counts.x-1||it.index.y==0||it.index.y==node_grid.counts.y-1) type.Append(1);
            else type.Append(10);}
        else if(n_old_component==1){
            node_index.Get_Or_Insert(it.index)=ID++;
            Xm.Append(node_grid.Node(it.index));
            type.Append(100);}}
    ID=0;
    elements.Resize(particle_grid.counts.Product());
    for(int i=0;i<node_grid.counts.x-1;i+=2){
        for(int j=0;j<node_grid.counts.y-1;j+=2){
            TV_INT aa(i,j),bb(i+1,j),cc(i+2,j),dd(i+2,j+1),ee(i+2,j+2),ff(i+1,j+2),gg(i,j+2),hh(i,j+1);
            int a=node_index.Get_Or_Insert(aa),b=node_index.Get_Or_Insert(bb),c=node_index.Get_Or_Insert(cc),d=node_index.Get_Or_Insert(dd),e=node_index.Get_Or_Insert(ee),f=node_index.Get_Or_Insert(ff),g=node_index.Get_Or_Insert(gg),h=node_index.Get_Or_Insert(hh);
            elements(ID).Append(a);elements(ID).Append(b);elements(ID).Append(c);elements(ID).Append(d);elements(ID).Append(e);elements(ID).Append(f);elements(ID).Append(g);elements(ID).Append(h);
            ID++;}}
    X=Xm;
    TIMING_END("Initialize polygon voronoi mesh");
    Initialize_Neighbor_Cells();
    Build_Association();
}
//#####################################################################
// Function Initialize_With_A_Triangulated_Area
//#####################################################################
template<class T> void VORONOI_2D<T>::
Initialize_With_A_Triangulated_Area(const TRIANGULATED_AREA<T>& ta)
{
    //*** Preprocessing
    HASHTABLE<TV_INT,bool> triedge2type; // true-boundary,false-interior
    ARRAY<bool> particle2type(ta.particles.number); // true-boundary,false-interior
    ARRAY<int> boundary_particles;
    ARRAY<int> interior_particles;
    HASHTABLE<int,ARRAY<int> > particle2neighbortris;
    
    HASHTABLE<TV_INT,int> edge_count;
    for(int i=0;i<ta.mesh.elements.m;i++){
        for(int j=0;j<3;j++){
            int jp1=(j==2)?0:(j+1);
            int a=ta.mesh.elements(i)(j);
            particle2neighbortris.Get_Or_Insert(a).Append(i);
            int b=ta.mesh.elements(i)(jp1);
            TV_INT edge(min(a,b),max(a,b));
            if(edge_count.Get_Pointer(edge)==NULL) edge_count.Get_Or_Insert(edge)=1;
            else edge_count.Get_Or_Insert(edge)++;}}
    for(typename HASHTABLE<TV_INT,int>::ITERATOR it(edge_count);it.Valid();it.Next()){
        if(it.Data()==1) triedge2type.Get_Or_Insert(it.Key())=true;
        else{PHYSBAM_ASSERT(it.Data()==2);triedge2type.Get_Or_Insert(it.Key())=false;}}
    particle2type.Fill(false);
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(triedge2type);it.Valid();it.Next())
        if(it.Data()==true) for(int d=0;d<2;d++) particle2type(it.Key()(d))=true;
    for(int i=0;i<ta.particles.number;i++){
        if(particle2type(i)) boundary_particles.Append(i);
        else interior_particles.Append(i);}
    for(int i=0;i<interior_particles.m;i++){
        int p=interior_particles(i);
        ARRAY<int> unsorted_neighbors=particle2neighbortris.Get_Or_Insert(p);
        int original_size=unsorted_neighbors.m;
        ARRAY<int> sorted_neighbors;
        int selection=unsorted_neighbors(0);
        sorted_neighbors.Append(selection);
        int id=unsorted_neighbors.Find(selection);
        unsorted_neighbors.Remove_Index_Lazy(id);
        while(unsorted_neighbors.m!=0){
            for(int k=0;k<unsorted_neighbors.m;k++){
                ARRAY<int> common;
                common.Find_Common_Elements(ta.mesh.elements(selection),ta.mesh.elements(unsorted_neighbors(k)));
                PHYSBAM_ASSERT(common.m<=2); // sanity check
                if(common.m==2){
                    selection=unsorted_neighbors(k);
                    sorted_neighbors.Append(selection);
                    unsorted_neighbors.Remove_Index_Lazy(k);
                    break;}}}
        PHYSBAM_ASSERT(sorted_neighbors.m==original_size); // sanity check
        // sanity check
        for(int j=0;j<sorted_neighbors.m;j++){
            int jp1=(j==sorted_neighbors.m-1)?0:(j+1);
            ARRAY<int> common;
            common.Find_Common_Elements(ta.mesh.elements(sorted_neighbors(j)),ta.mesh.elements(sorted_neighbors(jp1)));
            PHYSBAM_ASSERT(common.m==2);}
        particle2neighbortris.Get_Or_Insert(p)=sorted_neighbors;}

    //*** go through all tris, each tri insert an interior node, build map node<->tri on the fly.
    HASHTABLE<int,int> tri2node;
    for(int i=0;i<ta.mesh.elements.m;i++){
        TV barycenter=TV();
        for(int d=0;d<3;d++) barycenter+=ta.particles.X(ta.mesh.elements(i)(d));
        barycenter*=(T)0.3333333333333333;
        X.Append(barycenter);
        type.Append(10);
        tri2node.Get_Or_Insert(i)=X.m-1;}

    //*** go through all triangle edges: if it's an interior edge, insert a face node; if it's a boundary edge, insert a boundary node. build map node<->edge on the fly.
    HASHTABLE<TV_INT,int> edge2node;
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(triedge2type);it.Valid();it.Next()){
        TV_INT edge=it.Key();
        X.Append((T)0.5*(ta.particles.X(edge(0))+ta.particles.X(edge(1))));
        edge2node.Get_Or_Insert(edge)=X.m-1;
        if(it.Data()) type.Append(1);
        else type.Append(100);}
    
    //*** go through all boundary partilces, each particle insert an boundary node, build map node<->particle on the fly.
    HASHTABLE<int,int> particle2node;
    for(int i=0;i<boundary_particles.m;i++){
        int p=boundary_particles(i);
        X.Append(ta.particles.X(p));
        type.Append(1);
        particle2node.Get_Or_Insert(p)=X.m-1;}
            
    //*** resize polygons to the size of particles.
    elements.Resize(ta.particles.number);

    //*** for all interior particle p
    for(int i=0;i<interior_particles.m;i++){
        int p=interior_particles(i);
        ARRAY<int>& nt=particle2neighbortris.Get_Or_Insert(p);
        for(int aa=0;aa<nt.m;aa++){
            int bb=(aa==nt.m-1)?0:(aa+1);
            int a=nt(aa),b=nt(bb);
            ARRAY<int> common;common.Find_Common_Elements(ta.mesh.elements(a),ta.mesh.elements(b));
            PHYSBAM_ASSERT(common.m==2);
            TV_INT edge(min(common(0),common(1)),max(common(0),common(1)));
            elements(p).Append(tri2node.Get_Or_Insert(a));
            elements(p).Append(edge2node.Get_Or_Insert(edge));}}
    
    //*** for all boundary particle p
    for(int i=0;i<boundary_particles.m;i++){
        int p=boundary_particles(i);
        ARRAY<int>& nt=particle2neighbortris.Get_Or_Insert(p);
        HASHTABLE<TV_INT,bool> boundary_edges_share_p;
        for(int j=0;j<nt.m;j++){
            int a=ta.mesh.elements(nt(j))(0),b=ta.mesh.elements(nt(j))(1),c=ta.mesh.elements(nt(j))(2);
            TV_INT edge1=TV_INT(min(a,b),max(a,b)),edge2=TV_INT(min(b,c),max(b,c)),edge3=TV_INT(min(c,a),max(c,a));
            if(triedge2type.Get_Or_Insert(edge1) && ((a==p) || (b==p))) boundary_edges_share_p.Get_Or_Insert(edge1)=true;
            if(triedge2type.Get_Or_Insert(edge2) && ((b==p) || (c==p))) boundary_edges_share_p.Get_Or_Insert(edge2)=true;
            if(triedge2type.Get_Or_Insert(edge3) && ((c==p) || (a==p))) boundary_edges_share_p.Get_Or_Insert(edge3)=true;}
        PHYSBAM_ASSERT(boundary_edges_share_p.Size()==2);
        TV_INT e1,e2;
        int count=0;
        for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(boundary_edges_share_p);it.Valid();it.Next()){
            if(count++==0) e1=it.Key();
            else e2=it.Key();}
        elements(p).Append(particle2node.Get_Or_Insert(p));
        elements(p).Append(edge2node.Get_Or_Insert(e1));
        
        ARRAY<int> ntlocal=nt;
        int selection;
        for(int j=0;j<ntlocal.m;j++){
            int a=ta.mesh.elements(ntlocal(j))(0),b=ta.mesh.elements(ntlocal(j))(1),c=ta.mesh.elements(ntlocal(j))(2);
            TV_INT edge1=TV_INT(min(a,b),max(a,b)),edge2=TV_INT(min(b,c),max(b,c)),edge3=TV_INT(min(c,a),max(c,a));
            if(edge1==e1 || edge2==e1 || edge3==e1){
                selection=ntlocal(j);
                break;}}
        elements(p).Append(tri2node.Get_Or_Insert(selection));
        
        int id=ntlocal.Find(selection);
        ntlocal.Remove_Index_Lazy(id);
        while(ntlocal.m>0){
            for(int k=0;k<ntlocal.m;k++){
                ARRAY<int> common;common.Find_Common_Elements(ta.mesh.elements(selection),ta.mesh.elements(ntlocal(k)));
                PHYSBAM_ASSERT(common.m==1 || common.m==2);
                if(common.m==2){
                    TV_INT edge(min(common(0),common(1)),max(common(0),common(1)));
                    selection=ntlocal(k);
                    elements(p).Append(edge2node.Get_Or_Insert(edge));
                    elements(p).Append(tri2node.Get_Or_Insert(selection));
                    ntlocal.Remove_Index_Lazy(k);
                    break;}}}
        PHYSBAM_ASSERT(ta.mesh.elements(selection)(0)==e2(0) || ta.mesh.elements(selection)(1)==e2(0) || ta.mesh.elements(selection)(2)==e2(0));
        PHYSBAM_ASSERT(ta.mesh.elements(selection)(0)==e2(1) || ta.mesh.elements(selection)(1)==e2(1) || ta.mesh.elements(selection)(2)==e2(1));
        elements(p).Append(edge2node.Get_Or_Insert(e2));}

    //*** insert face nodes to fix stuff
    HASHTABLE<TV_INT,int> face_fixer;
    for(int i=0;i<boundary_particles.m;i++){
        int p=boundary_particles(i);
        for(int aa=0;aa<elements(p).m;aa++){
            int bb=(aa==elements(p).m-1)?0:(aa+1);
            TV_INT edge(min(elements(p)(aa),elements(p)(bb)),max(elements(p)(aa),elements(p)(bb)));
            if((type(edge(0))==1 && type(edge(1))==1) || (type(edge(0))==1 && type(edge(1))==10) || (type(edge(0))==10 && type(edge(1))==1)){
                if(face_fixer.Get_Pointer(edge)==NULL){
                    X.Append((T)0.5*(X(edge(0))+X(edge(1))));
                    type.Append(100);
                    face_fixer.Get_Or_Insert(edge)=X.m-1;}}}}

    for(int i=0;i<boundary_particles.m;i++){
        int p=boundary_particles(i);
        ARRAY<int> new_elements=elements(p);
        for(int aa=0;aa<elements(p).m;aa++){
            int bb=(aa==elements(p).m-1)?0:(aa+1);
            TV_INT edge(min(elements(p)(aa),elements(p)(bb)),max(elements(p)(aa),elements(p)(bb)));
            if((type(edge(0))==1 && type(edge(1))==1) || (type(edge(0))==1 && type(edge(1))==10) || (type(edge(0))==10 && type(edge(1))==1)){
                int f_index=face_fixer.Get_Or_Insert(edge);
                int id=new_elements.Find(elements(p)(aa));
                new_elements.Insert(f_index,id+1);}}
        elements(p)=new_elements;}

    Xm=X;
    Initialize_Neighbor_Cells();
    Build_Association();
}

//#####################################################################
// Function Initialize_With_And_As_A_Triangulated_Area_And_Relocate_Particles_To_Tri_Centers
//#####################################################################
template<class T> void VORONOI_2D<T>::
Initialize_With_And_As_A_Triangulated_Area_And_Relocate_Particles_To_Tri_Centers(const TRIANGULATED_AREA<T>& ta,MPM_PARTICLES<TV>& mpm_particles)
{
    HASHTABLE<TV_INT,bool> triedge2type; // true-boundary,false-interior
    ARRAY<bool> particle2type(ta.particles.number); // true-boundary,false-interior
    HASHTABLE<TV_INT,int> edge_count;
    for(int i=0;i<ta.mesh.elements.m;i++){
        for(int j=0;j<3;j++){
            int jp1=(j==2)?0:(j+1);
            int a=ta.mesh.elements(i)(j);
            int b=ta.mesh.elements(i)(jp1);
            TV_INT edge(min(a,b),max(a,b));
            if(edge_count.Get_Pointer(edge)==NULL) edge_count.Get_Or_Insert(edge)=1;
            else edge_count.Get_Or_Insert(edge)++;}}
    for(typename HASHTABLE<TV_INT,int>::ITERATOR it(edge_count);it.Valid();it.Next()){
        if(it.Data()==1) triedge2type.Get_Or_Insert(it.Key())=true;
        else{PHYSBAM_ASSERT(it.Data()==2);triedge2type.Get_Or_Insert(it.Key())=false;}}
    particle2type.Fill(false);
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(triedge2type);it.Valid();it.Next())
        if(it.Data()==true) for(int d=0;d<2;d++) particle2type(it.Key()(d))=true;
    
    mpm_particles.Clean_Memory();
    mpm_particles.Resize(ta.mesh.elements.m);
    elements.Resize(ta.mesh.elements.m);
    for(int i=0;i<ta.mesh.elements.m;i++){
        TV barycenter;
        for(int d=0;d<3;d++){
            barycenter+=ta.particles.X(ta.mesh.elements(i)(d));
            elements(i).Append(ta.mesh.elements(i)(d));}
        barycenter*=(T)0.333333333333333;
        mpm_particles.X(i)=barycenter;}
    mpm_particles.Xm=mpm_particles.X;

    for(int i=0;i<ta.particles.X.m;i++){
        X.Append(ta.particles.X(i));
        if(particle2type(i)) type.Append(1);
        else type.Append(10);}
    
    HASHTABLE<TV_INT,int> face_fixer;
    for(int i=0;i<ta.mesh.elements.m;i++){
        for(int aa=0;aa<elements(i).m;aa++){
            int bb=(aa==elements(i).m-1)?0:(aa+1);
            TV_INT edge(min(elements(i)(aa),elements(i)(bb)),max(elements(i)(aa),elements(i)(bb)));
            if((type(edge(0))==1 && type(edge(1))==1) || (type(edge(0))==10 && type(edge(1))==10) || (type(edge(0))==1 && type(edge(1))==10) || (type(edge(0))==10 && type(edge(1))==1)){
                if(face_fixer.Get_Pointer(edge)==NULL){
                    X.Append((T)0.5*(X(edge(0))+X(edge(1))));
                    type.Append(100);
                    face_fixer.Get_Or_Insert(edge)=X.m-1;}}}}

    for(int i=0;i<ta.mesh.elements.m;i++){
        ARRAY<int> new_elements=elements(i);
        for(int aa=0;aa<elements(i).m;aa++){
            int bb=(aa==elements(i).m-1)?0:(aa+1);
            TV_INT edge(min(elements(i)(aa),elements(i)(bb)),max(elements(i)(aa),elements(i)(bb)));
            if((type(edge(0))==1 && type(edge(1))==1) || (type(edge(0))==10 && type(edge(1))==10) || (type(edge(0))==1 && type(edge(1))==10) || (type(edge(0))==10 && type(edge(1))==1)){
                int f_index=face_fixer.Get_Or_Insert(edge);
                int id=new_elements.Find(elements(i)(aa));
                new_elements.Insert(f_index,id+1);}}
        elements(i)=new_elements;}

    // for(int i=0;i<elements.m;i++){
    //     PHYSBAM_ASSERT(elements(i).m==6);
    //     LOG::cout<<type(elements(i)(0))<<" "<<type(elements(i)(1))<<" "<<type(elements(i)(2))<<" "<<type(elements(i)(3))<<" "<<type(elements(i)(4))<<" "<<type(elements(i)(5))<<std::endl;
    // }
    // exit(0);

    Xm=X;
    Initialize_Neighbor_Cells();
    Build_Association();

}

//#####################################################################
// function Initialize_Neighbor_Cells
//#####################################################################
template<class T> void VORONOI_2D<T>::
Initialize_Neighbor_Cells()
{
    LOG::cout<<"Initializing Neighbor Cells for voronoi mesh..."<<std::endl;
    TIMING_START;
    HASHTABLE<int,ARRAY<int> > f2e;
    for(int i=0;i<elements.m;i++) for(int j=0;j<elements(i).m;j++) if(type(elements(i)(j))==100) f2e.Get_Or_Insert(elements(i)(j)).Append(i);
    for(typename HASHTABLE<int,ARRAY<int> >::ITERATOR it(f2e);it.Valid();it.Next())
        if(it.Data().m>1){
            PHYSBAM_ASSERT(it.Data().m==2);
            neighbor_cells.Append(TRIPLE<int,int,bool>(it.Data()(0),it.Data()(1),true));}
    TIMING_END("");
    // LOG::cout<<neighbor_cells<<std::endl;

}
//#####################################################################
// Function Build_Association
//#####################################################################
template<class T> void VORONOI_2D<T>::
Build_Association()
{
    association.Clean_Memory();
    association.Resize(Xm.m);
    for(int c=0;c<elements.m;c++) for(int v=0;v<elements(c).m;v++) association(elements(c)(v)).Append(c);
}
//#####################################################################
// Function Build_Segments
//#####################################################################
template<class T> void VORONOI_2D<T>::
Build_Segments()
{
    segments.Clean_Memory();
    HASHTABLE<TV_INT,bool> seg;
    for(int e=0;e<elements.m;e++){
        int total_nodes=elements(e).m;
        for(int a=0;a<total_nodes;a++){
            int first=elements(e)(a),second;
            if(a==total_nodes-1) second=elements(e)(0);
            else second=elements(e)(a+1);
            if(first<second) seg.Get_Or_Insert(TV_INT(first,second))=true;
            else seg.Get_Or_Insert(TV_INT(second,first))=true;}}
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(seg);it.Valid();it.Next()) segments.Append(it.Key());
}
//#####################################################################
// Function Build_Boundary_Segments
//#####################################################################
template<class T> void VORONOI_2D<T>::
Build_Boundary_Segments()
{
    boundary_segments.Clean_Memory();
    HASHTABLE<TV_INT,int> seg;
    for(int e=0;e<elements.m;e++){
        int total_nodes=elements(e).m;
        for(int a=0;a<total_nodes;a++){
            int first=elements(e)(a),second;
            if(a==total_nodes-1) second=elements(e)(0);
            else second=elements(e)(a+1);
            if(first<second) seg.Get_Or_Insert(TV_INT(first,second))++;
            else seg.Get_Or_Insert(TV_INT(second,first))++;}}
    for(typename HASHTABLE<TV_INT,int>::ITERATOR it(seg);it.Valid();it.Next()) if(it.Data()==1) boundary_segments.Append(it.Key());
}

//#####################################################################
// Function Deform_Mesh_Using_Particle_Deformation
//#####################################################################
template<class T> void VORONOI_2D<T>::
Deform_Mesh_Using_Particle_Deformation(const ARRAY_VIEW<TV>& particle_Xm,const ARRAY_VIEW<TV>& particle_X,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fe,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fp,const bool constrain_face_centers)
{
    X.Clean_Memory();X.Resize(Xm.m);
    ARRAY<TV> particle_b(particle_X.m);
    ARRAY<MATRIX<T,TV::m> > particle_F(particle_X.m);
    for(int i=0;i<particle_b.m;i++){
        particle_F(i)=particle_Fe(i)*particle_Fp(i);
        particle_b(i)=particle_X(i)-particle_F(i)*particle_Xm(i);}
    for(int i=0;i<X.m;i++){
        X(i)=TV();
        for(int p=0;p<association(i).m;p++){
            int particle=association(i)(p);
            X(i)+=particle_F(particle)*Xm(i)+particle_b(particle);}
        X(i)/=association(i).m;}
    if(constrain_face_centers){
        HASHTABLE<int,TV_INT> fc2parents;
        for(int e=0;e<elements.m;e++){
            for(int a=0;a<elements(e).m;a++){
                int me=elements(e)(a);
                int left=(a==0)?elements(e)(elements(e).m-1):elements(e)(a-1);
                int right=(a==elements(e).m-1)?elements(e)(0):elements(e)(a+1);
                if(type(me)==100) fc2parents.Get_Or_Insert(me)=TV_INT(left,right).Sorted();}}
        for(typename HASHTABLE<int,TV_INT>::ITERATOR it(fc2parents);it.Valid();it.Next())
            X(it.Key())=(T)0.5*(X(it.Data()(0))+X(it.Data()(1)));}
}
//#####################################################################
// Function Crack
//#####################################################################
template<class T> void VORONOI_2D<T>::
Crack(const ARRAY_VIEW<TV>& particle_X,const T threshold)
{
    LOG::cout<<"crack start"<<std::endl;

    T threshold_squared=sqr(threshold);
    for(int nc=0;nc<neighbor_cells.m;nc++){
        TRIPLE<int,int,bool>& tr=neighbor_cells(nc);

        // if(tr.z){
        //     ARRAY<int> shared_nodes;shared_nodes.Find_Common_Elements(elements(tr.x),elements(tr.y));
        //     if(shared_nodes.m!=3){
        //         LOG::cout<<shared_nodes.m<<std::endl;
        //         LOG::cout<<tr.x<<" "<<elements(tr.x)<<std::endl;
        //         LOG::cout<<tr.y<<" "<<elements(tr.y)<<std::endl;
        //         LOG::cout<<neighbor_cells(nc)<<std::endl;
        //     }
        //     PHYSBAM_ASSERT(shared_nodes.m==3);
        // }

        if(tr.z && (particle_X(tr.x)-particle_X(tr.y)).Magnitude_Squared()>threshold_squared){
            ARRAY<int> shared_nodes;shared_nodes.Find_Common_Elements(elements(tr.x),elements(tr.y));

            // if(shared_nodes.m!=3){
            //     LOG::cout<<shared_nodes.m<<std::endl;
            //     LOG::cout<<tr.x<<" "<<elements(tr.x)<<std::endl;
            //     LOG::cout<<tr.y<<" "<<elements(tr.y)<<std::endl;
            //     LOG::cout<<neighbor_cells(nc)<<std::endl;
            // }
            PHYSBAM_ASSERT(shared_nodes.m==3);

            VECTOR<int,3> types(type(shared_nodes(0)),type(shared_nodes(1)),type(shared_nodes(2)));
            int sum=types.Sum();
            switch(sum){
                case 120:{
                    // LOG::cout<<"case 120"<<std::endl;
                    int f1=0,i1=0,i2=0;
                    for(int i=0;i<3;i++){
                        if(types(i)==100) f1=shared_nodes(i);
                        else if(types(i)==10){
                            if(i1==0) i1=shared_nodes(i);
                            else i2=shared_nodes(i);}
                        else PHYSBAM_FATAL_ERROR();}
                    type(i1)=1;type(i2)=1;
                    type.Append(100);Xm.Append(Xm(f1));
                    int f2=Xm.m-1;
                    for(int i=0;i<elements(tr.y).m;i++) if(elements(tr.y)(i)==f1) elements(tr.y)(i)=f2;
                    tr.z=false;


                    // LOG::cout<<elements(6)<<std::endl;
                    // LOG::cout<<elements(11)<<std::endl;
                    // ARRAY<int> sn;sn.Find_Common_Elements(elements(6),elements(11));
                    // LOG::cout<<sn.m<<std::endl;

                    break;}
                case 111:{
                    // LOG::cout<<"case 111"<<std::endl;
                    int f1=0,i1=0,b1=0;
                    for(int i=0;i<3;i++){
                        if(types(i)==100) f1=shared_nodes(i);
                        else if(types(i)==10) i1=shared_nodes(i);
                        else if(types(i)==1) b1=shared_nodes(i);
                        else PHYSBAM_FATAL_ERROR();}
                    type(i1)=1;
                    type.Append(100);Xm.Append(Xm(f1));
                    int f2=Xm.m-1;
                    for(int i=0;i<elements(tr.y).m;i++) if(elements(tr.y)(i)==f1) elements(tr.y)(i)=f2;
                    tr.z=false;

                    type.Append(1);Xm.Append(Xm(b1));
                    int b2=Xm.m-1;

                    ARRAY<int> csb;
                    for(int i=0;i<elements.m;i++) if(elements(i).Contains(b1)) csb.Append(i);
                    ARRAY<int> groupA;
                    groupA.Append(csb(0));
                    csb.Remove_Index_Lazy(0);
                    for(int i=0;i<groupA.m;i++){
                        ARRAY<int> new_members;
                        for(int j=0;j<csb.m;j++){
                            int small_cell=(groupA(i)>csb(j))?csb(j):groupA(i);
                            int big_cell=(groupA(i)==small_cell)?csb(j):groupA(i);
                            if(neighbor_cells.Contains(TRIPLE<int,int,bool>(small_cell,big_cell,true)))
                                new_members.Append(csb(j));}
                        for(int j=0;j<new_members.m;j++){
                            groupA.Append(new_members(j));
                            csb.Remove_Index_Lazy(csb.Find(new_members(j)));}}

                    for(int i=0;i<groupA.m;i++){
                        int this_cell=groupA(i);
                        for(int j=0;j<elements(this_cell).m;j++) if(elements(this_cell)(j)==b1) elements(this_cell)(j)=b2;}

                    // LOG::cout<<elements(6)<<std::endl;
                    // LOG::cout<<elements(11)<<std::endl;
                    // ARRAY<int> sn;sn.Find_Common_Elements(elements(6),elements(11));
                    // LOG::cout<<sn.m<<std::endl;


                    break;}
                case 102:{
                    // LOG::cout<<"case 102"<<std::endl;

                    // {
                    //     LOG::cout<<elements(6)<<std::endl;
                    //     LOG::cout<<elements(11)<<std::endl;
                    //     ARRAY<int> sn;sn.Find_Common_Elements(elements(6),elements(11));
                    //     LOG::cout<<sn.m<<std::endl;
                    // }

                    int f1=0,bfirst1=0,bsecond1=0;
                    for(int i=0;i<3;i++){
                        if(types(i)==100) f1=shared_nodes(i);
                        else if(types(i)==1){
                            if(bfirst1==0) bfirst1=shared_nodes(i);
                            else bsecond1=shared_nodes(i);}
                        else PHYSBAM_FATAL_ERROR();}
                    type.Append(100);Xm.Append(Xm(f1));
                    int f2=Xm.m-1;
                    for(int i=0;i<elements(tr.y).m;i++) if(elements(tr.y)(i)==f1) elements(tr.y)(i)=f2;
                    tr.z=false;

                    // LOG::cout<<"CRACKING:"<<tr<<std::endl;
                    // {
                    //     LOG::cout<<elements(6)<<std::endl;
                    //     LOG::cout<<elements(11)<<std::endl;
                    //     ARRAY<int> sn;sn.Find_Common_Elements(elements(6),elements(11));
                    //     LOG::cout<<sn.m<<std::endl;
                    // }                    

                    type.Append(1);Xm.Append(Xm(bfirst1));
                    int bfirst2=Xm.m-1;

                    ARRAY<int> csbfirst;
                    for(int i=0;i<elements.m;i++) if(elements(i).Contains(bfirst1)) csbfirst.Append(i);
                    ARRAY<int> groupA;
                    groupA.Append(csbfirst(0));
                    csbfirst.Remove_Index_Lazy(0);
                    for(int i=0;i<groupA.m;i++){
                        ARRAY<int> new_members;
                        for(int j=0;j<csbfirst.m;j++){
                            int small_cell=(groupA(i)>csbfirst(j))?csbfirst(j):groupA(i);
                            int big_cell=(groupA(i)==small_cell)?csbfirst(j):groupA(i);
                            if(neighbor_cells.Contains(TRIPLE<int,int,bool>(small_cell,big_cell,true)))
                                new_members.Append(csbfirst(j));}
                        for(int j=0;j<new_members.m;j++){
                            groupA.Append(new_members(j));
                            csbfirst.Remove_Index_Lazy(csbfirst.Find(new_members(j)));}}

                    for(int i=0;i<groupA.m;i++){
                        int this_cell=groupA(i);
                        for(int j=0;j<elements(this_cell).m;j++) if(elements(this_cell)(j)==bfirst1) elements(this_cell)(j)=bfirst2;}

                    // {
                    //     LOG::cout<<elements(6)<<std::endl;
                    //     LOG::cout<<elements(11)<<std::endl;
                    //     ARRAY<int> sn;sn.Find_Common_Elements(elements(6),elements(11));
                    //     LOG::cout<<sn.m<<std::endl;
                    // }                    

                    type.Append(1);Xm.Append(Xm(bsecond1));
                    int bsecond2=Xm.m-1;

                    ARRAY<int> csbsecond;
                    for(int i=0;i<elements.m;i++) if(elements(i).Contains(bsecond1)) csbsecond.Append(i);
                    ARRAY<int> groupC;
                    groupC.Append(csbsecond(0));
                    csbsecond.Remove_Index_Lazy(0);
                    for(int i=0;i<groupC.m;i++){
                        ARRAY<int> new_members;
                        for(int j=0;j<csbsecond.m;j++){
                            int small_cell=(groupC(i)>csbsecond(j))?csbsecond(j):groupC(i);
                            int big_cell=(groupC(i)==small_cell)?csbsecond(j):groupC(i);
                            if(neighbor_cells.Contains(TRIPLE<int,int,bool>(small_cell,big_cell,true)))
                                new_members.Append(csbsecond(j));}
                        for(int j=0;j<new_members.m;j++){
                            groupC.Append(new_members(j));
                            csbsecond.Remove_Index_Lazy(csbsecond.Find(new_members(j)));}}

                    // LOG::cout<<"groupC:"<<groupC<<std::endl;
                    // LOG::cout<<"are 6 and 1  neighbors?"<<neighbor_cells<<std::endl;
                    for(int i=0;i<groupC.m;i++){
                        int this_cell=groupC(i);
                        for(int j=0;j<elements(this_cell).m;j++) if(elements(this_cell)(j)==bsecond1) elements(this_cell)(j)=bsecond2;}

                    // LOG::cout<<elements(6)<<std::endl;
                    // LOG::cout<<elements(11)<<std::endl;
                    // ARRAY<int> sn;sn.Find_Common_Elements(elements(6),elements(11));
                    // LOG::cout<<sn.m<<std::endl;
                    


                    break;}    
                default:
                    PHYSBAM_FATAL_ERROR();
            };
        }
    }
}
//#####################################################################
template class VORONOI_2D<float>;
template class VORONOI_2D<double>;
}
