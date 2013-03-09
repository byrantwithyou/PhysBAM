//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
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
    // Preprocessing
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
            

}
//#####################################################################
// Function Initialize_Neighbor_Cells
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
Deform_Mesh_Using_Particle_Deformation(const ARRAY_VIEW<TV>& particle_Xm,const ARRAY_VIEW<TV>& particle_X,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fe,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fp)
{
    X.Clean_Memory();X.Resize(Xm.m);
    bool use_deformation_gradient=true;
    bool use_average_world_space=false;
    if(use_deformation_gradient){
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
            X(i)/=association(i).m;}}
    else if(use_average_world_space){
        for(int i=0;i<X.m;i++){
            X(i)=TV();
            for(int p=0;p<association(i).m;p++){
                int particle=association(i)(p);
                X(i)+=particle_X(particle);}
            X(i)/=association(i).m;}}
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
        if(tr.z && (particle_X(tr.x)-particle_X(tr.y)).Magnitude_Squared()>threshold_squared){
            ARRAY<int> shared_nodes;shared_nodes.Find_Common_Elements(elements(tr.x),elements(tr.y));
            PHYSBAM_ASSERT(shared_nodes.m==3);
            VECTOR<int,3> types(type(shared_nodes(0)),type(shared_nodes(1)),type(shared_nodes(2)));
            int sum=types.Sum();
            switch(sum){
                case 120:{
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
                    break;}
                case 111:{
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
                    ARRAY<int> cells_share_b1;
                    for(int i=0;i<elements.m;i++) if(elements(i).Contains(b1)) cells_share_b1.Append(i);
                    ARRAY<int> groupA,groupB;
                    groupA.Append(cells_share_b1(0));
                    for(int i=1;i<cells_share_b1.m;i++){
                        int this_cell=cells_share_b1(i);
                        bool flag=false;
                        for(int j=0;j<groupA.m;j++){
                            int that_cell=groupA(j);
                            int small_cell=(this_cell>that_cell)?that_cell:this_cell;
                            int big_cell=(this_cell==small_cell)?that_cell:this_cell;
                            if(neighbor_cells.Contains(TRIPLE<int,int,bool>(small_cell,big_cell,true))){
                                flag=true;
                                break;}}
                        if(flag) groupA.Append(this_cell);
                        else groupB.Append(this_cell);}
                    for(int i=0;i<groupB.m;i++){
                        int this_cell=groupB(i);
                        for(int j=0;j<elements(this_cell).m;j++) if(elements(this_cell)(j)==b1) elements(this_cell)(j)=b2;}
                    break;}
                case 102:{
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
                    
                    type.Append(1);Xm.Append(Xm(bfirst1));
                    int bfirst2=Xm.m-1;
                    ARRAY<int> cells_share_bfirst1;
                    for(int i=0;i<elements.m;i++) if(elements(i).Contains(bfirst1)) cells_share_bfirst1.Append(i);
                    ARRAY<int> groupA,groupB;
                    groupA.Append(cells_share_bfirst1(0));
                    for(int i=1;i<cells_share_bfirst1.m;i++){
                        int this_cell=cells_share_bfirst1(i);
                        bool flag=false;
                        for(int j=0;j<groupA.m;j++){
                            int that_cell=groupA(j);
                            int small_cell=(this_cell>that_cell)?that_cell:this_cell;
                            int big_cell=(this_cell==small_cell)?that_cell:this_cell;
                            if(neighbor_cells.Contains(TRIPLE<int,int,bool>(small_cell,big_cell,true))){
                                flag=true;
                                break;}}
                        if(flag) groupA.Append(this_cell);
                        else groupB.Append(this_cell);}
                    for(int i=0;i<groupB.m;i++){
                        int this_cell=groupB(i);
                        for(int j=0;j<elements(this_cell).m;j++) if(elements(this_cell)(j)==bfirst1) elements(this_cell)(j)=bfirst2;}

                    type.Append(1);Xm.Append(Xm(bsecond1));
                    int bsecond2=Xm.m-1;
                    ARRAY<int> cells_share_bsecond1;
                    for(int i=0;i<elements.m;i++) if(elements(i).Contains(bsecond1)) cells_share_bsecond1.Append(i);
                    ARRAY<int> groupC,groupD;
                    groupC.Append(cells_share_bsecond1(0));
                    for(int i=1;i<cells_share_bsecond1.m;i++){
                        int this_cell=cells_share_bsecond1(i);
                        bool flag=false;
                        for(int j=0;j<groupC.m;j++){
                            int that_cell=groupC(j);
                            int small_cell=(this_cell>that_cell)?that_cell:this_cell;
                            int big_cell=(this_cell==small_cell)?that_cell:this_cell;
                            if(neighbor_cells.Contains(TRIPLE<int,int,bool>(small_cell,big_cell,true))){
                                flag=true;
                                break;}}
                        if(flag) groupC.Append(this_cell);
                        else groupD.Append(this_cell);}
                    for(int i=0;i<groupD.m;i++){
                        int this_cell=groupD(i);
                        for(int j=0;j<elements(this_cell).m;j++) if(elements(this_cell)(j)==bsecond1) elements(this_cell)(j)=bsecond2;}

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
