//#####################################################################
// Copyright 2016, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Intersections/BOX_POLYGON_INTERSECTION.h>
#include "VORONOI_DIAGRAM.h"
using namespace PhysBAM;
template<class V,class F>
inline void For_Each_Outgoing_Coedge(V* v,F func)
{
    auto* a=v->coedge,*b=a->pair->next,*c=b->pair->next;
    func(a);
    func(b);
    func(c);
}
template<class V,class F>
inline void For_Each_Outgoing_Coedge_Reverse(V* v,F func)
{
    auto* a=v->coedge,*b=a->pair->next,*c=b->pair->next;
    func(c);
    func(b);
    func(a);
}
//#####################################################################
// Function Criterion
//#####################################################################
template<class T> T VORONOI_DIAGRAM<T>::VERTEX::
Criterion(const TV& A,const TV& B,const TV& C,const TV& D) const
{
    TV X=A-D;
    TV Y=B-D;
    TV Z=C-D;
    T ax=X.Magnitude_Squared()*Y.Cross(Z).x;
    T ay=Y.Magnitude_Squared()*Z.Cross(X).x;
    T az=Z.Magnitude_Squared()*X.Cross(Y).x;
    return ax+ay+az;
}
//#####################################################################
// Function Criterion
//#####################################################################
template<class T> T VORONOI_DIAGRAM<T>::VERTEX::
Criterion(const TV& X) const
{
    COEDGE* ce=coedge;
    TV A=ce->cell->X;
    ce=ce->pair->next;
    TV B=ce->cell->X;
    ce=ce->pair->next;
    TV C=ce->cell->X;
    return Criterion(A,B,C,X);
}
//#####################################################################
// Function Compute_Point
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::VERTEX::
Compute_Point(const TV& A,const TV& B,const TV& C)
{
    TV Y=B-A,Z=C-A;
    T b=Y.Magnitude_Squared();
    T c=Z.Magnitude_Squared();
    X=TV(b*Z.y-Y.y*c,Y.x*c-b*Z.x)/(2*Y.Cross(Z).x)+A;
}
//#####################################################################
// Function Compute_Point
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::VERTEX::
Compute_Point()
{
    COEDGE* ce=coedge;
    TV A=ce->cell->X;
    ce=ce->pair->next;
    TV B=ce->cell->X;
    ce=ce->pair->next;
    TV C=ce->cell->X;
    Compute_Point(A,B,C);
}
//#####################################################################
// Function Select_First_In_Vertex
//#####################################################################
template<class CE,class TV> inline auto
Select_First_In_Vertex(CE* start,const TV& new_pt)
{
    if(!start->head) return start->tail;
    if(!start->tail) return start->head;
    auto crit0=start->head->Criterion(new_pt);
    auto crit1=start->tail->Criterion(new_pt);
    if(crit0<crit1) return start->head;
    return start->tail;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> VORONOI_DIAGRAM<T>::
VORONOI_DIAGRAM()
    :radius(0)
{
    random.Set_Seed(1235);
}
//#####################################################################
// Function Discover_Inside
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Discover_Inside(ARRAY<COEDGE*>& in_ce,ARRAY<COEDGE*>& adj_ce,
    ARRAY<VERTEX*>& in_v,COEDGE* ce,const TV& new_pt)
{
    if(!ce->head || ce->head->state!=unknown){
        adj_ce.Append(ce);
        return;}

    VERTEX* v = ce->head;
    COEDGE* ceL = ce->next;
    COEDGE* ceR = ceL->pair->next;

    if((ceL->head && ceL->head->state==in) // (C4)
        || (ceR->head && ceR->head->state==in) // (C4)
        || ceR->cell->state==incident // (C5)
        || v->Criterion(new_pt)>=0){ // geometric criterion
        v->state=out;
        adj_ce.Append(ce);
        return;}

    v->state=in;
    ceR->cell->state=incident;
    in_v.Append(v);
    in_ce.Append(ce);
    Discover_Inside(in_ce,adj_ce,in_v,ceR,new_pt);
    Discover_Inside(in_ce,adj_ce,in_v,ceL,new_pt);
}
//#####################################################################
// Function Insert_Point
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Point(int p,const TV& new_pt)
{
    if(p>=clipped_piece_offset)
        Insert_Point(clipped_pieces(p-clipped_piece_offset).coedge,new_pt);
    else Insert_Point(pieces(p).coedge,new_pt);
}
//#####################################################################
// Function Insert_Point
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Point(COEDGE* start,const TV& new_pt)
{
    ARRAY<VERTEX*> in_vertices;
    ARRAY<COEDGE*> in_coedges,adj_coedges,new_coedges;

    // Make sure at least one vertex is in in_vertices
    VERTEX* v=Select_First_In_Vertex(start,new_pt);
    in_vertices.Append(v);
    v->state=in;
    For_Each_Outgoing_Coedge(v,[](COEDGE* ce){ce->cell->state=incident;});
    For_Each_Outgoing_Coedge_Reverse(v,[&](COEDGE* ce){Discover_Inside(in_coedges,adj_coedges,in_vertices,ce,new_pt);});

    // Allocate new cell
    CELL* cell = new CELL;
    cells.Append(cell);
    cell->X = new_pt;
    
    adj_coedges.Append(+adj_coedges(0)); // + to avoid aliasing
    for(int i=0;i<adj_coedges.m-1;i++){
        COEDGE* ce=adj_coedges(i)->pair;
        if(ce->tail) ce->tail->state=unknown;
        ce->cell->coedge=ce;

        VERTEX* H=new VERTEX;
        ce->head=H;
        ce->pair->tail=H;
        H->state=unknown;
        H->coedge=ce->pair;
        H->Compute_Point(new_pt,ce->cell->X,ce->pair->cell->X);
        
        COEDGE* A=new COEDGE;
        new_coedges.Append(A);
        ce->next=A;
        A->prev=ce;
        A->next=adj_coedges(i+1);
        A->next->prev=A;
        A->cell=ce->cell;
        A->cell->coedge=A;
        A->tail=H;
        ce->cell->state=non_incident;

        COEDGE* B=new COEDGE;
        A->pair=B;
        B->pair=A;
        B->head=H;
        B->cell=cell;}

    for(int i=0;i<in_coedges.m;i++){
        Remove_Coedge(in_coedges(i));
        Remove_Coedge(in_coedges(i)->pair);
        delete in_coedges(i)->pair;
        delete in_coedges(i);}

    for(int i=0;i<adj_coedges.m-1;i++){
        COEDGE* D=adj_coedges(i+1)->pair->next;
        COEDGE* E=adj_coedges(i)->pair->next->pair;
        D->next=adj_coedges(i);
        D->next->prev=D;
        D->head=D->next->tail;
        E->next=D->pair;
        D->pair->prev=E;
        D->pair->tail=E->head;}

    cell->coedge=adj_coedges(0)->pair->next->pair;
    
    for(int i=0;i<adj_coedges.m-1;i++){
        Update_Coedge(adj_coedges(i));
        Update_Coedge(adj_coedges(i)->pair);}

    for(int i=0;i<new_coedges.m;i++){
        Insert_Coedge(new_coedges(i));
        Insert_Coedge(new_coedges(i)->pair);}

    for(int i=0;i<in_vertices.m;i++){
        delete in_vertices(i);}
}
//#####################################################################
// Function Visualize_State
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Visualize_State(const char* title) const
{
    VECTOR<T,3> color_r[2]={{1,1,0},{0,1,0}};
    VECTOR<T,3> color_rc(1,.5,0);
    VECTOR<T,3> color_c[3]={{1,1,0},{1,0,0},{0,1,0}};
    VECTOR<T,3> color_cv(0,1,0);
//    VECTOR<T,3> color_v[3]={{0,0,1},{1,0,0},{0,1,0}};
    T unbounded_length = .5;
//    T offset = 0;

//    VECTOR<T,3> color_p[]={{1,.5,.5},{.5,1,.5},{.5,.5,1},{0,1,0},{1,.5,0},{1,1,0},{0,.5,0}};
    
    // for(auto it:coedges)
    // {
    //     const VERTEX* h=it->head;
    //     const VERTEX* t=it->tail;
    //     assert(h || t);
    //     TV XH,XT;
    //     if(h) XH=h->X;
    //     else{
    //         TV dir = (it->cell->X - it->pair->cell->X).Rotate_Clockwise_90().Normalized();
    //         XH=XT=t->X+dir*unbounded_length;}

    //     if(t) XT=t->X;
    //     else{
    //         TV dir = (it->cell->X - it->pair->cell->X).Rotate_Clockwise_90().Normalized();
    //         XT=XH-dir*unbounded_length;}

    //     TV dir=(XH-XT).Normalized()*offset,dir_in=dir.Rotate_Counterclockwise_90();

    //     Add_Debug_Object(VECTOR<TV,2>(XH-dir,XT+dir),color_c[0]);}

    for(auto it:cells)
    {
        Add_Debug_Particle(it->X,color_r[it->state]);
        COEDGE* ce=it->coedge;
        PHYSBAM_ASSERT(ce);
        while(ce)
        {
            const VERTEX* h=ce->head;
            const VERTEX* t=ce->tail;
            assert(h || t);
            TV XH,XT;
            if(h) XH=h->X;
            else{
                TV dir = (ce->cell->X - ce->pair->cell->X).Rotate_Clockwise_90().Normalized();
                XH=XT=t->X+dir*unbounded_length;}

            if(t) XT=t->X;
            else{
                TV dir = (ce->cell->X - ce->pair->cell->X).Rotate_Clockwise_90().Normalized();
                XT=XH-dir*unbounded_length;}

            Add_Debug_Object(VECTOR<TV,2>(XH,XT),color_c[0]);

            ce=ce->next;
            if(ce==it->coedge) break;
        }
        if(!ce)
        {
            COEDGE* ce=it->coedge->prev;
            while(ce)
            {
                const VERTEX* h=ce->head;
                const VERTEX* t=ce->tail;
                assert(h || t);
                TV XH,XT;
                if(h) XH=h->X;
                else{
                    TV dir = (ce->cell->X - ce->pair->cell->X).Rotate_Clockwise_90().Normalized();
                    XH=XT=t->X+dir*unbounded_length;}

                if(t) XT=t->X;
                else{
                    TV dir = (ce->cell->X - ce->pair->cell->X).Rotate_Clockwise_90().Normalized();
                    XT=XH-dir*unbounded_length;}

                Add_Debug_Object(VECTOR<TV,2>(XH,XT),color_c[0]);

                ce=ce->prev;
            }
        }
    }

    // for(auto it:vertices)
    //     Add_Debug_Particle(it->X,color_v[it->state]);

    Flush_Frame<TV>(title);
}
//#####################################################################
// Function Sanity_Checks
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Sanity_Checks() const
{
    for(int i=0;i<cells.m;i++){
        COEDGE* ce=cells(i)->coedge;
        do{
            assert(ce->cell==cells(i));
            ce=ce->next;
        }while(ce && ce!=cells(i)->coedge);}

    for(int i=0;i<pieces.m;i++) assert(pieces(i).coedge->piece==i);
    for(int i=0;i<clipped_pieces.m;i++) assert(clipped_pieces(i).coedge->piece==i+clipped_piece_offset);
}
//#####################################################################
// Function Update_Subtree_Area_Safe
//#####################################################################
template<class A> inline void
Update_Subtree_Area_Safe(ARRAY<A>& pieces,int i)
{
    int j=2*i+1,k=j+1;
    auto area=pieces(i).area;
    if(j<pieces.m){
        area+=pieces(j).subtree_area;
        if(k<pieces.m)
            area+=pieces(k).subtree_area;}
    pieces(i).subtree_area=area;
}
//#####################################################################
// Function Update_Subtree_Area_Safe
//#####################################################################
template<class A> inline void
Update_Subtree_Area(ARRAY<A>& pieces,int i)
{
    int j=2*i+1,k=j+1;
    pieces(i).subtree_area=pieces(i).area+pieces(j).subtree_area+pieces(k).subtree_area;
}
//#####################################################################
// Function Update_Piece_Tree
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Update_Piece_Tree(int i)
{
    Update_Subtree_Area_Safe(pieces,i);
    if(i>0){
        i=(i-1)/2;
        Update_Subtree_Area_Safe(pieces,i);
        while(i>0){
            i=(i-1)/2;
            Update_Subtree_Area(pieces,i);}}
}
//#####################################################################
// Function Update_Clipped_Piece_Tree
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Update_Clipped_Piece_Tree(int i)
{
    Update_Subtree_Area_Safe(clipped_pieces,i);
    if(i>0){
        i=(i-1)/2;
        Update_Subtree_Area_Safe(clipped_pieces,i);
        while(i>0){
            i=(i-1)/2;
            Update_Subtree_Area(clipped_pieces,i);}}
}
//#####################################################################
// Function Insert_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Coedge(COEDGE* ce)
{
    if(ce->cell->type==type_outside) return;
    if(!ce->head || !ce->tail){ // Infinite edges should not intersect the bounding box.
        assert(ce->cell->outside);
        assert(ce->head || ce->tail);
        if(ce->tail) assert(!bounding_box.Lazy_Inside(ce->tail->X));
        if(ce->head) assert(!bounding_box.Lazy_Inside(ce->head->X));}
    else if(!bounding_box.Lazy_Inside(ce->head->X) || !bounding_box.Lazy_Inside(ce->tail->X))
        Insert_Clipped_Coedge(ce);
    else{
        PIECE p;
        p.coedge=ce;
        p.h.A=ce->cell->X;
        p.h.B=ce->tail->X;
        p.h.C=ce->head->X;
        p.area=p.h.Compute(ce,radius,false);
        if(p.area<=0) ce->piece=-1;
        else{
            int i=pieces.Append(p);
            Update_Piece_Tree(i);
            ce->piece=i;}}
}
//#####################################################################
// Function Insert_Clipped_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Clipped_Coedge(COEDGE* ce)
{
    assert(bounding_box.Lazy_Inside(ce->cell->X));

    ARRAY<TV> polygon;
    polygon.Append(ce->cell->X);
    polygon.Append(ce->tail->X);
    polygon.Append(ce->head->X);
    INTERSECTION::Box_Polygon_Intersection(bounding_box,polygon,false);
    
    assert(polygon.m<=6);
    CLIPPED_PIECE p;
    p.coedge=ce;
    p.num_sub_pieces=polygon.m-2;
    T area=0;
    for(int i=0;i<polygon.m-2;i++){
        p.sub_pieces[i].A=polygon(0);
        p.sub_pieces[i].B=polygon(i+1);
        p.sub_pieces[i].C=polygon(i+2);
        area+=p.sub_pieces[i].Compute(ce,radius,true);}

    p.area=area;
    if(p.area<=0) ce->piece=-1;
    else{
        int i=clipped_pieces.Append(p);
        Update_Clipped_Piece_Tree(i);
        ce->piece=clipped_piece_offset+i;}
}
//#####################################################################
// Function Remove_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Remove_Coedge(COEDGE* ce)
{
    assert((ce->cell->type==type_outside)==(ce->piece<0));
    if(ce->piece>=clipped_piece_offset)
        Remove_Clipped_Piece(ce->piece-clipped_piece_offset);
    else if(ce->piece>=0) Remove_Piece(ce->piece);
}
//#####################################################################
// Function Remove_Piece
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Remove_Piece(int p)
{
    PIECE& last=pieces.Last(), &rem=pieces(p);
    rem.coedge->piece=-1;
    last.coedge->piece=p;
    rem=last;
    last.area=0;
    Update_Piece_Tree(p);
    Update_Piece_Tree(pieces.m-1);
    pieces.Pop();
}
//#####################################################################
// Function Remove_Piece
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Remove_Clipped_Piece(int p)
{
    CLIPPED_PIECE& last=clipped_pieces.Last(), &rem=clipped_pieces(p);
    rem.coedge->piece=-1;
    last.coedge->piece=p+clipped_piece_offset;
    rem=last;
    last.area=0;
    Update_Clipped_Piece_Tree(p);
    Update_Clipped_Piece_Tree(clipped_pieces.m-1);
    clipped_pieces.Pop();
}
//#####################################################################
// Function Update_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Update_Coedge(COEDGE* ce)
{
    Remove_Coedge(ce);
    Insert_Coedge(ce);
}
//#####################################################################
// Function Choose_Piece
//#####################################################################
template<class T> int VORONOI_DIAGRAM<T>::
Choose_Piece() const
{
    T pieces_area=pieces.m?pieces(0).subtree_area:0;
    T clipped_pieces_area=clipped_pieces.m?clipped_pieces(0).subtree_area:0;
    T total_area=pieces_area+clipped_pieces_area;
    assert(total_area>=0);
    if(total_area==0) return -1;
    T r=random.Get_Uniform_Number(0,total_area);
    int i=0;
    if(r>pieces_area){
        r-=pieces_area;
        while(1){
            assert(r<=clipped_pieces(i).subtree_area+1e-5);
            if(r<=clipped_pieces(i).area) return clipped_piece_offset+i;
            int j=2*i+1,k=j+1;
            r-=clipped_pieces(i).area;
            if(j>=clipped_pieces.m){assert(r<1e-5);return clipped_piece_offset+i;}
            if(r<=clipped_pieces(j).subtree_area){i=j;continue;}
            r-=clipped_pieces(j).subtree_area;
            if(k>=clipped_pieces.m){assert(r<1e-5);return clipped_piece_offset+j;}
            i=k;}}
    else{
        while(1){
            assert(r<=pieces(i).subtree_area+1e-5);
            if(r<=pieces(i).area) return i;
            int j=2*i+1,k=j+1;
            r-=pieces(i).area;
            if(j>=pieces.m){assert(r<1e-5);return i;}
            if(r<=pieces(j).subtree_area){i=j;continue;}
            r-=pieces(j).subtree_area;
            if(k>=pieces.m){assert(r<1e-5);return j;}
            i=k;}}
}
//#####################################################################
// Function Choose_Feasible_Point
//#####################################################################
template<class T> auto VORONOI_DIAGRAM<T>::
Choose_Feasible_Point(int p) const -> TV
{
    if(p>=clipped_piece_offset){
        const CLIPPED_PIECE& cp=clipped_pieces(p-clipped_piece_offset);
        T r=random.Get_Uniform_Number(0,cp.area);
        int s=0;
        for(;s<cp.num_sub_pieces-1;s++){
            r-=cp.sub_pieces[s].area;
            if(r<=0) break;}
        return cp.sub_pieces[s].Choose_Feasible_Point(random,radius);}
    else return pieces(p).h.Choose_Feasible_Point(random,radius);
}
//#####################################################################
// Function Disk_Inside_Triangle
//#####################################################################
template<class TV,class T> inline T
Disk_Inside_Triangle(TV B,TV C,T radius)
{
    return (T).5*(B.Cross(C).x-TV::Angle_Between(B,C)*sqr(radius));
}
//#####################################################################
// Function Compute
//#####################################################################
template<class T> T VORONOI_DIAGRAM<T>::PIECE_HELPER::
Compute(COEDGE* ce,T radius,bool clipped)
{
    B-=A;
    C-=A;

    if(ce->cell->outside){
        type=no_disc;
        return area=(T).5*B.Cross(C).x;} // Just a triangle: case A

    T B_mag2_min_r2=B.Magnitude_Squared()-sqr(radius);
    T C_mag2_min_r2=C.Magnitude_Squared()-sqr(radius);
    bool B_in_disk=B_mag2_min_r2<=0;
    bool C_in_disk=C_mag2_min_r2<=0;
    if(B_in_disk && C_in_disk){
        type=empty;
        return area=0;} // entire triangle in circle: case B

    QUADRATIC<T> quad((B-C).Magnitude_Squared(),-2*B.Dot(B-C),B_mag2_min_r2);
    quad.Compute_Roots();
    if(quad.roots==1){quad.roots=2;quad.root2=quad.root1;}
    if(quad.roots==0 || quad.root1>=1 || quad.root2<=0){
        type=full_disc;
        return area=Disk_Inside_Triangle(B,C,radius);} // case C

    if(quad.root1<=0 && quad.root2>=1){
        type=empty;
        return area=0;} // case B (numerical error)

    if(quad.root1>0 && quad.root2<1){ // case E
        type=both_out;
        TV P=B+quad.root1*(C-B),Q=B+quad.root2*(C-B);
        aux0=quad.root1;
        aux1=quad.root2;
        T A0=Disk_Inside_Triangle(B,P,radius);
        T A1=Disk_Inside_Triangle(Q,C,radius);
        aux2=A0/(A0+A1);
        return area=A0+A1;}

    aux0=quad.root1>0?quad.root1:quad.root2;
    TV P=B+aux0*(C-B);
    if(B_in_disk){
        type=out1;
        return area=Disk_Inside_Triangle(P,C,radius);} // case D
    type=out0;
    return area=Disk_Inside_Triangle(B,P,radius); // case D
}
template<class TV,class T>
inline TV Random_Sample_In_Triangle(RANDOM_NUMBERS<T>& random,const TV& B,const TV& C)
{
    T a=random.Get_Uniform_Number(0,1);
    T b=random.Get_Uniform_Number(0,1);
    if(a+b>1){a=1-a;b=1-b;}
    return a*B+b*C;
}
// Triangle OBC, with circle (O,radius) removed
// B lies on circle
// C lies outside circle
// interior of BC is outside circle
// Each random sample succeeds with probability >= 2/3
template<class TV,class T>
inline TV Random_Sample_In_Cut_Triangle(RANDOM_NUMBERS<T>& random,const TV& B,const TV& C,T radius)
{
    TV D=C*(radius/C.Magnitude());
    while(1)
    {
        TV X=Random_Sample_In_Triangle(random,B-D,C-D)+D;
        if(X.Magnitude_Squared()>=sqr(radius)) return X;
    }
}
//#####################################################################
// Function Choose_Feasible_Point
//#####################################################################
template<class T> auto VORONOI_DIAGRAM<T>::PIECE_HELPER::
Choose_Feasible_Point(RANDOM_NUMBERS<T>& random,T radius) const -> TV
{
    assert(type!=unset && type!=empty);
    if(type==no_disc) return Random_Sample_In_Triangle(random,B,C)+A; // case A
    if(type==full_disc){
        TV u=C-B;
        TV P=B-clamp(B.Dot(u)/(u).Magnitude_Squared(),-(T)1,(T)0)*u;
        P*=radius/P.Magnitude();
        T p=random.Get_Uniform_Number(0,area);
        T AT=(T).5*(B-P).Cross(C-P).x;
        if(p<=AT) return Random_Sample_In_Triangle(random,B-P,C-P)+P+A;
        p-=AT;
        if(p<=Disk_Inside_Triangle(B,P,radius))
            return Random_Sample_In_Cut_Triangle(random,P,B,radius)+A;
        return Random_Sample_In_Cut_Triangle(random,P,C,radius)+A;}
    if(type==out0) return Random_Sample_In_Cut_Triangle(random,B+aux0*(C-B),B,radius)+A;
    if(type==out1) return Random_Sample_In_Cut_Triangle(random,B+aux0*(C-B),C,radius)+A;

    // case E
    if(random.Get_Uniform_Number(0,1)>aux2)
        return Random_Sample_In_Cut_Triangle(random,B+aux0*(C-B),B,radius)+A;
    return Random_Sample_In_Cut_Triangle(random,B+aux1*(C-B),C,radius)+A;
}
//#####################################################################
// Function Init
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Init(const RANGE<TV>& box,T radius_input)
{
    radius=radius_input;
    bounding_box=box;
    TV e=box.Edge_Lengths();

    // Four points far enough out from the corners that no infinite edge can reach the box.

    TV P[4];
    P[0]=box.Center()+e*(T)2;
    P[2]=box.Center()-e*(T)2;
    P[1]=TV(P[2].x,P[0].y);
    P[3]=TV(P[0].x,P[2].y);
    TV E=random.Get_Uniform_Vector(box);

    CELL *cellP[4]={};
    COEDGE *PE[4]={},*EP[4]={},*PQ[4]={},*QP[4]={};
    VERTEX* V[4]={};
    for(int i=0;i<4;i++){
        PE[i]=new COEDGE;
        EP[i]=new COEDGE;
        PQ[i]=new COEDGE;
        QP[i]=new COEDGE;
        PE[i]->pair=EP[i];
        EP[i]->pair=PE[i];
        PQ[i]->pair=QP[i];
        QP[i]->pair=PQ[i];
        cellP[i]=new CELL;
        cellP[i]->X=P[i];
        cellP[i]->outside=true;
        cellP[i]->coedge=PE[i];
        cellP[i]->type=type_outside;
        cells.Append(cellP[i]);
        V[i]=new VERTEX;
        V[i]->coedge=PE[i];}

    CELL* cellE=new CELL;
    cellE->X=E;
    cellE->outside=false;
    cellE->coedge=EP[0];
    cells.Append(cellE);

    for(int i=0;i<4;i++){
        int k=(i+1)%4;
        PE[i]->prev=PQ[i];
        PQ[i]->next=PE[i];
        QP[i]->prev=PE[k];
        PE[k]->next=QP[i];
        EP[i]->next=EP[k];
        EP[k]->prev=EP[i];
        PE[i]->cell=cellP[i];
        EP[i]->cell=cellE;
        PQ[i]->cell=cellP[i];
        QP[i]->cell=cellP[k];
        PE[i]->tail=V[i];
        EP[i]->head=V[i];
        PE[k]->head=V[i];
        EP[k]->tail=V[i];
        PQ[i]->head=V[i];
        QP[i]->tail=V[i];}

    for(int i=0;i<4;i++)
        V[i]->Compute_Point();
    Visualize_State("before tiles");
    
    for(int i=0;i<4;i++)
        Insert_Coedge(EP[i]);
}
//#####################################################################
// Function Sample_Fully
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Sample_Fully(const RANGE<TV>& box,T radius)
{
    Init(box,radius);
    while(1){
        int p=Choose_Piece();
        if(p<0) break;
        TV X=Choose_Feasible_Point(p);
        Insert_Point(p,X);}
}
//#####################################################################
// Function Get_Samples
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Get_Samples(ARRAY<TV>& X) const
{
    X.Resize(cells.m);
    for(int i=0;i<cells.m;i++)
        X(i)=cells(i)->X;
}
namespace PhysBAM{
template struct VORONOI_DIAGRAM<float>;
template struct VORONOI_DIAGRAM<double>;
}
