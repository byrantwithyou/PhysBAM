//#####################################################################
// Copyright 2016, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/LOG.h>
#include <Tools/Polynomials/QUADRATIC.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "VORONOI_DIAGRAM.h"
using namespace PhysBAM;
template<> char VORONOI_DIAGRAM<float>::next_cell = 'A';
template<> char VORONOI_DIAGRAM<float>::next_vertex = 'G';
template<> char VORONOI_DIAGRAM<float>::next_coedge = 'a';
template<> char VORONOI_DIAGRAM<double>::next_cell = 'A';
template<> char VORONOI_DIAGRAM<double>::next_vertex = 'G';
template<> char VORONOI_DIAGRAM<double>::next_coedge = 'a';
template<class V,class F>
inline void For_Each_Outgoing_Coedge(V* v,F func)
{
    auto* a=v->first_coedge,*b=a->pair->next,*c=b->pair->next;
    func(a);
    func(b);
    func(c);
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
    COEDGE* ce=first_coedge;
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
    COEDGE* ce=first_coedge;
    TV A=ce->cell->X;
    ce=ce->pair->next;
    TV B=ce->cell->X;
    ce=ce->pair->next;
    TV C=ce->cell->X;
    Compute_Point(A,B,C);
}
//#####################################################################
// Function Criterion
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::VERTEX::
Print() const
{
    const char * state_names[] = {"unknown","in","out"};
    LOG::printf("vertex name=%c X=%P state=%s first_coedge=%c\n",
        name,X,state_names[(int)state],first_coedge?first_coedge->name:'-');
}
//#####################################################################
// Function Criterion
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::COEDGE::
Print() const
{
    LOG::printf("coedge name=%c pair=%c next=%c prev=%c head=%c tail=%c cell=%c\n",
        name,pair?pair->name:'-',next?next->name:'-',prev?prev->name:'-',head?head->name:'-',tail?tail->name:'-',cell?cell->name:'-');
}
//#####################################################################
// Function Criterion
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::CELL::
Print() const
{
    const char * state_names[] = {"non_incident","incident"};
    LOG::printf("cell name=%c X=%P state=%s first_coedge=%c\n",
        name,X,state_names[(int)state],first_coedge?first_coedge->name:'-');
}
//#####################################################################
// Function Select_First_In_Vertex
//#####################################################################
template<class T> auto VORONOI_DIAGRAM<T>::
Select_First_In_Vertex(COEDGE* start,const TV& new_pt) const -> VERTEX*
{
    if(!start->head) return start->tail;
    if(!start->tail) return start->head;
    T crit0=start->head->Criterion(new_pt);
    T crit1=start->tail->Criterion(new_pt);
    if(crit0<crit1) return start->head;
    return start->tail;
}
//#####################################################################
// Function Discover_Inside
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Discover_Inside(ARRAY<COEDGE*>& in_ce,ARRAY<COEDGE*>& adj_ce,
    ARRAY<VERTEX*>& out_v,ARRAY<VERTEX*>& in_v,COEDGE* ce,const TV& new_pt)
{
    printf("discover %c\n",ce->name);
    if(!ce->head || ce->head->state!=unknown){
        printf("  early fail\n");
        adj_ce.Append(ce);
        return;}

    VERTEX* v = ce->head;
    COEDGE* ceL = ce->next;
    COEDGE* ceR = ceL->pair->next;
    printf("  v=%c ceL=%c ceR=%c\n",v->name,ceL->name,ceR->name);

    printf("  tests: %i %i %i %i\n",ceL->head && ceL->head->state==in,ceR->head && ceR->head->state==in,ceR->cell->state==incident,v->Criterion(new_pt)>=0);
    if((ceL->head && ceL->head->state==in) // (C4)
        || (ceR->head && ceR->head->state==in) // (C4)
        || ceR->cell->state==incident // (C5)
        || v->Criterion(new_pt)>=0){ // geometric criterion
        printf("  failed - out (%c %c)\n",v->name,ce->name);
        v->state=out;
        out_v.Append(v);
        adj_ce.Append(ce);
        return;}

    printf("  success - in (%c %c %c)\n",v->name,ce->name,ceR->name);
    v->state=in;
    ceR->cell->state=incident;
    in_v.Append(v);
    in_ce.Append(ce);
    if(ceR->head) printf("try %c -> %c\n",ce->name,ceR->name);
    if(ceR->head) Discover_Inside(in_ce,adj_ce,out_v,in_v,ceR,new_pt);
    if(ceL->head) printf("try %c -> %c\n",ce->name,ceL->name);
    if(ceL->head) Discover_Inside(in_ce,adj_ce,out_v,in_v,ceL,new_pt);
    printf("end %c\n",ce->name);
}
//#####################################################################
// Function Insert_Point
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Point(COEDGE* start,const TV& new_pt)
{
    puts("before Insert_Point");
    Print();
//    Add_Debug_Particle(new_pt,VECTOR<T,3>(1,0,1));
    ARRAY<VERTEX*> in_vertices,out_vertices; // out_vertices not needed anymore.
    ARRAY<COEDGE*> in_coedges,adj_coedges,new_coedges;
    
    // Make sure at least one vertex is in in_vertices
    VERTEX* v=Select_First_In_Vertex(start,new_pt);
    in_vertices.Append(v);
    v->state=in;
    For_Each_Outgoing_Coedge(v,[](COEDGE* ce){ce->cell->state=incident;});
    For_Each_Outgoing_Coedge(v,[&](COEDGE* ce){Discover_Inside(in_coedges,adj_coedges,out_vertices,in_vertices,ce,new_pt);});

    puts("after discovery");
    Print();

    // Allocate new cell
    CELL* cell = new CELL;
    cells.Append(cell);
    cell->X = new_pt;

    printf("in_vertices:");
    for(int i=0;i<in_vertices.m;i++) printf(" %c",in_vertices(i)->name);
    printf("\n");

    printf("out_vertices:");
    for(int i=0;i<out_vertices.m;i++) printf(" %c",out_vertices(i)->name);
    printf("\n");

    printf("in_coedges:");
    for(int i=0;i<in_coedges.m;i++) printf(" %c",in_coedges(i)->name);
    printf("\n");

    printf("adj_coedges:");
    for(int i=0;i<adj_coedges.m;i++) printf(" %c",adj_coedges(i)->name);
    printf("\n");

    adj_coedges.Append(+adj_coedges(0)); // + to avoid aliasing
    for(int i=0;i<adj_coedges.m-1;i++){
        COEDGE* ce=adj_coedges(i)->pair;
        if(ce->tail) puts("should not happen ce->tail");
        if(ce->tail) ce->tail->state=unknown;

        VERTEX* H=new VERTEX;
        vertices.insert(H);
        ce->head=H;
        ce->pair->tail=H;
        H->state=unknown;
        H->first_coedge=ce->pair;

        COEDGE* A=new COEDGE;
        coedges.insert(A);
        new_coedges.Append(A);
        ce->next=A;
        A->prev=ce;
        A->next=adj_coedges(i+1);
        A->next->prev=A;
        A->cell=ce->cell;
        A->tail=H;
        ce->cell->state=non_incident;

        COEDGE* B=new COEDGE;
        coedges.insert(B);
        A->pair=B;
        B->pair=A;
        B->head=H;
        B->cell=cell;}

    for(int i=0;i<adj_coedges.m-1;i++){
        COEDGE* A=adj_coedges(i)->prev;
        COEDGE* B=A->pair;
        B->tail=A->head=adj_coedges(i+1)->tail;
        B->prev=A->next->pair->next->pair;
        B->prev->next=B;}

    for(int i=0;i<adj_coedges.m-1;i++){
        COEDGE* A=adj_coedges(i)->prev;
        printf("compute on %c\n",A->head->name);
        A->head->Compute_Point();}

    for(int i=0;i<adj_coedges.m-1;i++)
        Update_Coedge(adj_coedges(i));

    for(int i=0;i<new_coedges.m;i++)
        Insert_Coedge(new_coedges(i));

    for(int i=0;i<in_coedges.m;i++){
        Remove_Coedge(in_coedges(i));
        coedges.erase(in_coedges(i));
        coedges.erase(in_coedges(i)->pair);
        delete in_coedges(i)->pair;
        delete in_coedges(i);}

    for(int i=0;i<in_vertices.m;i++){
        vertices.erase(in_vertices(i));
        delete in_vertices(i);}
    puts("after Insert_Point");
    Print();
}
//#####################################################################
// Function First_Three_Points
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
First_Three_Points(TV A,TV B,TV C)
{
    puts("before First_Three_Points");
    Print();
    if((A-C).Cross(B-C).x<0) exchange(B,C);
    CELL* cellA = new CELL;
    CELL* cellB = new CELL;
    CELL* cellC = new CELL;
    cellA->X=A;
    cellB->X=B;
    cellC->X=C;
    cellA->outside=true;
    cellB->outside=true;
    cellC->outside=true;
    
    COEDGE* AB = new COEDGE;
    COEDGE* BA = new COEDGE;
    COEDGE* AC = new COEDGE;
    COEDGE* CA = new COEDGE;
    COEDGE* BC = new COEDGE;
    COEDGE* CB = new COEDGE;
    AB->cell=AC->cell=cellA;
    BA->cell=BC->cell=cellB;
    CA->cell=CB->cell=cellC;
    cellA->first_coedge=AB;
    cellB->first_coedge=BC;
    cellC->first_coedge=CA;

    AB->pair=BA;
    BA->pair=AB;
    AC->pair=CA;
    CA->pair=AC;
    BC->pair=CB;
    CB->pair=BC;

    AB->next=AC;
    AC->prev=AB;
    CA->next=CB;
    CB->prev=CA;
    BC->next=BA;
    BA->prev=BC;

    VERTEX* V = new VERTEX;
    V->first_coedge=AC;
    For_Each_Outgoing_Coedge(V,[=](COEDGE* ce){ce->tail=V;ce->pair->head=V;});
    V->Compute_Point();

    cells.Append(cellA);
    cells.Append(cellB);
    cells.Append(cellC);
    coedges.insert(AB);
    coedges.insert(BA);
    coedges.insert(AC);
    coedges.insert(CA);
    coedges.insert(BC);
    coedges.insert(CB);
    vertices.insert(V);
    puts("after First_Three_Points");
    Print();
}
//#####################################################################
// Function Visualize_State
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Visualize_State(const char* title) const
{
    VECTOR<T,3> color_r(1,0,0);
    VECTOR<T,3> color_rc(1,.5,0);
    VECTOR<T,3> color_c(1,1,0);
    VECTOR<T,3> color_cv(0,1,0);
    VECTOR<T,3> color_v(0,0,1);
    T unbounded_length = .5;
    T offset = .01;
        
    for(int i=0;i<cells.m;i++)
        Add_Debug_Particle(cells(i)->X,color_r);
        
    for(auto it:vertices)
        Add_Debug_Particle(it->X,color_v);

    for(auto it:coedges)
    {
        const VERTEX* h=it->head;
        const VERTEX* t=it->tail;
        assert(h || t);
        TV XH,XT;
        if(h) XH=h->X;
        else
        {
            TV dir = (it->cell->X - it->pair->cell->X).Rotate_Clockwise_90().Normalized();
            XH=XT=t->X+dir*unbounded_length;
        }

        if(t) XT=t->X;
        else
        {
            TV dir = (it->cell->X - it->pair->cell->X).Rotate_Clockwise_90().Normalized();
            XT=XH-dir*unbounded_length;
        }

        TV dir=(XH-XT).Normalized()*offset,dir_in=dir.Rotate_Counterclockwise_90();
            
        Add_Debug_Object(VECTOR<TV,2>(XH-dir+dir_in,XT+dir+dir_in),color_c);
    }

    Flush_Frame<TV>(title);
}
//#####################################################################
// Function Sanity_Checks
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Sanity_Checks() const
{
    for(int i=0;i<cells.m;i++){
        COEDGE* ce=cells(i)->first_coedge;
        do{
            assert(ce->cell==cells(i));
            ce=ce->next;
        }while(ce && ce!=cells(i)->first_coedge);}

    for(auto ce:coedges){
        assert(ce->pair->pair==ce);
        if(ce->next) assert(ce->next->prev==ce);
        assert(ce->head==ce->pair->tail);
        if(ce->head && ce->tail) assert((ce->tail->X-ce->cell->X).Cross(ce->head->X-ce->cell->X).x>0);
        assert(!ce->next == !ce->head);
        assert(!ce->prev == !ce->tail);
        if(ce->next) assert(ce->next->pair->next->pair->next->pair==ce);}

    for(auto v:vertices){
        if(v) assert(v->first_coedge->tail==v);}
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Print() const
{
    for(int i=0;i<cells.m;i++) cells(i)->Print();
    for(auto ce:coedges) ce->Print();
    for(auto v:vertices) v->Print();
}
//#####################################################################
// Function Update_Piece_Tree
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Update_Piece_Tree(int i,T diff_area)
{
    pieces(i).subtree_area+=diff_area;
    while(i>0){
        i=(i-1)/2;
        pieces(i).subtree_area+=diff_area;}
}
//#####################################################################
// Function Update_Clipped_Piece_Tree
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Update_Clipped_Piece_Tree(int i,T diff_area)
{
    clipped_pieces(i).subtree_area+=diff_area;
    while(i>0){
        i=(i-1)/2;
        clipped_pieces(i).subtree_area+=diff_area;}
}
//#####################################################################
// Function Insert_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Coedge(COEDGE* ce)
{
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
        T area=p.h.Compute(ce,radius,false);
        int i=pieces.Append(p);
        p.this_area=area;
        Update_Piece_Tree(i,area);
        ce->piece=i;}
}
//#####################################################################
// Function Insert_Clipped_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Clipped_Coedge(COEDGE* ce)
{
    assert(bounding_box.Lazy_Inside(ce->cell->X));

    ARRAY<TV> vert;
    vert.Append(ce->cell->X);
    vert.Append(ce->tail->X);
    vert.Append(ce->head->X);
    vert.Append(ce->cell->X);

    TV bound[2]={bounding_box.min_corner,bounding_box.max_corner};
    for(int s=0;s<2;s++){
        int sign=2*s-1;
        for(int a=0;a<2;a++){
            int i=1,j;
            for(;i<vert.m-1;i++)
                if(vert(i)(a)*sign<bound[s](a)*sign)
                    break;
            if(i==vert.m-1) continue;
            for(j=i+1;j<vert.m-1;j++)
                if(vert(j)(a)*sign>=bound[s](a)*sign)
                    break;

            T u=(bound[s](a)-vert(i)(a))/(vert(i-1)(a)-vert(i)(a));
            T v=(bound[s](a)-vert(j)(a))/(vert(j-1)(a)-vert(j)(a));

            TV A=vert(i)+u*(vert(i-1)-vert(i));
            TV B=vert(j)+v*(vert(j-1)-vert(j));

            ARRAY<TV> tmp(vert.Array_View(0,i));
            tmp.Append(A);
            tmp.Append(B);
            tmp.Append_Elements(vert.Array_View(j,vert.m-j));
            tmp.Exchange(vert);}}
    vert.Pop();

    assert(vert.m<=6);
    CLIPPED_PIECE p;
    p.coedge=ce;
    p.num_sub_pieces=vert.m-2;
    T area=0;
    for(int i=0;i<vert.m-2;i++){
        p.sub_pieces[i].A=vert(0);
        p.sub_pieces[i].B=vert(i+1);
        p.sub_pieces[i].C=vert(i+2);
        area+=p.sub_pieces[i].Compute(ce,radius,true);}

    int i=clipped_pieces.Append(p);
    p.this_area=area;
    Update_Clipped_Piece_Tree(i,area);
    ce->piece=first_clipped_piece_index+i;
}
//#####################################################################
// Function Remove_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Remove_Coedge(COEDGE* ce)
{
    if(ce->piece<first_clipped_piece_index) Remove_Piece(ce->piece);
    else Remove_Clipped_Piece(ce->piece-first_clipped_piece_index);
}
//#####################################################################
// Function Remove_Piece
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Remove_Piece(int p)
{
    PIECE& last=pieces.Last(), &rem=pieces(p);
    if(p==pieces.m-1)
        Update_Piece_Tree(p,-rem.this_area);
    else{
        Update_Piece_Tree(pieces.m-1,-last.this_area);
        T diff=last.this_area-rem.this_area;
        last.subtree_area=rem.subtree_area;
        rem=last;
        Update_Piece_Tree(p,diff);}
    pieces.Pop();
}
//#####################################################################
// Function Remove_Piece
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Remove_Clipped_Piece(int p)
{
    CLIPPED_PIECE& last=clipped_pieces.Last(), &rem=clipped_pieces(p);
    if(p==clipped_pieces.m-1)
        Update_Piece_Tree(p,-rem.this_area);
    else{
        Update_Piece_Tree(clipped_pieces.m-1,-last.this_area);
        T diff=last.this_area-rem.this_area;
        last.subtree_area=rem.subtree_area;
        rem=last;
        Update_Piece_Tree(p,diff);}
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
Choose_Piece()
{
    T total_area=0;
    if(pieces.m) total_area+=pieces(0).subtree_area;
    if(clipped_pieces.m) total_area+=clipped_pieces(0).subtree_area;
    T r=random.Get_Uniform_Number(0,total_area);
    int i=0;
    if(r<=clipped_pieces(0).subtree_area){
        r-=clipped_pieces(0).subtree_area;
        while(1){
            if(r<=clipped_pieces(i).this_area) return first_clipped_piece_index+i;
            int j=2*i+1,k=j+1;
            r-=clipped_pieces(i).this_area;
            if(j>=clipped_pieces.m){assert(r<1e-5);return first_clipped_piece_index+i;}
            if(r<=clipped_pieces(j).subtree_area){i=j;continue;}
            r-=clipped_pieces(j).subtree_area;
            if(k>=clipped_pieces.m){assert(r<1e-5);return first_clipped_piece_index+j;}
            i=k;}}
    else{
        while(1){
            if(r<=pieces(i).this_area) return i;
            int j=2*i+1,k=j+1;
            r-=pieces(i).this_area;
            if(j>=pieces.m){assert(r<1e-5);return i;}
            if(r<=pieces(j).subtree_area){i=j;continue;}
            r-=pieces(j).subtree_area;
            if(k>=pieces.m){assert(r<1e-5);return j;}
            i=k;}}
}
template<class TV,class T>
inline T Disk_Inside_Triangle(TV B,TV C,T radius)
{
    return B.Cross(C).x - TV::Angle_Between(B,C)*sqr(radius)/2;
}
//#####################################################################
// Function Compute_Available_Area
//#####################################################################
template<class T> T VORONOI_DIAGRAM<T>::PIECE_HELPER::
Compute(COEDGE* ce,T radius,bool clipped)
{
    if(!clipped){
        A=ce->cell->X;
        B=ce->tail->X-A;
        C=ce->head->X-A;}

    if(ce->cell->outside){
        type=no_disc;
        return B.Cross(C).x;} // Just a triangle: case A

    T B_mag2_min_r2=B.Magnitude_Squared()-sqr(radius);
    T C_mag2_min_r2=C.Magnitude_Squared()-sqr(radius);
    bool B_in_disk=B_mag2_min_r2<=0;
    bool C_in_disk=C_mag2_min_r2<=0;
    if(B_in_disk && C_in_disk){
        type=empty;
        return 0;} // entire triangle in circle: case B

    QUADRATIC<T> quad((B-C).Magnitude_Squared(),-2*B.Dot(B-C),B_mag2_min_r2);
    quad.Compute_Roots();
    if(quad.roots==1){quad.roots=2;quad.root2=quad.root1;}
    if(quad.roots==0 || quad.root1>=1 || quad.root2<=0){
        type=full_disc;
        aux0=(quad.root1+quad.root2)/2; // for point sampling
        aux1=Disk_Inside_Triangle(B,C,radius);
        return aux1;} // case C

    if(quad.root1<=0 && quad.root2>=1){
        type=empty;
        return 0;} // case B (numerical error)

    if(quad.root1>0 && quad.root2<1){ // case E
        type=both_out;
        TV P=B+quad.root1*(C-B),Q=B+quad.root2*(C-B);
        aux0=quad.root1;
        aux1=quad.root2;
        T A0=Disk_Inside_Triangle(B,P,radius);
        T A1=Disk_Inside_Triangle(Q,C,radius);
        aux2=A0/(A0+A1);
        return A0+A1;}

    aux0=quad.root1>0?quad.root1:quad.root2;
    TV P=B+aux0*(C-B);
    if(B_in_disk){
        type=out1;
        return Disk_Inside_Triangle(P,C,radius);} // case D
    type=out0;
    return Disk_Inside_Triangle(B,P,radius); // case D
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
// Function Compute_Available_Area
//#####################################################################
template<class T> auto VORONOI_DIAGRAM<T>::PIECE_HELPER::
Choose_Feasible_Point(RANDOM_NUMBERS<T>& random,T radius) const -> TV
{
    assert(type!=unset && type!=empty);
    if(type==no_disc) return Random_Sample_In_Triangle(random,B,C)+A; // case A
    if(type==full_disc){
        TV P=B+aux0*(C-B);
        P*=radius/P.Magnitude();
        T p=random.Get_Uniform_Number(0,aux1);
        T AT=(B-P).Cross(C-P).x;
        if(p<AT) return Random_Sample_In_Triangle(random,B-P,C-P)+P+A;
        p-=AT;
        if(p<Disk_Inside_Triangle(P,C,radius))
            return Random_Sample_In_Cut_Triangle(random,P,B,radius);
        return Random_Sample_In_Cut_Triangle(random,P,C,radius);}
    if(type==out0) return Random_Sample_In_Cut_Triangle(random,B+aux0*(C-B),C,radius);
    if(type==out1) return Random_Sample_In_Cut_Triangle(random,B+aux0*(C-B),B,radius);

    // case E
    if(random.Get_Uniform_Number(0,1)>aux2)
        return Random_Sample_In_Cut_Triangle(random,B+aux0*(C-B),B,radius);
    return Random_Sample_In_Cut_Triangle(random,B+aux1*(C-B),C,radius);
}
//#####################################################################
// Function Init
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Init(const RANGE<TV>& box)
{
    bounding_box=box;
    TV e=box.Edge_Lengths();

    // Four points far enough out from the corners that no infinite edge can reach the box.
    TV P[4];
    P[0]=box.Center()+e*2;
    P[2]=box.Center()-e*2;
    P[1]=TV(P[0].x,P[2].y);
    P[3]=TV(P[2].x,P[0].y);

    TV Q=random.Get_Uniform_Vector(box);

    // TODO: Construct initial diagram.


    
}
namespace PhysBAM{
template struct VORONOI_DIAGRAM<float>;
template struct VORONOI_DIAGRAM<double>;
}
