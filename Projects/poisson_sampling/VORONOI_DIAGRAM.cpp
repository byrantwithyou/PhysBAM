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
template<class V,class F>
inline void For_Each_Outgoing_Coedge_Reverse(V* v,F func)
{
    auto* a=v->first_coedge,*b=a->pair->next,*c=b->pair->next;
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
// Function Print
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::VERTEX::
Print() const
{
    const char * state_names[] = {"unknown","in","out"};
    LOG::printf("vertex name=%c X=%P state=%s first_coedge=%c\n",
        name,X,state_names[(int)state],first_coedge?first_coedge->name:'-');
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::COEDGE::
Print() const
{
    LOG::printf("coedge name=%c pair=%c next=%c prev=%c head=%c tail=%c cell=%c\n",
        name,pair?pair->name:'-',next?next->name:'-',prev?prev->name:'-',head?head->name:'-',tail?tail->name:'-',cell?cell->name:'-');
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::CELL::
Print() const
{
    const char * state_names[] = {"non_incident","incident"};
    LOG::printf("cell name=%c X=%P state=%s first_coedge=%c\n",
        name,X,state_names[(int)state],first_coedge?first_coedge->name:'-');
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::PIECE_HELPER::
Print() const
{
    const char* type_name[]={"unset","empty","no_disc","full_disc","out0","out1","both_out"};

    LOG::printf(" type=%s A=%P B=%P C=%P aux0=%P aux1=%P aux2=%P",
        type_name[type],A,B,C,aux0,aux1,aux2);
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::PIECE::
Print() const
{
    LOG::printf("piece coedge=%c this_area=%P subtree_area=%P",
        coedge?coedge->name:'-',this_area,subtree_area);
    h.Print();
    LOG::printf("\n");
}
//#####################################################################
// Function Print
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::CLIPPED_PIECE::
Print() const
{
    LOG::printf("piece coedge=%c this_area=%P subtree_area=%P [",
        coedge?coedge->name:'-',this_area,subtree_area);
    for(int i=0;i<num_sub_pieces;i++)
        sub_pieces[i].Print();
    LOG::printf(" ]\n");
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
    printf("try %c -> %c\n",ce->name,ceR->name);
    Discover_Inside(in_ce,adj_ce,out_v,in_v,ceR,new_pt);
    printf("try %c -> %c\n",ce->name,ceL->name);
    Discover_Inside(in_ce,adj_ce,out_v,in_v,ceL,new_pt);
    printf("end %c\n",ce->name);
}
//#####################################################################
// Function Insert_Point
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Insert_Point(int p,const TV& new_pt)
{
    if(p>=first_clipped_piece_index)
        Insert_Point(clipped_pieces(p-first_clipped_piece_index).coedge,new_pt);
    else Insert_Point(pieces(p).coedge,new_pt);
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
    For_Each_Outgoing_Coedge_Reverse(v,[&](COEDGE* ce){Discover_Inside(in_coedges,adj_coedges,out_vertices,in_vertices,ce,new_pt);});

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

    for(int i=0;i<adj_coedges.m;i++) adj_coedges(i)->state=1;
    for(int i=0;i<in_coedges.m;i++) in_coedges(i)->state=2;
    Visualize_State("adj_coedges in_coedges");
    for(int i=0;i<adj_coedges.m;i++) adj_coedges(i)->state=0;
    for(int i=0;i<in_coedges.m;i++) in_coedges(i)->state=0;
    
    adj_coedges.Append(+adj_coedges(0)); // + to avoid aliasing
    for(int i=0;i<adj_coedges.m-1;i++){
        COEDGE* ce=adj_coedges(i)->pair;
        if(ce->tail) puts("should not happen ce->tail");
        if(ce->tail) ce->tail->state=unknown;
        ce->cell->first_coedge=ce;

        VERTEX* H=new VERTEX;
        vertices.insert(H);
        ce->head=H;
        ce->pair->tail=H;
        H->state=unknown;
        H->first_coedge=ce->pair;
        H->Compute_Point(new_pt,ce->cell->X,ce->pair->cell->X);
        
        COEDGE* A=new COEDGE;
        coedges.insert(A);
        new_coedges.Append(A);
        ce->next=A;
        A->prev=ce;
        A->next=adj_coedges(i+1);
        A->next->prev=A;
        A->cell=ce->cell;
        A->cell->first_coedge=A;
        A->tail=H;
        A->state=2;
        ce->cell->state=non_incident;

        COEDGE* B=new COEDGE;
        coedges.insert(B);
        new_coedges.Append(B);
        A->pair=B;
        B->pair=A;
        B->head=H;
        B->cell=cell;
        B->state=1;
        Visualize_State("pass");}

    LOG::printf("before connect\n");
    Print();

    for(int i=0;i<in_coedges.m;i++){
        Remove_Coedge(in_coedges(i));
        coedges.erase(in_coedges(i));
        coedges.erase(in_coedges(i)->pair);
        delete in_coedges(i)->pair;
        delete in_coedges(i);}

    for(int i=0;i<adj_coedges.m-1;i++){
        COEDGE* D=adj_coedges(i+1)->pair->next;
        COEDGE* E=adj_coedges(i)->pair->next->pair;
        LOG::printf("vertex %c %P\n",adj_coedges(i)->tail->name,adj_coedges(i)->tail->X);
        LOG::printf("base: %c %c %c   %c %c\n",
            adj_coedges(i)->name,
            adj_coedges(i+1)->name,
            D->name,
            E->name,
            D->pair->name
        );
//        COEDGE* A=adj_coedges(i)->pair->next;
//        COEDGE* B=A->pair;
//        COEDGE* C=A->next;
        D->next=adj_coedges(i);
        D->next->prev=D;
        D->head=D->next->tail;
        E->next=D->pair;
        D->pair->prev=E;
        D->pair->tail=E->head;
//        B->tail=A->head=C->tail;
//        B->prev=C->pair->next->pair;
//        B->prev->next=B;
    }

    cell->first_coedge=adj_coedges(0)->pair->next->pair;

    Visualize_State("connections");
    
    for(int i=0;i<adj_coedges.m-1;i++)
        Update_Coedge(adj_coedges(i));

    for(int i=0;i<new_coedges.m;i++)
        Insert_Coedge(new_coedges(i));

    for(int i=0;i<in_vertices.m;i++){
        vertices.erase(in_vertices(i));
        delete in_vertices(i);}
    Sanity_Checks();
    Visualize_State("done");
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
    VECTOR<T,3> color_r[2]={{1,1,0},{0,1,0}};
    VECTOR<T,3> color_rc(1,.5,0);
    VECTOR<T,3> color_c[3]={{1,1,0},{1,0,0},{0,1,0}};
    VECTOR<T,3> color_cv(0,1,0);
    VECTOR<T,3> color_v[3]={{0,0,1},{1,0,0},{0,1,0}};
    T unbounded_length = .5;
    T offset = .01;

    VECTOR<T,3> color_p[]={{1,.5,.5},{.5,1,.5},{.5,.5,1},{0,1,0},{1,.5,0},{1,1,0},{0,.5,0}};
    for(const auto& it:pieces)
    {
        Add_Debug_Object(VECTOR<TV,3>(TV(),it.h.B,it.h.C)+it.h.A,color_p[it.h.type]);
    }
    for(const auto& it:clipped_pieces)
    {
        for(int i=0;i<it.num_sub_pieces;i++)
            Add_Debug_Object(VECTOR<TV,3>(TV(),it.sub_pieces[i].B,it.sub_pieces[i].C)+it.sub_pieces[i].A,color_p[it.sub_pieces[i].type]);
    }
    
    for(auto it:coedges)
    {
        const VERTEX* h=it->head;
        const VERTEX* t=it->tail;
        assert(h || t);
        TV XH,XT;
        if(h) XH=h->X;
        else{
            TV dir = (it->cell->X - it->pair->cell->X).Rotate_Clockwise_90().Normalized();
            XH=XT=t->X+dir*unbounded_length;}

        if(t) XT=t->X;
        else{
            TV dir = (it->cell->X - it->pair->cell->X).Rotate_Clockwise_90().Normalized();
            XT=XH-dir*unbounded_length;}

        TV dir=(XH-XT).Normalized()*offset,dir_in=dir.Rotate_Counterclockwise_90();

        Add_Debug_Object(VECTOR<TV,2>(XH-dir+dir_in,XT+dir+dir_in),color_c[it->state]);
        Add_Debug_Object(VECTOR<TV,2>((XH+XT)/2+dir_in,it->cell->X),VECTOR<T,3>(1,0,1));
        char name[2]={it->name};
        Add_Debug_Text((XH+XT)/2+dir_in*4,name,color_c[it->state]);
    }

    for(auto it:cells){
        Add_Debug_Particle(it->X,color_r[it->state]);
        char name[2]={it->name};
        Add_Debug_Text(it->X,name,color_r[it->state]);}

    for(auto it:vertices){
        Add_Debug_Particle(it->X,color_v[it->state]);
        char name[2]={it->name};
        Add_Debug_Text(it->X,name,color_v[it->state]);}

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
    LOG::printf("sizes %i %zi %zi %i\n",cells.m,coedges.size(),vertices.size(),pieces.m);
    for(auto c:cells) c->Print();
    for(auto ce:coedges) ce->Print();
    for(auto v:vertices) v->Print();
    for(auto p:pieces) p.Print();
    for(auto p:clipped_pieces) p.Print();
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
    printf("insert coedge %c %c %i\n",ce->name,ce->cell->name,ce->cell->type);
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

    ARRAY<TV> polygon;
    polygon.Append(ce->cell->X);
    polygon.Append(ce->tail->X);
    polygon.Append(ce->head->X);
    LOG::printf("before clip %P\n",polygon);
    for(int i=0;i<polygon.m;i++) Add_Debug_Object(VECTOR<TV,2>(polygon(i),polygon((i+1)%polygon.m)),VECTOR<T,3>(1,1,0));
    INTERSECTION::Box_Polygon_Intersection(bounding_box,polygon,false);
    LOG::printf("after clip %P\n",polygon);
    for(int i=0;i<polygon.m;i++) Add_Debug_Object(VECTOR<TV,2>(polygon(i),polygon((i+1)%polygon.m)),VECTOR<T,3>(0,1,0));
    Flush_Frame<TV>("clip");
    
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

    p.this_area=area;
    VECTOR<T,3> color_p[]={{1,.5,.5},{.5,1,.5},{.5,.5,1},{0,1,0},{1,.5,0},{1,1,0},{0,.5,0}};
    for(int i=0;i<p.num_sub_pieces;i++)
        Add_Debug_Object(VECTOR<TV,3>(TV(),p.sub_pieces[i].B,p.sub_pieces[i].C)+p.sub_pieces[i].A,color_p[p.sub_pieces[i].type]);
    Flush_Frame<TV>("Insert");
    
    int i=clipped_pieces.Append(p);
    Update_Clipped_Piece_Tree(i,area);
    ce->piece=first_clipped_piece_index+i;
}
//#####################################################################
// Function Remove_Coedge
//#####################################################################
template<class T> void VORONOI_DIAGRAM<T>::
Remove_Coedge(COEDGE* ce)
{
    assert((ce->cell->type==type_outside)==(ce->piece<0));
    if(ce->piece>=first_clipped_piece_index)
        Remove_Clipped_Piece(ce->piece-first_clipped_piece_index);
    else if(ce->piece>=0) Remove_Piece(ce->piece);
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
        Update_Clipped_Piece_Tree(p,-rem.this_area);
    else{
        Update_Clipped_Piece_Tree(clipped_pieces.m-1,-last.this_area);
        T diff=last.this_area-rem.this_area;
        last.subtree_area=rem.subtree_area;
        rem=last;
        Update_Clipped_Piece_Tree(p,diff);}
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
//#####################################################################
// Function Choose_Feasible_Point
//#####################################################################
template<class T> auto VORONOI_DIAGRAM<T>::
Choose_Feasible_Point(int p) const -> TV
{
    if(p>=first_clipped_piece_index){
        const CLIPPED_PIECE& cp=clipped_pieces(p-first_clipped_piece_index);
        T r=random.Get_Uniform_Number(0,cp.this_area);
        int s=0;
        for(;s<cp.num_sub_pieces-1;s++){
            r-=cp.sub_pieces[s].this_area;
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
        return this_area=(T).5*B.Cross(C).x;} // Just a triangle: case A

    T B_mag2_min_r2=B.Magnitude_Squared()-sqr(radius);
    T C_mag2_min_r2=C.Magnitude_Squared()-sqr(radius);
    bool B_in_disk=B_mag2_min_r2<=0;
    bool C_in_disk=C_mag2_min_r2<=0;
    if(B_in_disk && C_in_disk){
        type=empty;
        return this_area=0;} // entire triangle in circle: case B

    QUADRATIC<T> quad((B-C).Magnitude_Squared(),-2*B.Dot(B-C),B_mag2_min_r2);
    quad.Compute_Roots();
    if(quad.roots==1){quad.roots=2;quad.root2=quad.root1;}
    if(quad.roots==0 || quad.root1>=1 || quad.root2<=0){
        type=full_disc;
        aux0=(quad.root1+quad.root2)/2; // for point sampling
        aux1=Disk_Inside_Triangle(B,C,radius);
        return this_area=aux1;} // case C

    if(quad.root1<=0 && quad.root2>=1){
        type=empty;
        return this_area=0;} // case B (numerical error)

    if(quad.root1>0 && quad.root2<1){ // case E
        type=both_out;
        TV P=B+quad.root1*(C-B),Q=B+quad.root2*(C-B);
        aux0=quad.root1;
        aux1=quad.root2;
        T A0=Disk_Inside_Triangle(B,P,radius);
        T A1=Disk_Inside_Triangle(Q,C,radius);
        aux2=A0/(A0+A1);
        return this_area=A0+A1;}

    aux0=quad.root1>0?quad.root1:quad.root2;
    TV P=B+aux0*(C-B);
    if(B_in_disk){
        type=out1;
        return this_area=Disk_Inside_Triangle(P,C,radius);} // case D
    type=out0;
    return this_area=Disk_Inside_Triangle(B,P,radius); // case D
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
        TV P=B+aux0*(C-B);
        P*=radius/P.Magnitude();
        T p=random.Get_Uniform_Number(0,aux1);
        T AT=(T).5*(B-P).Cross(C-P).x;
        if(p<=AT) return Random_Sample_In_Triangle(random,B-P,C-P)+P+A;
        p-=AT;
        if(p<=Disk_Inside_Triangle(P,B,radius))
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
        coedges.insert(PE[i]);
        coedges.insert(EP[i]);
        coedges.insert(PQ[i]);
        coedges.insert(QP[i]);
        cellP[i]=new CELL;
        cellP[i]->X=P[i];
        cellP[i]->outside=true;
        cellP[i]->first_coedge=PE[i];
        cellP[i]->type=type_outside;
        cells.Append(cellP[i]);
        V[i]=new VERTEX;
        V[i]->first_coedge=PE[i];
        vertices.insert(V[i]);}

    CELL* cellE=new CELL;
    cellE->X=E;
    cellE->outside=false;
    cellE->first_coedge=EP[0];
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

    Print();
    Sanity_Checks();
}
namespace PhysBAM{
template struct VORONOI_DIAGRAM<float>;
template struct VORONOI_DIAGRAM<double>;
}
