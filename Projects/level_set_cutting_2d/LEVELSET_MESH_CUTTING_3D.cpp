#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Log/LOG.h>
#include <Tools/Math_Tools/cyclic_shift.h>
#include <Tools/Matrices/MATRIX.h>
#include "LEVELSET_MESH_CUTTING_3D.h"

using namespace PhysBAM;

static bool Orientations_Match(VECTOR<int,3> a,VECTOR<int,3> b)
{
    while(a(2)!=b(2)) cyclic_shift(b(0),b(1),b(2));
    return a(1)==b(1);
}

static bool Orientations_Match(VECTOR<int,4> a,VECTOR<int,4> b)
{
    int n=0;
    for(;a(3)!=b(3);n++) cyclic_shift(b);
    return (n&1)!=Orientations_Match(a.Remove_Index(3),b.Remove_Index(3));
}

static void Intersection_Point_Face(const VECTOR<int,4>& e,int l,int a,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    for(int i=0;i<4;i++)
        if(e(i)!=l){
            LEVELSET_MESH_CUTTING_3D::TET t={e,e};
            t.indices(i)=a;
            cut_mesh.Append(t);}
}

static void Intersection_Edge_Edge(const VECTOR<int,4>& e,int i,int j,int k,int l,int a,int b,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    typedef VECTOR<int,4> E;
    if(!Orientations_Match(E(i,j,k,l),e)) exchange(i,j);
    LEVELSET_MESH_CUTTING_3D::TET r={e,E(a,j,k,b)},s={e,E(a,k,i,b)},t={e,E(a,i,l,b)},u={e,E(a,l,j,b)};
    cut_mesh.Append(r);
    cut_mesh.Append(s);
    cut_mesh.Append(t);
    cut_mesh.Append(u);
}

static void Intersection_Edge_Face(const VECTOR<int,4>& e,VECTOR<int,3> f,VECTOR<int,2> s,int a,int b,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    typedef VECTOR<int,4> E;
    int l=e.Sum()-f.Sum(),i=s(0)+s(1)-l;
    if(!Orientations_Match(f.Append(l),e)) exchange(f(0),f(1));
    while(f(0)!=i) cyclic_shift(f);
    int j=f(1),k=f(2);
    LEVELSET_MESH_CUTTING_3D::TET r={e,E(l,a,j,b)},t={e,E(a,i,j,b)},u={e,E(k,a,l,b)},v={e,E(k,i,a,b)},w={e,E(l,k,j,b)};
    cut_mesh.Append(r);
    cut_mesh.Append(t);
    cut_mesh.Append(u);
    cut_mesh.Append(v);
    cut_mesh.Append(w);
}

static void Intersection_Face_Face(const VECTOR<int,4>& e,VECTOR<int,3> f,VECTOR<int,3> g,int a,int b,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    typedef VECTOR<int,4> E;
    int l=e.Sum()-f.Sum(),i=e.Sum()-g.Sum();
    if(!Orientations_Match(f.Append(l),e)) exchange(f(0),f(1));
    while(f(0)!=i) cyclic_shift(f);
    int j=f(1),k=f(2);
    LEVELSET_MESH_CUTTING_3D::TET r={e,E(b,k,j,a)},t={e,E(l,b,j,a)},u={e,E(l,k,b,a)},v={e,E(i,l,j,a)},w={e,E(i,k,l,a)};
    cut_mesh.Append(r);
    cut_mesh.Append(t);
    cut_mesh.Append(u);
    cut_mesh.Append(v);
    cut_mesh.Append(w);
}

static void Double_Cuts_1(const LEVELSET_MESH_CUTTING_3D::TET& tet,VECTOR<int,2> s,int a,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    for(int p=0;p<4;p++)
        if(s.Contains(tet.indices(p))){
            LEVELSET_MESH_CUTTING_3D::TET t=tet;
            t.indices(p)=a;
            cut_mesh.Append(t);}
}

static void Double_Cuts_2(const LEVELSET_MESH_CUTTING_3D::TET& tet,VECTOR<int,2> s0,VECTOR<int,2> s1,int a,int b,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    typedef VECTOR<int,4> E;
    int i=s0.Contains(s1.x)?s1.x:s1.y,j=s0.Sum()-i,k=s1.Sum()-i,l=tet.indices.Sum()-i-j-k;
    if(!Orientations_Match(E(i,j,k,l),tet.indices)){
        exchange(j,k);
        exchange(a,b);}
    assert(i!=l);
    if(j<k){
        LEVELSET_MESH_CUTTING_3D::TET r={tet.parent,E(i,a,b,l)},t={tet.parent,E(b,a,j,l)},u={tet.parent,E(k,b,j,l)};
        cut_mesh.Append(r);
        cut_mesh.Append(t);
        cut_mesh.Append(u);}
    else{
        LEVELSET_MESH_CUTTING_3D::TET r={tet.parent,E(i,a,b,l)},t={tet.parent,E(b,a,k,l)},u={tet.parent,E(a,j,k,l)};
        cut_mesh.Append(r);
        cut_mesh.Append(t);
        cut_mesh.Append(u);}
}

static void Double_Cuts_3(const LEVELSET_MESH_CUTTING_3D::TET& tet,VECTOR<int,2> s0,VECTOR<int,2> s1,VECTOR<int,2> s2,int a,int b,int c,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    typedef VECTOR<int,4> E;
    int l=s0.Contains(s1.x)?s1.x:s1.y,i=s0.Sum()-l,j=s1.Sum()-l,k=s2.Sum()-l;
    if(i>j){exchange(i,j);exchange(a,b);}
    if(i>k){exchange(i,k);exchange(a,c);}
    if(!Orientations_Match(E(i,j,k,l),tet.indices)){
        exchange(j,k);
        exchange(b,c);}
    if(j<k){
        LEVELSET_MESH_CUTTING_3D::TET r={tet.parent,E(a,b,c,l)},t={tet.parent,E(a,b,i,c)},u={tet.parent,E(b,j,i,c)},v={tet.parent,E(i,j,k,c)};
        cut_mesh.Append(r);
        cut_mesh.Append(t);
        cut_mesh.Append(u);
        cut_mesh.Append(v);}
    else{
        LEVELSET_MESH_CUTTING_3D::TET r={tet.parent,E(a,b,c,l)},t={tet.parent,E(a,b,i,c)},u={tet.parent,E(b,k,i,c)},v={tet.parent,E(i,j,k,b)};
        cut_mesh.Append(r);
        cut_mesh.Append(t);
        cut_mesh.Append(u);
        cut_mesh.Append(v);}
}

static void Double_Cuts_4(const LEVELSET_MESH_CUTTING_3D::TET& tet,VECTOR<int,2> s0,VECTOR<int,2> s1,VECTOR<int,2> s2,VECTOR<int,2> s3,
    int a,int b,int c,int d,ARRAY<LEVELSET_MESH_CUTTING_3D::TET>& cut_mesh)
{
    typedef VECTOR<int,4> E;
    if(!s0.Contains(s1.x) && !s0.Contains(s1.y)){exchange(s1,s2);exchange(b,c);}
    if(!s0.Contains(s2.x) && !s0.Contains(s2.y)){exchange(s2,s3);exchange(c,d);}
    int i=s0.Contains(s1.x)?s1.x:s1.y,k=s0.Sum()-i,l=s1.Sum()-i,j=s2.Sum()-k;
    assert((s3==VECTOR<int,2>(j,l) || s3==VECTOR<int,2>(l,j)));
    if(k>l){
        exchange(i,j);
        exchange(k,l);
        exchange(a,d);
        exchange(b,c);}
    if(!Orientations_Match(E(i,j,k,l),tet.indices)){
        exchange(i,j);
        exchange(a,c);
        exchange(b,d);}
    if(i<j){
        LEVELSET_MESH_CUTTING_3D::TET r={tet.parent,E(i,j,c,d)},t={tet.parent,E(i,c,a,d)},u={tet.parent,E(i,a,b,d)},v={tet.parent,E(a,c,k,d)},w={tet.parent,E(a,k,b,d)},x={tet.parent,E(k,l,b,d)};
        cut_mesh.Append(r);
        cut_mesh.Append(t);
        cut_mesh.Append(u);
        cut_mesh.Append(v);
        cut_mesh.Append(w);
        cut_mesh.Append(x);}
    else{
        LEVELSET_MESH_CUTTING_3D::TET r={tet.parent,E(i,b,j,a)},t={tet.parent,E(b,d,j,a)},u={tet.parent,E(d,c,j,a)},v={tet.parent,E(a,c,k,d)},w={tet.parent,E(a,k,b,d)},x={tet.parent,E(k,l,b,d)};
        cut_mesh.Append(r);
        cut_mesh.Append(t);
        cut_mesh.Append(u);
        cut_mesh.Append(v);
        cut_mesh.Append(w);
        cut_mesh.Append(x);}
}

void LEVELSET_MESH_CUTTING_3D::Subdivide(const ARRAY<TV_INT4>& mesh,ARRAY<T>& phi0,ARRAY<T>& phi1,ARRAY<TET>& cut_mesh)
{
    enum CUT_TYPE {cut_none,cut_a,cut_b,cut_ab,cut_ba,cut_int};
    int flip[6]={cut_none,cut_a,cut_b,cut_ba,cut_ab,cut_int};

    int num_nodes=phi0.m;
    HASHTABLE<VECTOR<int,2>,PAIR<int,int> > edge_hash;
    HASHTABLE<VECTOR<int,3>,int> tri_hash;
    ARRAY<PAIR<TV_INT,TV> > weights(num_nodes);

    for(int i=0;i<mesh.m;i++){
        TV_INT4 e=mesh(i);
        for(int j=0;j<4;j++)
            for(int k=j+1;k<4;k++){
                VECTOR<int,2> s(e(j),e(k));
                if(edge_hash.Contains(s)) continue;
                VECTOR<T,2> p0(phi0.Subset(s)),p1(phi1.Subset(s));
                bool cross0=(p0(0)>0 && p0(1)<0) || (p0(0)<0 && p0(1)>0);
                bool cross1=(p1(0)>0 && p1(1)<0) || (p1(0)<0 && p1(1)>0);
                PAIR<int,int> z(-1,0);
                if(!p0(0) && !p0(1) && cross1){
                    T interp=p1(0)/(p1(0)-p1(1));
                    z.y=cut_int;
                    z.x=phi0.Append(0);
                    phi1.Append(0);
                    weights.Append(PAIR<TV_INT,TV>(TV_INT(e(j),e(k),-1),TV(1-interp,interp,0)));}
                else if(!p1(0) && !p1(1) && cross0){
                    T interp=p0(0)/(p0(0)-p0(1));
                    z.y=cut_int;
                    z.x=phi0.Append(0);
                    phi1.Append(0);
                    weights.Append(PAIR<TV_INT,TV>(TV_INT(e(j),e(k),-1),TV(1-interp,interp,0)));}
                else if(!cross0) z.y=cross1?cut_b:cut_none;
                else if(!cross1) z.y=cut_a;
                else{
                    T interp=p0(0)/(p0(0)-p0(1));
                    T p1_at_c0=p1(0)+interp*(p1(1)-p1(0));
                    if(p1_at_c0==0) z.y=cut_int;
                    else if((p1_at_c0>0)==(p1(0)>0)) z.y=cut_ab;
                    else z.y=cut_ba;
                    z.x=phi0.Append(0);
                    phi1.Append(p1_at_c0);
                    weights.Append(PAIR<TV_INT,TV>(TV_INT(e(j),e(k),-1),TV(1-interp,interp,0)));}
                edge_hash.Set(s,z);
                edge_hash.Set(s.Reversed(),PAIR<int,int>(z.x,flip[z.y]));}}

    for(int i=0;i<mesh.m;i++){
        TV_INT4 e=mesh(i);
        for(int j=0;j<4;j++){
            TV_INT t=e.Remove_Index(j),st=t.Sorted();
            if(tri_hash.Contains(st)) continue;
//            LOG::cout<<"NEW"<<std::endl;
            bool not_inside=false;
            int ring=0,corner0=0,corner1=0;
            for(int k=0;k<3;k++){
                if(phi0(t(k))==0){
                    corner0++;
                    if(phi1(t(k))==0){not_inside=true;break;}
                    ring=ring*16+1;}
                else if(phi1(t(k))==0){
                    corner1++;
                    ring=ring*16+2;}
                VECTOR<int,2> s(t(k),t((k+1)%3));
                switch(edge_hash.Get(s).y){
                    case cut_none:break;
                    case cut_ba:ring=ring*16+2;
                    case cut_a:ring=ring*16+1;break;
                    case cut_ab:ring=ring*16+1;
                    case cut_b:ring=ring*16+2;break;
                    case cut_int:not_inside=true;break;}
//                LOG::cout<<phi0(t(k))<<" "<<phi1(t(k))<<"   "<<edge_hash.Get(s).y<<"  "<<ring<<"  "<<t<<" -> "<<phi0.m<<std::endl;

}
            if(!not_inside && corner0<2 && corner1<2 && (ring==0x1212 || ring==0x2121)){
                TV p0(phi0.Subset(st)),p1(phi1.Subset(st));
                MATRIX<T,2> M(p0(0)-p0(2),p1(0)-p1(2),p0(1)-p0(2),p1(1)-p1(2));
                VECTOR<T,2> R(-p0(2),-p1(2)),L(M.Inverse_Times(R));
                tri_hash.Set(st,phi0.Append(0));
                phi1.Append(0);
                weights.Append(PAIR<TV_INT,TV>(st,L.Append(1-L.Sum())));}
            else tri_hash.Set(st,-1);}}

    ARRAY<TET> temp_mesh;
    ARRAY<PAIR<int,TV_INT> > hits;
    
    for(int i=0;i<mesh.m;i++){
        TV_INT4 e=mesh(i);
        TV4 p0(phi0.Subset(e)),p1(phi1.Subset(e));
        TET t={e,e};
        hits.Remove_All();

        if(phi0.Min()>=0 || phi0.Max()<=0 || phi1.Min()>=0 || phi1.Max()<=0){
            temp_mesh.Append(t);
            continue;}

        for(int j=0;j<4;j++)
            if(!p0(j) && !p1(j))
                hits.Append(PAIR<int,TV_INT>(e(j),TV_INT(e(j),-1,-1)));

        for(int j=0;j<4;j++)
            for(int k=j+1;k<4;k++){
                VECTOR<int,2> s(e(j),e(k));
                PAIR<int,int> z;
                edge_hash.Get(s,z);
                if(z.y==cut_int)
                    hits.Append(PAIR<int,TV_INT>(z.x,TV_INT(s(0),s(1),-1)));}

        for(int j=0;j<4;j++){
            TV_INT t=e.Remove_Index(j),st=t.Sorted();
            int z=tri_hash.Get(st);
            if(z>=0)
                hits.Append(PAIR<int,TV_INT>(z,st));}

        // Must penetrate exactly twice
        if(hits.m!=2){
            if(hits.m!=0) LOG::cout<<"hits "<<hits<<"   "<<e<<std::endl;

            temp_mesh.Append(t);
            continue;}

        // Penetrates through node
        if(hits(0).y.y==-1){
            if(hits(1).y.z==-1){
                temp_mesh.Append(t);
                continue;}
            Intersection_Point_Face(e,hits(0).y.x,hits(1).x,temp_mesh);
            continue;}

        // Penetrates through edge
        if(hits(0).y.z==-1){
            if(hits(1).y.z==-1){
                Intersection_Edge_Edge(e,hits(0).y.x,hits(0).y.y,hits(1).y.x,hits(1).y.y,hits(0).x,hits(1).x,temp_mesh);
                continue;}

            Intersection_Edge_Face(e,hits(1).y,hits(0).y.Remove_Index(2),hits(1).x,hits(0).x,temp_mesh);
            continue;}

        Intersection_Face_Face(e,hits(0).y,hits(1).y,hits(0).x,hits(1).x,temp_mesh);}

    ARRAY<PAIR<int,VECTOR<int,2> > > dce;
    for(int i=0;i<temp_mesh.m;i++){
        dce.Remove_All();
        for(int j=0;j<4;j++)
            for(int k=j+1;k<4;k++){
                VECTOR<int,2> s(temp_mesh(i).indices(j),temp_mesh(i).indices(k));
                PAIR<int,int> z;
                if(edge_hash.Get(s,z) && z.x>=0)
                    dce.Append(PAIR<int,VECTOR<int,2> >(z.x,s));}

        PHYSBAM_ASSERT(dce.m<=4);
        if(dce.m==0){
            cut_mesh.Append(temp_mesh(i));
            continue;}

        if(dce.m==1){
            Double_Cuts_1(temp_mesh(i),dce(0).y,dce(0).x,cut_mesh);
            continue;}

        if(dce.m==2){
            Double_Cuts_2(temp_mesh(i),dce(0).y,dce(1).y,dce(0).x,dce(1).x,cut_mesh);
            continue;}

        if(dce.m==3){
            Double_Cuts_3(temp_mesh(i),dce(0).y,dce(1).y,dce(2).y,dce(0).x,dce(1).x,dce(2).x,cut_mesh);
            continue;}

        Double_Cuts_4(temp_mesh(i),dce(0).y,dce(1).y,dce(2).y,dce(3).y,dce(0).x,dce(1).x,dce(2).x,dce(3).x,cut_mesh);}

    for(int i=0;i<cut_mesh.m;i++){
        for(int j=0;j<4;j++){
            int k=cut_mesh(i).indices(j);
            int pk=cut_mesh(i).parent.Find(k);
            if(pk>=0){cut_mesh(i).weights(j)(pk)=1;continue;}
            PAIR<TV_INT,TV>& pr=weights(k);
            for(int l=0;l<3;l++)
                if(pr.x(l)>=0)
                    cut_mesh(i).weights(j)(cut_mesh(i).parent.Find(pr.x(l)))=pr.y(l);}}
}

