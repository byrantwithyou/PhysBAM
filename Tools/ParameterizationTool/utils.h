#ifndef DK_UTILS_H
#define DK_UTILS_H
#include "../../Public_Library/Geometry/Triangulated_surface.h"
#include "Cyclic_list_array.h"
#include <math.h>
#include <Particles/SOLIDS_PARTICLES.h>

using namespace PhysBAM;
#define Pi 3.1415926
#define EPSILON 0.00001
#define MAX_VALUE 99999999
#define FLAG_REMOVED 0
#define FLAG_ACTIVE  1
#define FLAG_ALIVE     1
#define FLAG_DEAD    2
#define VECTOR VECTOR_3D <double>
template <class T0, class T1> class DK_VERTEX_LIST;

bool Compute_Barycentric(const VECTOR p1, const VECTOR p2, const VECTOR p3, const VECTOR p, double &u1, double &u2);
int SameTriples(int ti, int tj, int tk, int tti, int ttj, int ttk, int &notused);

template <class T> class DK_VERTEX
{
public:
    int index, flag;
    bool tag;
    T area, curvature, weight;
    T u, v, real_u, real_v;
    VECTOR_3D <T> temp_P;
    DK_VERTEX() {flag=FLAG_ACTIVE; tag=false; real_u=-1; real_v=-1;}
};
template <class T>
class DK_TRIANGLE
{
public:
    int ti, tj, tk;
    double ctan_i, ctan_j, ctan_k;
    VECTOR normal;
    double area;
    int flag;
    ARRAY <int> inside_vertices;
    void Arrange(int i, int &j, int &k, double &cot_i, double &cot_j, double &cot_k)
    {
        if( i==ti ) {        j=tj; k=tk;    cot_i=ctan_i; cot_j=ctan_j; cot_k=ctan_k; }
        else if( i==tj ) {  j=ti; k=tk; cot_i=ctan_j; cot_j=ctan_i; cot_k=ctan_k; }
        else if( i==tk ) {    j=ti; k=tj; cot_i=ctan_k; cot_j=ctan_i; cot_k=ctan_j; }
        else printf("ERROR: failed to arrange.\n");
    }
};

template <class T0, class T1>
class DK_TRIANGLE_LIST
{
    ARRAY< DK_TRIANGLE<T1> *> old_triangles_info;
    ARRAY< DK_TRIANGLE<T1> *> new_triangles_info;
    SOLIDS_PARTICLES<T0,VECTOR_3D<T0> > *particles;
    TRIANGLE_MESH *mesh;
    DK_VERTEX_LIST<T0, T1> *vlist;
public:
    DK_TRIANGLE_LIST(    SOLIDS_PARTICLES<T0,VECTOR_3D<T0> > *emb_boundary_particles, 
                    TRIANGLE_MESH *currentmesh, 
                    DK_VERTEX_LIST<T0, T1> *vl, 
                    ARRAY < ARRAY <int > > &vertices_info)
    {
        vlist=vl;
        mesh=currentmesh;
        particles=emb_boundary_particles;
        for(int i=1; i<=mesh->triangles.m; i++)
        {
            DK_TRIANGLE<T1> *t=new DK_TRIANGLE<T1>;
            t->ti=mesh->triangles(1, i);
            t->tj=mesh->triangles(2, i);
            t->tk=mesh->triangles(3, i);
            t->flag=FLAG_ACTIVE;

            VECTOR e1, e2, e3;
            e1=particles->X(t->tj)-particles->X(t->ti);
            e2=particles->X(t->tk)-particles->X(t->ti);
            e3=particles->X(t->tk)-particles->X(t->tj);
            t->normal=VECTOR::Cross_Product(e1, e2);
            if( t->normal.Magnitude_Squared()==0)    t->area=0;
            else                                    t->area=t->normal.Normalize();
            e1.Normalize();
            e2.Normalize();
            e3.Normalize();
            double dot;
            dot=VECTOR::Dot_Product(e1, e2);
            t->ctan_i=dot / sqrt( 1-dot*dot);
            dot=VECTOR::Dot_Product(-e1, e3);
            t->ctan_j=dot / sqrt( 1-dot*dot);
            dot=VECTOR::Dot_Product(-e2, -e3);
            t->ctan_k=dot / sqrt( 1-dot*dot);
            old_triangles_info.Append(t);
        }
        if( mesh->triangles.m == vertices_info.m)
            for( int i=1; i<=vertices_info.m; i++)
                old_triangles_info(i)->inside_vertices=vertices_info(i);
    }
    
    ~DK_TRIANGLE_LIST() 
    {for(int i=1; i<=old_triangles_info.m; i++) delete old_triangles_info(i);for( i=1; i<=new_triangles_info.m; i++) delete new_triangles_info(i);}
    
    ARRAY <int> *remove(int i, int &ri, int &rj, int &rk)
    {
        assert(i>0 && i<=old_triangles_info.m); 
        old_triangles_info(i)->flag=FLAG_REMOVED;
        ri=old_triangles_info(i)->ti;
        rj=old_triangles_info(i)->tj;
        rk=old_triangles_info(i)->tk;
        return &(old_triangles_info(i)->inside_vertices);
    }

    bool Add(int tti, int ttj, int ttk, ARRAY <DK_VERTEX<double> *> *free_vlist,
        CYCLIC_ARRAY <DK_VERTEX <double> *> &maps, int p0)
    {
        int ti=maps(tti)->index;
        int tj=maps(ttj)->index;
        int tk=maps(ttk)->index;
        int notused;

        if( ExistTriangle(tti, ttj, ttk, maps)==true)
        {
            printf("Warning: The triangle %d %d %d exists.\n", ti, tj, tk);
            return false;
        }
        DK_TRIANGLE<T1> *t=new     DK_TRIANGLE<T1>;
        t->ti=ti; t->tj=tj; t->tk=tk;
    
        for(int i=1; i<=free_vlist->m; i++)
        {
            if( (*free_vlist)(i)->tag==true)
            {
                double u, v;
                if ( Compute_Barycentric( maps(tti)->temp_P, maps(ttj)->temp_P,
                    maps(ttk)->temp_P, (*free_vlist)(i)->temp_P, u, v))
                {
                    (*free_vlist)(i)->u=u;
                    (*free_vlist)(i)->v=v;
                    t->inside_vertices.Append( (*free_vlist)(i)->index);
                    (*free_vlist)(i)->tag=false;                    
                }
            }
        }
        new_triangles_info.Append(t);
        //next, we need to update the neighbornode list. correct.
        mesh->triangles.Append(ti, tj, tk);
        int tindex=mesh->triangles.m;
        (*mesh->neighbor_nodes)(ti).Append_Unique(tj);
        (*mesh->neighbor_nodes)(ti).Append_Unique(tk);
        (*mesh->neighbor_nodes)(tj).Append_Unique(ti);
        (*mesh->neighbor_nodes)(tj).Append_Unique(tk);
        (*mesh->neighbor_nodes)(tk).Append_Unique(ti);
        (*mesh->neighbor_nodes)(tk).Append_Unique(tj);
        (*mesh->incident_triangles)(ti).Append(tindex);
        (*mesh->incident_triangles)(tj).Append(tindex);
        (*mesh->incident_triangles)(tk).Append(tindex);
        return true;
    }

    void Getlist( ARRAYS<int> &ret, ARRAY< ARRAY <int> > &vertices_info)
    {
        for(int i=1; i<=old_triangles_info.m; i++)
        {
            if( old_triangles_info(i)->flag!=FLAG_REMOVED )
            {
                ret.Append( old_triangles_info(i)->ti, 
                    old_triangles_info(i)->tj,     old_triangles_info(i)->tk);
                vertices_info.Append( old_triangles_info(i)->inside_vertices );
            }
        }
        //it may not work here, check it. the bottleneck
        ARRAYS<int> secondret(3,0);
        ARRAY< ARRAY <int> > second_vinfo;
        secondret.Clean_Up_Memory();
        second_vinfo.Clean_Up_Memory();
        for( i=1; i<=new_triangles_info.m; i++)
        {
            secondret.Append( new_triangles_info(i)->ti, 
                new_triangles_info(i)->tj,     new_triangles_info(i)->tk);
            second_vinfo.Append( new_triangles_info(i)->inside_vertices );
        }
        ret.Append_Elements(secondret);
        vertices_info.Append_Elements(second_vinfo);
    }

    DK_TRIANGLE<T1> * GetInfo(int t){    return old_triangles_info(t);    }

    bool InvalidTriangle(int ti, int tj, int tk, CYCLIC_ARRAY <DK_VERTEX <double> *> &maps, int ddd=0)
    {
        ti=maps(ti)->index;
        tj=maps(tj)->index;
        tk=maps(tk)->index;
        
        //test if the triangle ti, tj, tk has already existed in the current mesh.
        int nNeighbor=0;
        int Neighbor[10][2];
        for( int i=1; i<=(*mesh->incident_triangles)(ti).m; i++)
        {
            int triangle=(*mesh->incident_triangles)(ti)(i);
            int tti=mesh->triangles(1, triangle);
            int ttj=mesh->triangles(2, triangle);
            int ttk=mesh->triangles(3, triangle);
            int notused;
            int degree=SameTriples(ti, tj, tk, tti, ttj, ttk, notused);
            if(degree==3)    return true;
            else if(degree==2)
            {
                for(int j=0; j<nNeighbor; j++)
                    if( Neighbor[j][0]==notused) break;
                if (j==nNeighbor)
                {
                    Neighbor[j][0]=notused;
                    Neighbor[j][1]=tti+ttj+ttk-ti-notused;
                    nNeighbor++;
                }
                else 
                {
                    if( ExistTriangle(notused, tj, tk)) return true;
                }
            }
        }
        for(i=0; i<nNeighbor; i++)
            for(int j=i+1; j<nNeighbor; j++)
            {
                if( Neighbor[i][1]==Neighbor[j][1])    continue;
                if( ExistTriangle( Neighbor[i][0], Neighbor[j][0], ti))
                {
                    if(ExistTriangle( Neighbor[i][0], Neighbor[j][0], tj))
                        if( ExistTriangle( Neighbor[i][0], tj, tk)||ExistTriangle(Neighbor[j][0], tj, tk))
                            return true;
                    if(ExistTriangle( Neighbor[i][0], Neighbor[j][0], tk))
                        if( ExistTriangle( Neighbor[i][0], tj, tk)||ExistTriangle(Neighbor[j][0], tj, tk))
                            return true;
                }
            }
        return false;
    }

    bool ExistTriangle(int ti, int tj, int tk, CYCLIC_ARRAY <DK_VERTEX <double> *> &maps)
    {
        ti=maps(ti)->index;
        tj=maps(tj)->index;
        tk=maps(tk)->index;
        //test if the triangle ti, tj, tk has already existed in the current mesh.
        for( int i=1; i<=(*mesh->incident_triangles)(ti).m; i++)
        {
            int triangle=(*mesh->incident_triangles)(ti)(i);
            int tti=mesh->triangles(1, triangle);
            int ttj=mesh->triangles(2, triangle);
            int ttk=mesh->triangles(3, triangle);
            int notused;
            if( SameTriples(ti, tj, tk, tti, ttj, ttk, notused)==3)    return true;
        }
        return false;
    }
    bool ExistTriangle(int ti, int tj, int tk)
    {
        //test if the triangle ti, tj, tk has already existed in the current mesh.
        for( int i=1; i<=(*mesh->incident_triangles)(ti).m; i++)
        {
            int triangle=(*mesh->incident_triangles)(ti)(i);
            int tti=mesh->triangles(1, triangle);
            int ttj=mesh->triangles(2, triangle);
            int ttk=mesh->triangles(3, triangle);
            int notused;
            if( SameTriples(ti, tj, tk, tti, ttj, ttk, notused)==3)    return true;
        }
        return false;
    }

};

template <class T0, class T1>
class DK_VERTEX_LIST
{
    T1 MAX_area;
    T1 MAX_curvature;
    SOLIDS_PARTICLES<T0,VECTOR_3D<T0> > *particles;
    DK_TRIANGLE_LIST<T0, T1> *tl;
    TRIANGLE_MESH *mesh;
public:    
    ARRAY <DK_VERTEX <T1> *> *vertex_sequence;
    int nVertices;
    DK_VERTEX<T1> *vertices;
    DK_VERTEX<T1> * operator()(int i) { assert( i>=1 && i<=nVertices); return &vertices[i]; }
    DK_VERTEX_LIST(    SOLIDS_PARTICLES<T0,VECTOR_3D<T0> > *emb_boundary_particles) 
    {
        particles=emb_boundary_particles;
        nVertices=particles->number;
        vertices=new DK_VERTEX <T1> [nVertices+1];
        for( int i=1; i<=nVertices; i++)
        {
            vertices[i].index=i;
            vertices[i].real_u=-1;
            vertices[i].real_v=-1;
            vertices[i].tag=false;
        }
        vertex_sequence=NULL;
    }
    ~DK_VERTEX_LIST(){ delete []vertices; if( vertex_sequence!=NULL) delete vertex_sequence;}
    void refresh(DK_TRIANGLE_LIST<T0, T1> *tlist, TRIANGLE_MESH *currentmesh)
    {
        tl=tlist;
        mesh=currentmesh;
        MAX_area=0;
        MAX_curvature=0;
        //First, compute the area and the curvature of each vertex.
        for(int i=1; i<=nVertices; i++)
            if(vertices[i].flag!=FLAG_REMOVED)
            {
                vertices[i].flag=FLAG_ACTIVE;
                //compute the area and curvature.
                double mixed_area=0;
                vertices[i].area=0;
                VECTOR mean_K=VECTOR(0, 0, 0);
                double Gaussian_K=0;
                for(int j=1; j<=(*mesh->incident_triangles)(i).m; j++)
                {
                    int t=(*mesh->incident_triangles)(i)(j);
                    DK_TRIANGLE<T1> *info=tl->GetInfo(t);
                    vertices[i].area+=info->area;

                    int tj, tk;
                    double ctan_i, ctan_j, ctan_k;
                    info->Arrange(i, tj, tk, ctan_i, ctan_j, ctan_k);
                    VECTOR e1=particles->X(tj)-particles->X(i);
                    VECTOR e2=particles->X(tk)-particles->X(i);
                    
                    if( ctan_i > 0 && ctan_j > 0 && ctan_k > 0) //voronoi area used.
                        mixed_area+=(e1.Magnitude_Squared()*ctan_k + e2.Magnitude_Squared()*ctan_j)/8.0f;
                    else if(ctan_i < 0 )
                        mixed_area+=info->area/2;
                    else
                        mixed_area+=info->area/4;
                    
                    mean_K+=e1*ctan_k+e2*ctan_j;
                    Gaussian_K+=VECTOR::Angle_Between(e1, e2);
                }
                //decide the curvature computing method.
                vertices[i].curvature=mean_K.Magnitude()/mixed_area;
                //vertices[i].curvature=fabsf(2*Pi-Gaussian_K)/mixed_area;
                if( MAX_area < vertices[i].area)    MAX_area=vertices[i].area;
                if( MAX_curvature < vertices[i].curvature) MAX_curvature=vertices[i].curvature;
            }
        //Then, set the weight.
        if( vertex_sequence!=NULL)    delete vertex_sequence;
        vertex_sequence=new ARRAY <DK_VERTEX <double> *>;
        for( i=1; i<=nVertices; i++)
            if( vertices[i].flag==FLAG_ACTIVE )
            {
                vertices[i].weight=
                    vertices[i].area/MAX_area+vertices[i].curvature/MAX_curvature;
                vertex_sequence->Append(&vertices[i]);
            }
        //Finally, sort the sequence.
        QuickSort(*vertex_sequence, 1, vertex_sequence->m);
    }
    void remove(int index) 
    {
        vertices[index].flag=FLAG_REMOVED; 
        for(int i=1; i<= (*mesh->neighbor_nodes)(index).m; i++)
        {
            if( vertices[(*mesh->neighbor_nodes)(index)(i)].flag==FLAG_REMOVED)
                printf("ERROR: removed vertex %d appears again!\n", (*mesh->neighbor_nodes)(index)(i));
            vertices[(*mesh->neighbor_nodes)(index)(i)].flag=FLAG_DEAD;
        }
    }
    int Partition(ARRAY <DK_VERTEX <double> *>&A, int p, int r)
    {
        double x=A(r)->weight;
        int i=p-1;
        for( int j=p; j<=r-1; j++)
        {
            if( A(j)->weight <=x )
            {
                i++;
                DK_VERTEX <double> *temp=A(i);
                A(i)=A(j);
                A(j)=temp;
            }
        }
        DK_VERTEX <double> *temp=A(i+1);
        A(i+1)=A(r);
        A(r)=temp;
        return i+1;
    }
    void QuickSort(ARRAY <DK_VERTEX <double> *>&A, int p, int r)
    {
        if( p < r)
        {
            int q=Partition(A, p, r);
            QuickSort(A, p, q-1);
            QuickSort(A, q+1, r);
        }
    }
};

#endif
