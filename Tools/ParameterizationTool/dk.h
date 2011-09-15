#ifndef DK_H
#define DK_H

#include <cstring>
#include <fstream>
#include <iostream> // writing to screen: cin, cout
#include <math.h> // standard math library
//#include <stdlib>

#include "../../Public_Library/Geometry/Triangulated_Surface.h"
#include "../../Public_Library/Grids/Segment_Mesh.h"
#include "utils.h"
using namespace PhysBAM;
using std::fstream;
using std::ios;


class DK
{
public:
    TRIANGLE_MESH backup_boundary_mesh;
    TRIANGLE_MESH emb_boundary_mesh;
    SOLIDS_PARTICLES<float,VECTOR_3D<float> > emb_boundary_particles;
    TRIANGULATED_SURFACE<float> emb_boundary_surf;
    TRIANGLE_MESH *currentmesh;
    DK_VERTEX_LIST <float, double> *vlist;
    DK_TRIANGLE_LIST <float, double> *tlist;
    ARRAY < ARRAY <int > > vertices_info;
    double aveX, aveY, aveZ;

public:
    DK():emb_boundary_surf(emb_boundary_mesh,emb_boundary_particles)
    {
        aveX=aveY=aveZ=0;
        tlist=NULL;
    }
    ~DK()
    {
        if(tlist!=NULL)        delete tlist;
        if(vlist!=NULL)        delete vlist;
    }

    void Open( char * filename)
    {
        fstream input;
        input.open(filename,ios::in|ios::binary);
        assert(input.is_open());
        emb_boundary_surf.Clean_Up_Memory();
        emb_boundary_surf.Read<float>(input);   

        for(int t=1; t<=emb_boundary_surf.particles.number; t++)
        {
              VECTOR_3D<double> P=emb_boundary_particles.X(t);
            aveX+=P.x/emb_boundary_particles.number;
            aveY+=P.y/emb_boundary_particles.number;
            aveZ+=P.z/emb_boundary_particles.number;
        }
        vlist=new DK_VERTEX_LIST <float, double> (&emb_boundary_particles);
        currentmesh=&emb_boundary_mesh;
        backup_boundary_mesh=emb_boundary_mesh;
        backup_boundary_mesh.Initialize_Topologically_Sorted_Neighbor_Nodes();
        backup_boundary_mesh.Initialize_Incident_Triangles();
        printf("Open file %s successfully.\n", filename);
        printf("%d particles and %d triangles.\n", emb_boundary_particles.number,
            currentmesh->triangles.m);
    }

    void Init_Run()
    {
        currentmesh->Initialize_Topologically_Sorted_Neighbor_Nodes();
        currentmesh->Initialize_Incident_Triangles();
        if(tlist!=NULL)    delete tlist;
        tlist=new DK_TRIANGLE_LIST <float, double> (&emb_boundary_particles, currentmesh, vlist, vertices_info);
        vlist->refresh(tlist, currentmesh);
    }
    void Run()
    {
        Init_Run();    
        for(int i=1; i<= (*vlist->vertex_sequence).m; i++)
        {
            int index=(*vlist->vertex_sequence)(i)->index;
            if( vlist->vertices[index].flag==FLAG_ACTIVE) 
            {
                vlist->remove(index);
                Build(index, *vlist);
            }
        }

        ARRAYS<int> tlistarray(3,1,0);
        vertices_info.Clean_Up_Memory();
        tlist->Getlist( tlistarray , vertices_info);
        if( currentmesh!=NULL && currentmesh!=&emb_boundary_mesh)
            delete currentmesh;
        currentmesh=new TRIANGLE_MESH(tlistarray);
    }
    void flip(int t0, int t1, int p0, int p1, int p2, int p3)
    {
        VECTOR_3D <double> e1=emb_boundary_particles.X(p1)-emb_boundary_particles.X(p0);
        VECTOR_3D <double> e2=emb_boundary_particles.X(p2)-emb_boundary_particles.X(p0);
        VECTOR_3D <double> e3=emb_boundary_particles.X(p3)-emb_boundary_particles.X(p0);

        VECTOR_3D <double> P0=VECTOR_3D<double> (0, 0, 0);
        VECTOR_3D <double> P1=VECTOR_3D<double> (e1.Magnitude(), 0, 0);
        double x_costheta1=VECTOR_3D<double>::Dot_Product(e2, e1)/e1.Magnitude();
        VECTOR_3D <double> P2=VECTOR_3D<double> (x_costheta1, sqrt(e2.Magnitude_Squared()-
            x_costheta1*x_costheta1), 0);
        double x_costheta2=VECTOR_3D<double>::Dot_Product(e3, e1)/e1.Magnitude();
        VECTOR_3D <double> P3=VECTOR_3D<double> (x_costheta2, -sqrt(e3.Magnitude_Squared()-
            x_costheta2*x_costheta2), 0);
        
        ARRAY <DK_VERTEX <double> *> free_vlist;
        for( int k=1; k<=vertices_info(t0).m; k++)
        {
            double u=(*vlist)(vertices_info(t0)(k))->u;
            double v=(*vlist)(vertices_info(t0)(k))->v;
            VECTOR_3D<double> PP0, PP1, PP2;
            if ( p0==currentmesh->triangles(1, t0))
            {
                PP0=P0; PP1=P1; PP2=P2;
            }
            if ( p1==currentmesh->triangles(1, t0))
            {
                PP0=P1; PP1=P2; PP2=P0;
            }
            if ( p2==currentmesh->triangles(1, t0))
            {
                PP0=P2; PP1=P0; PP2=P1;
            }
            (*vlist)(vertices_info(t0)(k))->temp_P=PP0*u+PP1*v+PP2*(1-u-v);
            free_vlist.Append( (*vlist)( vertices_info(t0)(k) ) );
        }
        for(  k=1; k<=vertices_info(t1).m; k++)
        {
            double u=(*vlist)(vertices_info(t1)(k))->u;
            double v=(*vlist)(vertices_info(t1)(k))->v;
            VECTOR_3D<double> PP0, PP1, PP2;
            if ( p0==currentmesh->triangles(1, t1))
            {
                PP0=P0; PP1=P3; PP2=P1;
            }
            if ( p3==currentmesh->triangles(1, t1))
            {
                PP0=P3; PP1=P1; PP2=P0;
            }
            if ( p1==currentmesh->triangles(1, t1))
            {
                PP0=P1; PP1=P0; PP2=P3;
            }
            (*vlist)(vertices_info(t1)(k))->temp_P=PP0*u+PP1*v+PP2*(1-u-v);
            free_vlist.Append( (*vlist)( vertices_info(t1)(k) ) );
        }
        
        currentmesh->triangles(1, t0)=p0;
        currentmesh->triangles(2, t0)=p3;
        currentmesh->triangles(3, t0)=p2;
        currentmesh->triangles(1, t1)=p1;
        currentmesh->triangles(2, t1)=p2;
        currentmesh->triangles(3, t1)=p3;

        vertices_info(t0).Clean_Up_Memory();
        vertices_info(t1).Clean_Up_Memory();

        for(int i=1; i<=free_vlist.m; i++)
        {
            double u, v;
            if( Compute_Barycentric(P0, P3, P2, free_vlist(i)->temp_P, u, v)==true)
            {
                free_vlist(i)->u=u;
                free_vlist(i)->v=v;
                vertices_info(t0).Append( free_vlist(i)->index );
                continue;
            }
            else if( Compute_Barycentric(P1, P2, P3, free_vlist(i)->temp_P, u, v)==true)
            {
                free_vlist(i)->u=u;
                free_vlist(i)->v=v;
                vertices_info(t1).Append( free_vlist(i)->index );
                continue;
            }
            else printf("error when flipping\n");
        }
            
    
    }

    void Build(int p0, DK_VERTEX_LIST<float, double> &vlist);
    VECTOR_3D<double> GetPosition(int p0, CYCLIC_ARRAY <DK_VERTEX <double> *> &maps, int target);
    double GetAngle(int center, int i, int j, CYCLIC_ARRAY<DK_VERTEX <double> *> &maps, int p0);
    bool Is_Covered( CYCLIC_ARRAY <DK_VERTEX <double> *> &maps, int except);
    bool Output(char *filename);
};


#endif