#pragma once
#include <vector>
#include "dk.h"
#include <gl/glut.h>
#include <math.h>
using std::vector;
#define PI 3.1415926
struct Vertex
{
    int vi;
    double x, y, z;    //position
    double vx, vy, vz; //velocity, let assume no velocity accumulation ( infinite damping here)
    double ax, ay, az; // Acceleration
    int count;
    int index;
    void Normalize() { double invlength=1/sqrt(x*x+y*y+z*z); x=x*invlength; y=y*invlength; z=z*invlength;}
};
struct Edge
{
    int v0, v1;
    double InitialLength;

};
class SphericalPM
{
    vector <Edge> edges;
    vector <Vertex> vertices;
    int sv0,sv1;
public:
    SphericalPM(void);
    ~SphericalPM(void);
    void Draw(double swingangle, double elevateangle);
    void Initialze(DK &dk, int stv0, int stv1)
    {
        edges.clear();
        vertices.clear();
        dk.currentmesh->Initialize_Segment_Mesh();
        for( int i=1; i<=dk.currentmesh->segment_mesh->segments.m; i++)
        {
            int vv0 = dk.currentmesh->segment_mesh->segments(1, i);
            int vv1 = dk.currentmesh->segment_mesh->segments(2, i);
            int v0=AddVertex(dk.emb_boundary_particles.X(vv0).x, dk.emb_boundary_particles.X(vv0).y, dk.emb_boundary_particles.X(vv0).z, vv0);
            int v1=AddVertex(dk.emb_boundary_particles.X(vv1).x, dk.emb_boundary_particles.X(vv1).y, dk.emb_boundary_particles.X(vv1).z, vv1);
    
            Edge e;
            e.v0=v0;
            e.v1=v1;
            e.InitialLength=Distance(v0, v1)*5; //set up the inital length of the edge spring
            edges.push_back(e);
        }

        sv0=sv1=-1;
        if( stv0!=0)
        {
            stv0=(*dk.vlist->vertex_sequence)(stv0)->index;
            sv0=FindVertex(stv0);
        //    vertices[sv0].x=0;
        //    vertices[sv0].y=1;
        //    vertices[sv0].z=0;
        }
        if( stv1!=0)
        {
            stv1=(*dk.vlist->vertex_sequence)(stv1)->index;
            sv1=FindVertex(stv1);
        }
        printf("%d, %d\n", sv0, sv1);

        Centralize();
        Normalize();


        double backupx=vertices[sv0].x;
        double backupy=vertices[sv0].y;
        double backupz=vertices[sv0].z;

        double a, b, c, d;
        a= cos( acos(backupy)/2 );
        double bb=-backupz;
        double dd=backupx;
        b=bb * sqrt((1-a*a)/(bb*bb+dd*dd));
        c=0;
        d=dd * sqrt((1-a*a)/(bb*bb+dd*dd));
    //    printf("before %f, %f, %f\n", vertices[staticv].x, vertices[staticv].y, vertices[staticv].z);
        double M[3][3];
        M[0][0]=a*a+b*b-c*c-d*d;
        M[0][1]=2*b*c-2*a*d;
        M[0][2]=2*b*d+2*a*c;
        M[1][0]=2*b*c+2*a*d;
        M[1][1]=a*a-b*b+c*c-d*d;
        M[1][2]=2*c*d-2*a*b;
        M[2][0]=2*b*d-2*a*c;
        M[2][1]=2*c*d+2*a*b;
        M[2][2]=a*a-b*b-c*c+d*d;

        for(int i=0; i<vertices.size(); i++)
        {
            double xx=vertices[i].x;
            double yy=vertices[i].y;
            double zz=vertices[i].z;

            vertices[i].x=M[0][0]*xx+M[0][1]*yy+M[0][2]*zz;
            vertices[i].y=M[1][0]*xx+M[1][1]*yy+M[1][2]*zz;
            vertices[i].z=M[2][0]*xx+M[2][1]*yy+M[2][2]*zz;
        }

        if(sv1!=-1)
        {
            vertices[sv1].x=0;
            vertices[sv1].y=-1;
            vertices[sv1].z=0;
        }

        //printf("after %f, %f, %f\n", vertices[staticv].x, vertices[staticv].y, vertices[staticv].z);
    }
    void Centralize()
    {
        double sumx=0;
        double sumy=0;
        double sumz=0;
        for( int i=0; i<vertices.size(); i++)
        {
            sumx+=vertices[i].x;
            sumy+=vertices[i].y;
            sumz+=vertices[i].z;
        }
        sumx/=vertices.size();
        sumy/=vertices.size();
        sumz/=vertices.size();
        for( int i=0; i<vertices.size(); i++)
        {
            vertices[i].x-=sumx;
            vertices[i].y-=sumy;
            vertices[i].z-=sumz;
        }
    }
    void Normalize()
    {
        for( int i=0; i<vertices.size(); i++)
            vertices[i].Normalize();
    }
    void Interations(int i)
    {
        for(;i>0; i--)
            Iteration();
    }
    void Iteration()
    {
        static int tt=0;
        printf("iteration %d\n", tt);
        tt++;
        for( int i=0; i<edges.size(); i++)
        {
            vertices[edges[i].v0].ax=0;
            vertices[edges[i].v0].ay=0;
            vertices[edges[i].v0].az=0;
            vertices[edges[i].v1].ax=0;
            vertices[edges[i].v1].ay=0;
            vertices[edges[i].v1].az=0;
            vertices[edges[i].v0].count=0;
            vertices[edges[i].v1].count=0;
        }
        for( i=0; i<edges.size(); i++)
        {
            /*double f=(Distance(edges[i].v0, edges[i].v1)-edges[i].InitialLength)/Distance(edges[i].v0, edges[i].v1);
            double angle=    vertices[edges[i].v0].x * vertices[edges[i].v1].x+
                            vertices[edges[i].v0].y * vertices[edges[i].v1].y+
                            vertices[edges[i].v0].z * vertices[edges[i].v1].z;
            angle=acos(angle)/2;
            f=f*cos(angle);

            vertices[edges[i].v0].ax  = vertices[edges[i].v0].ax + f * ( vertices[edges[i].v1].x - vertices[edges[i].v0].x);
            vertices[edges[i].v0].ay  = vertices[edges[i].v0].ay + f * ( vertices[edges[i].v1].y - vertices[edges[i].v0].y);
            vertices[edges[i].v0].az  = vertices[edges[i].v0].az + f * ( vertices[edges[i].v1].z - vertices[edges[i].v0].z);
        
            vertices[edges[i].v1].ax  = vertices[edges[i].v1].ax + f * ( vertices[edges[i].v0].x - vertices[edges[i].v1].x);
            vertices[edges[i].v1].ay  = vertices[edges[i].v1].ay + f * ( vertices[edges[i].v0].y - vertices[edges[i].v1].y);
            vertices[edges[i].v1].az  = vertices[edges[i].v1].az + f * ( vertices[edges[i].v0].z - vertices[edges[i].v1].z);*/
            vertices[edges[i].v0].ax = vertices[edges[i].v0].ax + vertices[edges[i].v1].x;
            vertices[edges[i].v0].ay = vertices[edges[i].v0].ay + vertices[edges[i].v1].y;
            vertices[edges[i].v0].az = vertices[edges[i].v0].az + vertices[edges[i].v1].z;
            vertices[edges[i].v0].count ++;
            vertices[edges[i].v1].ax = vertices[edges[i].v1].ax + vertices[edges[i].v0].x;
            vertices[edges[i].v1].ay = vertices[edges[i].v1].ay + vertices[edges[i].v0].y;
            vertices[edges[i].v1].az = vertices[edges[i].v1].az + vertices[edges[i].v0].z;
            vertices[edges[i].v1].count ++;
        }
        //double backupx=vertices[staticv].x;
        //double backupy=vertices[staticv].y;
        //double backupz=vertices[staticv].z;
        for( i=0; i<vertices.size(); i++)
        {
            //vertices[i].x += vertices[i].ax *0.1;
            //vertices[i].y += vertices[i].ay *0.1;
            //vertices[i].z += vertices[i].az *0.1;
            vertices[i].x = vertices[i].ax / vertices[i].count;
            vertices[i].y = vertices[i].ay / vertices[i].count;
            vertices[i].z = vertices[i].az / vertices[i].count;
            //vertices[i].x = vertices[i].x/2 + vertices[i].ax / (2 * vertices[i].count);
            //vertices[i].y = vertices[i].y/2 + vertices[i].ay / (2 * vertices[i].count);
            //vertices[i].z = vertices[i].z/2 + vertices[i].az / (2 * vertices[i].count);
            //vertices[i].Normalize();
        }
        Centralize();
        Normalize();
        if( sv0!=-1)
        {
            vertices[sv0].x=0;
            vertices[sv0].y=1;
            vertices[sv0].z=0;
        }
        if( sv1!=-1)
        {
            vertices[sv1].x=0;
            vertices[sv1].y=-1;
            vertices[sv1].z=0;
        }

    }
    void Output(DK &dk)
    {
        ARRAY <DK_VERTEX <double> *> &vertex_sequence=*dk.vlist->vertex_sequence;
        for(int i=0;i<vertex_sequence.m;i++)
        {
            int newindex=FindVertex(vertex_sequence(i)->index);
            //printf("newindex %d\n", newindex);
            if(newindex<0)    printf("unexpected ERROR???\n");
            double x=vertices[newindex].x;
            double y=vertices[newindex].y;
            double z=vertices[newindex].z;

            double d=sqrtf(vertices[newindex].x*vertices[newindex].x + vertices[newindex].z*vertices[newindex].z);
            double u, v;
            //compute u
            if( x==0 )
            {
                if( z>0) u=PI/2;
                if( z<0) u=-PI/2;
            }
            else u=acos( x/sqrt(z*z+x*x));
            if( z<0 )    u=2*PI-u;

            //compute v
            if( d==0 )
            {
                if( y>=0 )    v=PI/2;
                if( y<0 )    v=-PI/2;
            }
            else v=atan(y/d);
            v+=PI/2;
            vertex_sequence(i)->real_u=u/(2*PI);
            vertex_sequence(i)->real_v=v/PI;
        }
    }
    int AddVertex(double x, double y, double z, int vi)
    {
        for( int i=0; i<vertices.size(); i++)
            if( vertices[i].vi==vi)    return i;
        Vertex v;
        v.vi=vi;
        v.x=x; v.y=y; v.z=z;
        v.vx=v.vy=v.vz=0;
        v.ax=v.ay=v.az=0;
        vertices.push_back(v);
        return i;
    }
    int FindVertex(int vi)
    {
        for(int i=0; i<vertices.size(); i++)
            if( vertices[i].vi==vi)    return i;
        return -1;
    }
    double Distance(int v0, int v1)
    {
        double x=vertices[v0].x-vertices[v1].x;
        double y=vertices[v0].y-vertices[v1].y;
        double z=vertices[v0].z-vertices[v1].z;
        return sqrt(x*x+y*y+z*z);
    }

};
