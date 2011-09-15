#include ".\sphericalpm.h"

SphericalPM::SphericalPM(void)
{
}

SphericalPM::~SphericalPM(void)
{
}

void SphericalPM::Draw(double swingangle, double elevateangle)
{
    gluLookAt(0, 0, 20, 0, 0, 0, 0, 1, 0);
    glRotated(swingangle, 1, 0, 0);
    glRotated(elevateangle, 0, 1, 0);
    glColor3f(0,0,1);
    glutWireSphere(1, 20, 20);

    glColor3f(1,1, 1);
    for(int i=0; i<edges.size(); i++)
    {
        int v0=edges[i].v0;
        int v1=edges[i].v1;
        glBegin(GL_LINES);
        glVertex3d( vertices[v0].x, vertices[v0].y, vertices[v0].z);
        glVertex3d( vertices[v1].x, vertices[v1].y, vertices[v1].z);
    //    printf("%d %d: %f, %f, %f\n", v0, v1, vertices[v0].x, vertices[v0].y, vertices[v0].z);
        glEnd();
    }
    for(int i=0; i<vertices.size(); i++)
    {
        if( sv0==i || sv1==i)
        {
            glColor3f(1, 1, 0);
            glPushMatrix();
            glTranslatef(vertices[i].x, vertices[i].y, vertices[i].z);
            glutWireSphere(0.05, 20, 20);
            glPopMatrix();
        }
    }

    glutSwapBuffers();
}
