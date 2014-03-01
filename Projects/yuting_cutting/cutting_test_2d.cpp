//  Created by Yuting Wang on 2/25/14.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.

#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <X11/Xlib.h>

#ifdef __APPLE__//Mac OS
//#  include <GL/glew.h>
#  include <OpenGL/OpenGL.h>
#  include <GLUT/glut.h>
#else
//#  include <GL/glew.h>
#  include <GL/freeglut.h>
#  include <GL/freeglut_ext.h>
#endif  // __APPLE__

#include "CUTTING_2D.h"
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Tools/Data_Structures/HASHTABLE.h>

using namespace PhysBAM;
using namespace std;

typedef double T;
typedef VECTOR<T, 2> TV;
typedef VECTOR<int,2> I2;
typedef VECTOR<int,3> I3;

//global variables
int argc1;
char **argv1;
int window_width=600;
int window_height=600;
TV starting_position;
T trans_speed=0.1;
T scale_speed=1.1;
TRIANGULATED_AREA<T>* ta=NULL;
SEGMENTED_CURVE_2D<T>* sc=NULL;
ARRAY<int> labels;
HASHTABLE<int> dragging_particles;
bool dragging=false,cutting=false;

void Run_Cutter()
{
    CUTTING<TV> cutter(*ta,*sc);
    if(sc->mesh.elements.m>0){
        cutter.Run(.01);
        ta->mesh.Identify_Connected_Components(labels);
        cutting=false;
        
        //reinitialize cutter->ta
        TRIANGULATED_AREA<T>* nta=TRIANGULATED_AREA<T>::Create();
        nta->particles.Resize(ta->particles.X.m);
        for(int i=0;i<ta->particles.X.m;++i)
            nta->particles.X(i)=ta->particles.X(i);
        nta->mesh.elements=ta->mesh.elements;
        nta->Update_Number_Nodes();
        delete ta;
        ta=nta;
    }
    delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
}

static void Key(unsigned char key, int x, int y)
{
    switch( key ) {
        case 033: // Escape Key
            exit(EXIT_SUCCESS);
            break;
        case 'w':
            for(int i=0;i<ta->particles.X.m;++i)
                ta->particles.X(i)*=scale_speed;
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)*=scale_speed;
            glutPostRedisplay();
            break;
        case 's':
            for(int i=0;i<ta->particles.X.m;++i)
                ta->particles.X(i)/=scale_speed;
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)/=scale_speed;
            glutPostRedisplay();
            break;
    }
}

void Mouse(int button, int state, int x, int y)
{
    TV location(2*x/T(window_width)-1,1-2*y/T(window_height));
    if(button==4){
        for(int i=0;i<ta->particles.X.m;++i)
            ta->particles.X(i)/=scale_speed;
        for(int i=0;i<sc->particles.X.m;++i)
            sc->particles.X(i)/=scale_speed;
        glutPostRedisplay();
        return;
    }
    else if(button==3){
        for(int i=0;i<ta->particles.X.m;++i)
            ta->particles.X(i)*=scale_speed;
        for(int i=0;i<sc->particles.X.m;++i)
            sc->particles.X(i)*=scale_speed;
        glutPostRedisplay();
        return;
    }
    if(state==GLUT_DOWN){
        if(button==GLUT_RIGHT_BUTTON){
            dragging=true;
            cutting=false;
            starting_position=location;
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            glBegin(GL_TRIANGLES);
            for(int t=0;t<ta->mesh.elements.m;t++){
                I3 tri=ta->mesh.elements(t);
                glColor4d(labels(t)/255.,0,0,0);
                for(int j=0;j<3;++j)
                    glVertex2d(ta->particles.X(tri(j))(0),ta->particles.X(tri(j))(1));
            }
            glEnd();
            unsigned char val;
            glReadPixels(x, window_height-y, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &val);
            dragging_particles.Clean_Memory();
            for(int t=0;t<ta->mesh.elements.m;t++)
                if(labels(t)==(int)val)
                    for(int i=0;i<3;++i)
                        dragging_particles.Set(ta->mesh.elements(t)(i));
        }
        else if(button==GLUT_LEFT_BUTTON){
            cutting=true;
            dragging=false;
            int p=sc->particles.Add_Element();
            sc->Update_Number_Nodes();
            sc->particles.X(p)=location;
        }
    }
    else if(state==GLUT_UP){
        if(button==GLUT_LEFT_BUTTON){
            Run_Cutter();
            glutPostRedisplay();
        }
        else if(button==GLUT_RIGHT_BUTTON)
            dragging=false;
    }
}

void Motion(int x, int y)
{
    TV location(2*x/T(window_width)-1,1-2*y/T(window_height));
    if (cutting){
        int p=sc->particles.Add_Element();
        sc->particles.X(p)=location;
        sc->mesh.elements.Append(I2(p-1,p));
        sc->Update_Number_Nodes();
    }
    else if(dragging){
        for(typename HASHTABLE<int>::ITERATOR it(dragging_particles);it.Valid();it.Next())
            ta->particles.X(it.Key())+=(location-starting_position);
        starting_position=location;
    }
    glutPostRedisplay();
}

void Render(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnableClientState(GL_VERTEX_ARRAY);
    
    ARRAY<TV> vertices;
    for(int t=0;t<sc->mesh.elements.m;t++){
        I2 edge=sc->mesh.elements(t);
        for(int i=0;i<2;++i)
            vertices.Append(sc->particles.X(edge(i)));
    }
    glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
    glColor4f(0.0, 0.0, 1.0, 1.0);
    glDrawArrays(GL_LINES,0,vertices.m);
    
    vertices.Remove_All();
    //    ta->mesh.Initialize_Boundary_Mesh();
    //    for(int t=0;t<ta->mesh.boundary_mesh->elements.m;t++){
    //        I2 edge=ta->mesh.boundary_mesh->elements(t);
    //        for(int i=0;i<2;++i)
    //            vertices.Append(ta->particles.X(edge(i)));
    //    }
    for(int t=0;t<ta->mesh.elements.m;t++){
        I3 tri=ta->mesh.elements(t);
        for(int i=0;i<3;++i){
            vertices.Append(ta->particles.X(tri(i)));
            vertices.Append(ta->particles.X(tri((i+1)%3)));
        }
    }
    glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
    glColor4f(0.0, 1.0, 0.0, 1.0);
    glDrawArrays(GL_LINES,0,vertices.m);
    
    vertices.Remove_All();
    for(int t=0;t<ta->mesh.elements.m;t++){
        I3 tri=ta->mesh.elements(t);
        for(int i=0;i<3;++i)
            vertices.Append(ta->particles.X(tri(i)));
    }
    glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
    glColor4f(1.0, 0.0, 0.0, 0.0);
    glDrawArrays(GL_TRIANGLES,0,vertices.m);
    
    glutSwapBuffers();
    glDisableClientState(GL_VERTEX_ARRAY);
}

void Reshape(GLint newWidth,GLint newHeight) {
    glViewport(0,0,newWidth,newHeight);
    window_width=newWidth;
    window_height=newHeight;
    glutPostRedisplay();
}

void Initialize_Meshes()
{
    //ta=TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(),.5),20);
    ta=TRIANGULATED_AREA<T>::Create();
    ta->particles.Add_Elements(4);
    ta->particles.X(0)=TV(0,.5);
    ta->particles.X(1)=TV(-.5,0);
    ta->particles.X(2)=TV(.5,0);
    ta->particles.X(3)=TV(0,-.5);
    ta->mesh.elements.Append(I3(0,1,2));
    ta->mesh.elements.Append(I3(2,1,3));
    ta->Update_Number_Nodes();
    
    ta->mesh.Identify_Connected_Components(labels);
    
    sc=SEGMENTED_CURVE_2D<T>::Create();
    
    sc->particles.Add_Elements(3);
    sc->particles.X(0)=TV(-.8,0);
    sc->particles.X(1)=TV(.1,0);
    sc->particles.X(2)=TV(.8,0);
    sc->mesh.elements.Append(I2(0,1));
    sc->mesh.elements.Append(I2(1,2));
    sc->Update_Number_Nodes();
    
    Run_Cutter();
    
    cout << "initialized mesh\n";
}

int main(int argc, char **argv)
{
    argc1 = argc;
    argv1 = argv;
    glutInit( &argc, argv );
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
    glutInitWindowSize(window_width,window_height);
    string window_name="cutting";
    glutCreateWindow(window_name.c_str());
    glClearColor(0.0,0.0,0.0,1.0);
    glEnable(GL_DEPTH_TEST);
    
    Initialize_Meshes();
    
    //glutSpecialFunc(Special_Key);
    glutKeyboardFunc(Key);
    glutMouseFunc(Mouse);
    glutMotionFunc(Motion);
    glutDisplayFunc(Render);
    glutReshapeFunc(Reshape);
    
    glutMainLoop();
    
    return 0;
}
