//  Created by Yuting Wang on 2/25/14.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <X11/Xlib.h>

#ifdef __APPLE__//Mac OS
//#include <GL/glew.h>
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#else
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
//#include <GL/glew.h>
#endif  // __APPLE__

#include <Tools/Data_Structures/HASHTABLE.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include "CUTTING_2D.h"

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
TRIANGULATED_AREA<T>* sim_ta=NULL;
SEGMENTED_CURVE_2D<T>* sc=NULL;
CUTTING<TV>* cutter;
ARRAY<int> labels;
HASHTABLE<int> dragging_particles;
bool dragging=false,cutting=false,draw_sim=true,draw_cutting_edge=true;

void Run_Cutter()
{
    if(sc->mesh.elements.m>0){
        cutter->Run(.01);
        cutter->ta->mesh.Identify_Connected_Components(labels);
        
        cutting=false;
        
        //reinitialize ta
        TRIANGULATED_AREA<T>* nta=TRIANGULATED_AREA<T>::Create();
        nta->particles.Resize(cutter->ta->particles.X.m);
        for(int i=0;i<cutter->ta->particles.X.m;++i)
            nta->particles.X(i)=cutter->ta->particles.X(i);
        nta->mesh.elements=cutter->ta->mesh.elements;
        nta->Update_Number_Nodes();
        delete cutter->ta;
        cutter->ta=nta;
    }
}

static void Key(unsigned char key, int x, int y)
{
    switch( key ) {
        case 033: // Escape Key
            exit(EXIT_SUCCESS);
            break;
        case 'z':
            draw_sim=!draw_sim;
            break;
        case 'x':
            draw_cutting_edge=!draw_cutting_edge;
            break;
        case 'q':
            for(int i=0;i<sim_ta->particles.X.m;++i)
                sim_ta->particles.X(i)*=scale_speed;
            cutter->Update_Material_Particles();
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)*=scale_speed;
            break;
        case 'e':
            for(int i=0;i<sim_ta->particles.X.m;++i)
                sim_ta->particles.X(i)/=scale_speed;
            cutter->Update_Material_Particles();
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)/=scale_speed;
            break;
        case 'w':
            for(int i=0;i<sim_ta->particles.X.m;++i)
                sim_ta->particles.X(i)(1)+=trans_speed;
            cutter->Update_Material_Particles();
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)(1)+=trans_speed;
            break;
        case 's':
            for(int i=0;i<sim_ta->particles.X.m;++i)
                sim_ta->particles.X(i)(1)-=trans_speed;
            cutter->Update_Material_Particles();
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)(1)-=trans_speed;
            break;
        case 'a':
            for(int i=0;i<sim_ta->particles.X.m;++i)
                sim_ta->particles.X(i)(0)-=trans_speed;
            cutter->Update_Material_Particles();
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)(0)-=trans_speed;
            break;
        case 'd':
            for(int i=0;i<sim_ta->particles.X.m;++i)
                sim_ta->particles.X(i)(0)+=trans_speed;
            cutter->Update_Material_Particles();
            for(int i=0;i<sc->particles.X.m;++i)
                sc->particles.X(i)(0)+=trans_speed;
            break;
    }
    glutPostRedisplay();
}

void Mouse(int button, int state, int x, int y)
{
    TV location(2*x/T(window_width)-1,1-2*y/T(window_height));
    if(button==4){
        for(int i=0;i<sim_ta->particles.X.m;++i)
            sim_ta->particles.X(i)/=scale_speed;
        cutter->Update_Material_Particles();
        for(int i=0;i<sc->particles.X.m;++i)
            sc->particles.X(i)/=scale_speed;
        glutPostRedisplay();
        return;
    }
    else if(button==3){
        for(int i=0;i<sim_ta->particles.X.m;++i)
            sim_ta->particles.X(i)*=scale_speed;
        cutter->Update_Material_Particles();
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
            for(int t=0;t<cutter->ta->mesh.elements.m;t++){
                I3 tri=cutter->ta->mesh.elements(t);
                glColor4d(labels(t)/255.,0,0,0);
                for(int j=0;j<3;++j)
                    glVertex2d(cutter->ta->particles.X(tri(j))(0),cutter->ta->particles.X(tri(j))(1));
            }
            glEnd();
            unsigned char val;
            glReadPixels(x, window_height-y, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &val);
            dragging_particles.Clean_Memory();
            for(int t=0;t<cutter->ta->mesh.elements.m;t++)
                if(labels(t)==(int)val)
                    for(int i=0;i<3;++i){
                        int f=cutter->particle_in_sim(cutter->ta->mesh.elements(t)(i)).x;
                        for(int j=0;j<3;++j)
                            dragging_particles.Set(sim_ta->mesh.elements(f)(j));
                    }
        }
        else if(button==GLUT_LEFT_BUTTON){
            cutting=true;
            dragging=false;
            //reinitialize cutting curve
            delete sc;
            sc=SEGMENTED_CURVE_2D<T>::Create();
            cutter->sc=sc;
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
            sim_ta->particles.X(it.Key())+=(location-starting_position);
        cutter->Update_Material_Particles();
        starting_position=location;
    }
    glutPostRedisplay();
}

void Render(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnableClientState(GL_VERTEX_ARRAY);
    
    //cutting curve
    ARRAY<TV> vertices;
    for(int t=0;t<sc->mesh.elements.m;t++){
        I2 edge=sc->mesh.elements(t);
        for(int i=0;i<2;++i)
            vertices.Append(sc->particles.X(edge(i)));
    }
    glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
    glColor4f(1.0, 1.0, 1.0, 1.0);
    glDrawArrays(GL_LINES,0,vertices.m);
    
    //edges of sim_ta
    if(draw_sim){
        vertices.Remove_All();
        for(int t=0;t<sim_ta->mesh.elements.m;t++){
            I3 tri=sim_ta->mesh.elements(t);
            for(int i=0;i<3;++i){
                vertices.Append(sim_ta->particles.X(tri(i)));
                vertices.Append(sim_ta->particles.X(tri((i+1)%3)));
            }
        }
        glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
        glColor4f(0.0, 0.0, 1.0, 1.0);
        glDrawArrays(GL_LINES,0,vertices.m);
    }
    
    //edges of ta
    if(draw_cutting_edge){
        vertices.Remove_All();
        for(int t=0;t<cutter->ta->mesh.elements.m;t++){
            I3 tri=cutter->ta->mesh.elements(t);
            for(int i=0;i<3;++i){
                vertices.Append(cutter->ta->particles.X(tri(i)));
                vertices.Append(cutter->ta->particles.X(tri((i+1)%3)));
            }
        }
        glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
        glColor4f(0.0, 1.0, 0.0, 1.0);
        glDrawArrays(GL_LINES,0,vertices.m);
    }
    
    //material elements
    vertices.Remove_All();
    for(int t=0;t<cutter->ta->mesh.elements.m;t++){
        I3 tri=cutter->ta->mesh.elements(t);
        for(int i=0;i<3;++i)
            vertices.Append(cutter->ta->particles.X(tri(i)));
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
    sim_ta=TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(),.5),10);
    
    //four triangles
//    sim_ta=TRIANGULATED_AREA<T>::Create();
//    sim_ta->particles.Add_Elements(4);
//    sim_ta->particles.X(0)=TV(0,.5);
//    sim_ta->particles.X(1)=TV(-.5,0);
//    sim_ta->particles.X(2)=TV(.5,0);
//    sim_ta->particles.X(3)=TV(0,-.5);
//    sim_ta->mesh.elements.Append(I3(0,1,2));
//    sim_ta->mesh.elements.Append(I3(2,1,3));
//    sim_ta->Update_Number_Nodes();
    
    
    //six triangles
//    sim_ta=TRIANGULATED_AREA<T>::Create();
//    sim_ta->particles.Add_Elements(5);
//    sim_ta->particles.X(0)=TV(0,.5);
//    sim_ta->particles.X(1)=TV(-.5,0);
//    sim_ta->particles.X(2)=TV(0,-.5);
//    sim_ta->particles.X(3)=TV(.5,0);
//    sim_ta->particles.X(4)=TV(0,0);
//    for(int i=0;i<4;++i)
//        sim_ta->mesh.elements.Append(I3(i,(i+1)%4,4));
//    sim_ta->Update_Number_Nodes();
    
    sc=SEGMENTED_CURVE_2D<T>::Create();
    
    sc->particles.Add_Elements(12);
    sc->particles.X(0)=TV(0,0.8);
    sc->particles.X(1)=TV(.0,0);
    sc->particles.X(2)=TV(0,-0.8);
    sc->particles.X(3)=TV(-0.8,0);
    sc->particles.X(4)=TV(.0,0);
    sc->particles.X(5)=TV(0.8,0);
    sc->particles.X(6)=TV(-0.5,-0.5);
    sc->particles.X(7)=TV(0,0);
    sc->particles.X(8)=TV(0.5,0.5);
    sc->particles.X(9)=TV(-0.5,0.5);
    sc->particles.X(10)=TV(0,0);
    sc->particles.X(11)=TV(0.5,-0.5);
    sc->mesh.elements.Append(I2(0,1));
    sc->mesh.elements.Append(I2(1,2));
    sc->mesh.elements.Append(I2(3,4));
    sc->mesh.elements.Append(I2(4,5));
    sc->mesh.elements.Append(I2(6,7));
    sc->mesh.elements.Append(I2(7,8));
    sc->mesh.elements.Append(I2(9,10));
    sc->mesh.elements.Append(I2(10,11));
    sc->Update_Number_Nodes();
    
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
    cout << "initialized mesh\n";
}

//one triangle degeneracy tests
void Initialize_Meshes0()
{
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(4);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(0,0);
    sim_ta->particles.X(2)=TV(.5,0);
    sim_ta->particles.X(3)=TV(0,-0.5);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->mesh.elements.Append(I3(2,1,3));
    sim_ta->Update_Number_Nodes();
    
    sc=SEGMENTED_CURVE_2D<T>::Create();
    
    sc->particles.Add_Elements(12);
    sc->particles.X(0)=TV(0,0.8);
    sc->particles.X(1)=TV(0,0.1);
    sc->particles.X(2)=TV(0,-0.1);
    sc->particles.X(3)=TV(-0.1,0);
    sc->particles.X(4)=TV(0.1,0);
    sc->particles.X(5)=TV(0.8,0);
    sc->particles.X(6)=TV(0.5,0);
    sc->particles.X(7)=TV(0,0.25);
    sc->particles.X(8)=TV(0,0.5);
    sc->particles.X(9)=TV(0.25,0);
    sc->particles.X(10)=TV(-0.1,-0.1);
    sc->particles.X(11)=TV(0.3,0.3);
    sc->mesh.elements.Append(I2(0,1));
    sc->mesh.elements.Append(I2(1,2));
    sc->mesh.elements.Append(I2(3,4));
    sc->mesh.elements.Append(I2(4,5));
    sc->mesh.elements.Append(I2(6,7));
    sc->mesh.elements.Append(I2(8,9));
    sc->mesh.elements.Append(I2(10,11));
    sc->Update_Number_Nodes();
    
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
    cout << "initialized mesh\n";
}

int main(int argc, char **argv)
{
    argc1 = argc;
    argv1 = argv;
    
    Initialize_Meshes();
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
    glutInitWindowSize(window_width,window_height);
    string window_name="cutting";
    glutCreateWindow(window_name.c_str());
    glClearColor(0.0,0.0,0.0,1.0);
    glEnable(GL_DEPTH_TEST);
    
    //glutSpecialFunc(Special_Key);
    glutKeyboardFunc(Key);
    glutMouseFunc(Mouse);
    glutMotionFunc(Motion);
    glutDisplayFunc(Render);
    glutReshapeFunc(Reshape);
    
    glutMainLoop();
    
    return 0;
}
