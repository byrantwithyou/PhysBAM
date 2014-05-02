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
#include "CONSISTENT_INTERSECTIONS.h"
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>

using namespace PhysBAM;
using namespace std;

typedef double T;
typedef VECTOR<T, 2> TV;
typedef VECTOR<int,2> I2;
typedef VECTOR<int,3> I3;

namespace VS{
    void read_tsc(){__asm__("rdtsc");}
    struct timeval starttime,stoptime;
    void start_timer(){gettimeofday(&starttime,NULL);read_tsc();}
    void stop_timer(){gettimeofday(&stoptime,NULL);read_tsc();}
    double get_time(){return ((double)stoptime.tv_sec-(double)starttime.tv_sec)*1000+(double)1e-3*(double)stoptime.tv_usec-(double)1e-3*(double)starttime.tv_usec;}
}

//global variables
int argc1;
char **argv1;
int window_width=600;
int window_height=600;
TV starting_position;
T trans_speed=0.05;
T scale_speed=1.2;
TRIANGULATED_AREA<T>* sim_ta=NULL;
SEGMENTED_CURVE_2D<T>* sc=NULL;
CUTTING<TV>* cutter=NULL;
ARRAY<int> labels;
HASHTABLE<int> dragging_particles;
bool dragging=false,cutting=false,draw_sim=true,draw_cutting_edge=false;

bool recording = false;
T timestamp = 0;
int fi = 0;
string writing_directory = "zoom_in_degeneracy";
void Write_Recording_Info(bool cut = false, bool new_time = true) {
    stringstream ss;
    ss << fi;
    string fis = ss.str();
    while (fis.length() < 6) {
        fis = "0" + fis;
    }
    
    if (new_time) {
        VS::stop_timer();
        timestamp = VS::get_time();
    }
    
    ofstream fs;
    fs.open(writing_directory+"/info-"+string(fis)+string(".txt"));
    fs << "TR " << 1 << " " << 1 << endl;
    fs << "BL " << -1 << " " << -1 << endl;
    fs << "TIME(millisec) " << timestamp << endl;
    fs << "CUT " << cut << endl;
    fs.close();
    
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/parent-"+fis+string(".tri2d.gz"), *sim_ta);
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/child-"+fis+string(".tri2d.gz"), *(cutter->ta));
    
    if (cut) {
        //write degeneracy
        CONSISTENT_INTERSECTIONS<TV> intersections(*(cutter->ta),*sc);
        intersections.Set_Tol();
        intersections.Compute();
        FILE_UTILITIES::Write_To_File<T>(writing_directory+"/hash-"+fis+string(".data"),intersections.hash_vv,intersections.hash_ve,intersections.hash_ev,intersections.hash_ee,intersections.hash_fv);
        
        //reinitialize ta
        TRIANGULATED_AREA<T>* nta=TRIANGULATED_AREA<T>::Create();
        nta->particles.Resize(cutter->ta->particles.X.m);
        for(int i=0;i<cutter->ta->particles.X.m;++i)
            nta->particles.X(i)=cutter->ta->particles.X(i);
        nta->mesh.elements=cutter->ta->mesh.elements;
        nta->Update_Number_Nodes();
        delete cutter->ta;
        cutter->ta=nta;
        //reinitialize sc
        SEGMENTED_CURVE_2D<T>* nsc=SEGMENTED_CURVE_2D<T>::Create();
        nsc->particles.Resize(sc->particles.X.m);
        for(int i=0;i<sc->particles.X.m;++i)
            nsc->particles.X(i)=sc->particles.X(i);
        nsc->mesh.elements=sc->mesh.elements;
        nsc->Update_Number_Nodes();
        delete sc;
        sc=nsc;
        cutter->sc = sc;
        
        FILE_UTILITIES::Write_To_File<T>(writing_directory+"/curve-"+fis+string(".curve2d.gz"), *sc);
    }
    
    ++fi;
}

void Run_Cutter()
{
    if (recording) {
        Write_Recording_Info(true);
    }
    
    if(sc->mesh.elements.m>0){
        cutter->Run(.1);
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
    cutter->ta->mesh.Identify_Connected_Components(labels);
    
    if (recording) {
        Write_Recording_Info(false, false);
    }
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
        if (recording) {
            return;
        }
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
    if (recording) {
        Write_Recording_Info();
    }
    
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

//1.one triangle, cut by mouse to show split.
//2.one triangle, degenerate cases.
//3.two triangles, not sharing node, cut by mouse to show effect without merging.
//4.two triangles, with merging, cut by mouse.
//5.two triangles, cut on the sharing edge, cut on nodes.
//6.two triangles, show sim triangles.
//7.disk with many triangles.

//one triangle split
void Step1()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(3);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(-0.5/sqrt(3),0);
    sim_ta->particles.X(2)=TV(.5/sqrt(3),0);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->Update_Number_Nodes();
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}

//one triangle one segment cut on node
void Step2()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(3);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(-0.5/sqrt(3),0);
    sim_ta->particles.X(2)=TV(.5/sqrt(3),0);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->Update_Number_Nodes();
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    sc->particles.Add_Elements(2);
    sc->particles.X(0)=TV(0,0.5);
    sc->particles.X(1)=TV(0,0);
    sc->mesh.elements.Append(I2(0,1));
    sc->Update_Number_Nodes();
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}

//one triangle 3 segments cut on node
void Step3()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(3);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(-0.5/sqrt(3),0);
    sim_ta->particles.X(2)=TV(.5/sqrt(3),0);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->Update_Number_Nodes();
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    sc->particles.Add_Elements(6);
    sc->particles.X(0)=TV(0,0.5);
    sc->particles.X(1)=TV(0,0);
    sc->particles.X(2)=TV(0.5/sqrt(3),0);
    sc->particles.X(3)=TV(-0.25/sqrt(3),0.25);
    sc->particles.X(4)=TV(-0.5/sqrt(3),0);
    sc->particles.X(5)=TV(0.25/sqrt(3),0.25);
    sc->mesh.elements.Append(I2(0,1));
    sc->mesh.elements.Append(I2(2,3));
    sc->mesh.elements.Append(I2(4,5));
    sc->Update_Number_Nodes();
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}

void Step4()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(6);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(-0.5,0);
    sim_ta->particles.X(2)=TV(0,-0.5);
    sim_ta->particles.X(3)=TV(0,0.5);
    sim_ta->particles.X(4)=TV(0,-0.5);
    sim_ta->particles.X(5)=TV(0.5,0);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->mesh.elements.Append(I3(3,4,5));
    sim_ta->Update_Number_Nodes();
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}

void Step5()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(4);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(-0.5,0);
    sim_ta->particles.X(2)=TV(0,-0.5);
    sim_ta->particles.X(3)=TV(0.5,0);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->mesh.elements.Append(I3(0,2,3));
    sim_ta->Update_Number_Nodes();
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}
void Step6()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(4);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(-0.5,0);
    sim_ta->particles.X(2)=TV(0,-0.5);
    sim_ta->particles.X(3)=TV(0.5,0);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->mesh.elements.Append(I3(0,2,3));
    sim_ta->Update_Number_Nodes();
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    sc->particles.Add_Elements(2);
    sc->particles.X(0)=TV(0,0.5);
    sc->particles.X(1)=TV(0,-0.5);
    sc->mesh.elements.Append(I2(0,1));
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}

void Step7()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TRIANGULATED_AREA<T>::Create();
    sim_ta->particles.Add_Elements(4);
    sim_ta->particles.X(0)=TV(0,.5);
    sim_ta->particles.X(1)=TV(-0.5,0);
    sim_ta->particles.X(2)=TV(0,-0.5);
    sim_ta->particles.X(3)=TV(0.5,0);
    sim_ta->mesh.elements.Append(I3(0,1,2));
    sim_ta->mesh.elements.Append(I3(0,2,3));
    sim_ta->Update_Number_Nodes();
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    sc->particles.Add_Elements(2);
    sc->particles.X(0)=TV(-0.5,0);
    sc->particles.X(1)=TV(0.5,0);
    sc->mesh.elements.Append(I2(0,1));
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}
void Step8()
{
    if(sim_ta) delete sim_ta;
    sim_ta=TESSELLATION::Generate_Triangles(SPHERE<TV>(TV(),.5),50);
    
    if(sc) delete sc;
    sc=SEGMENTED_CURVE_2D<T>::Create();
    sc->particles.Add_Elements(4);
    sc->Update_Number_Nodes();
    sc->particles.X(0)=TV(-1,0);
    sc->particles.X(1)=TV(0,0);
    sc->particles.X(2)=TV(1,0);
    sc->mesh.elements.Append(I2(0,1));
    sc->mesh.elements.Append(I2(1,2));
    RANDOM_NUMBERS<T> rn;
    T scale=1e-7;
    for (int i = 0; i < sim_ta->particles.X.m; ++i) {
        sim_ta->particles.X(i) += TV(rn.Get_Uniform_Number(-1,1)*scale, rn.Get_Uniform_Number(-1,1)*scale);
    }
    
    if(cutter) delete cutter;
    cutter=new CUTTING<TV>(sim_ta,sc);
    Run_Cutter();
}

void Step9()
{
    int t = 0;
    RANDOM_NUMBERS<T> rn;
    //cout << rn.Get_Uniform_Number(-1,1) << endl;
    while (t++ < 1e6) {
        if (0) {
            if(sim_ta) delete sim_ta;
            sim_ta=TRIANGULATED_AREA<T>::Create();
            FILE_UTILITIES::Read_From_File<T>("failDisk.tri2d.gz",*sim_ta);
            
            if(sc) delete sc;
            sc=SEGMENTED_CURVE_2D<T>::Create();
            FILE_UTILITIES::Read_From_File<T>("failCircle.seg2d.gz",*sc);
            
            if(cutter) delete cutter;
            cutter=new CUTTING<TV>(sim_ta,sc);
            Run_Cutter();
            
            if(labels.Max()<2) {
                return;
            }
        }
        else {
            TV shift;
            
            if(sc) delete sc;
            sc=TESSELLATION::Tessellate_Boundary(SPHERE<VECTOR<T,2> >(TV(),0.3),5);
            int np=sc->particles.X.m;
            
            if(sim_ta) delete sim_ta;
            sim_ta=TRIANGULATED_AREA<T>::Create();
            sim_ta->particles.Add_Elements(2*sc->particles.X.m+1);
            sim_ta->Update_Number_Nodes();
            sim_ta->particles.X(0)=TV();
            for (int i=0; i<np; ++i) {
                sim_ta->particles.X(i+1)=sc->particles.X(i);
                sim_ta->particles.X(i+1+np)=sc->particles.X(i)*2;
                sim_ta->mesh.elements.Append(I3(0,i+1,(i+1)%np+1));
                //i+1,(i+2)%np,
                sim_ta->mesh.elements.Append(I3(i+1,i+1+np,(i+1)%np+1));
                sim_ta->mesh.elements.Append(I3((i+1)%np+1,i+1+np,(i+1)%np+1+np));
            }
            
            for (int i=0; i<np; ++i) {
                sc->particles.X(i)+=shift;
            }
            for (int i=0; i<sim_ta->particles.X.m; ++i) {
                sim_ta->particles.X(i)+=shift;
            }
            //perturb
            T scale = 1e-6;
            cout.precision(15);
            cout << "ta particles: " ;
            for (int i = 0; i < sim_ta->particles.X.m; ++i) {
                //cout << sim_ta->particles.X(i) << " ";
                //cout << r*scale << endl;
                sim_ta->particles.X(i) += TV(rn.Get_Uniform_Number(-1,1)*scale, rn.Get_Uniform_Number(-1,1)*scale);
                //cout << sim_ta->particles.X(i) << " ";
            }
            cout << endl;
            
            TRIANGULATED_AREA<T> *tac = TRIANGULATED_AREA<T>::Create();
            tac->particles.Add_Elements(sim_ta->particles.X.m);
            tac->Update_Number_Nodes();
            for (int i = 0; i < sim_ta->particles.X.m; ++i) {
                tac->particles.X(i) = sim_ta->particles.X(i);
            }
            tac->mesh.elements = sim_ta->mesh.elements;
            
            cout << sc->particles.X.m << " sc particles: " ;
            for (int i = 0; i < sc->particles.X.m; ++i) {
                sc->particles.X(i) += TV(rn.Get_Uniform_Number(-1,1)*scale, rn.Get_Uniform_Number(-1,1)*scale);
                //cout << sc->particles.X(i) << " ";
            }
            if(cutter) delete cutter;
            cutter=new CUTTING<TV>(sim_ta,sc);
            Run_Cutter();
            cout << t << " " << labels.Max() << endl;
            if(labels.Max()==2) {
                FILE_UTILITIES::Write_To_File<T>("failDisk.tri2d.gz",*tac);
                FILE_UTILITIES::Write_To_File<T>("failCircle.seg2d.gz",*sc);
                return;
            }
        }
    }
}

void Step10()
{
    int t = 0;
    RANDOM_NUMBERS<T> rn;
    while (t++ < 1e6) {
        if (0) {
            if(sim_ta) delete sim_ta;
            sim_ta=TRIANGULATED_AREA<T>::Create();
            FILE_UTILITIES::Read_From_File<T>("failTri.tri2d.gz",*sim_ta);
            
            if(sc) delete sc;
            sc=SEGMENTED_CURVE_2D<T>::Create();
            FILE_UTILITIES::Read_From_File<T>("failTri.seg2d.gz",*sc);
            
            if(cutter) delete cutter;
            cutter=new CUTTING<TV>(sim_ta,sc);
            Run_Cutter();
            
            if(labels.Max()<2) {
                return;
            }
        }
        else {
            if(sim_ta) delete sim_ta;
            sim_ta=TRIANGULATED_AREA<T>::Create();
            sim_ta->particles.Add_Elements(4);
            sim_ta->Update_Number_Nodes();
            sim_ta->particles.X(0) = TV(0,0.5);
            sim_ta->particles.X(1) = TV(-0.5,-0.5);
            sim_ta->particles.X(2) = TV(0.5,-0.5);
            sim_ta->particles.X(3) = TV(0,0);
            sim_ta->mesh.elements.Append(I3(0,1,3));
            sim_ta->mesh.elements.Append(I3(1,2,3));
            sim_ta->mesh.elements.Append(I3(2,0,3));
            
            if(sc) delete sc;
            sc=SEGMENTED_CURVE_2D<T>::Create();
            sc->particles.Add_Elements(3);
            sc->Update_Number_Nodes();
            sc->particles.X(0) = TV(-1,0);
            sc->particles.X(1) = TV(0,0);
            sc->particles.X(2) = TV(1,0);
            sc->mesh.elements.Append(I2(0,1));
            sc->mesh.elements.Append(I2(1,2));
            
            //perturb
            T scale = 1e-6;
            cout.precision(15);
            for (int i = 0; i < sim_ta->particles.X.m; ++i) {
                //cout << sim_ta->particles.X(i) << endl;
                //cout << r*scale << endl;
                T r = rn.Get_Uniform_Number(-1,1);
                sim_ta->particles.X(i) += TV(r*scale, r*scale);
                cout << sim_ta->particles.X(i) << " ";
            }
            cout << endl;
            TRIANGULATED_AREA<T> *tac = TRIANGULATED_AREA<T>::Create();
            tac->particles.Add_Elements(sim_ta->particles.X.m);
            tac->Update_Number_Nodes();
            for (int i = 0; i < sim_ta->particles.X.m; ++i) {
                tac->particles.X(i) = sim_ta->particles.X(i);
            }
            tac->mesh.elements = sim_ta->mesh.elements;
            cout << "sc particles: " ;
            for (int i = 0; i < sc->particles.X.m; ++i) {
                T r = rn.Get_Uniform_Number(-1,1);
                sc->particles.X(i) += TV(r*scale, r*scale);
                cout << sc->particles.X(i) << " ";
            }
            cout << endl;
            if(cutter) delete cutter;
            cutter=new CUTTING<TV>(sim_ta,sc);
            Run_Cutter();
            cout << t << " " << labels.Max() << endl;
            if(labels.Max()!=2) {
                FILE_UTILITIES::Write_To_File<T>("failTri.tri2d.gz",*tac);
                FILE_UTILITIES::Write_To_File<T>("failTri.seg2d.gz",*sc);
                return;
            }
        }
    }
}

void Step11()
{
    int t = 0;
    RANDOM_NUMBERS<T> rn;
    while (t++ < 1e6) {
        if (0) {
            if(sim_ta) delete sim_ta;
            sim_ta=TRIANGULATED_AREA<T>::Create();
            FILE_UTILITIES::Read_From_File<T>("failTri.tri2d.gz",*sim_ta);
            
            if(sc) delete sc;
            sc=SEGMENTED_CURVE_2D<T>::Create();
            FILE_UTILITIES::Read_From_File<T>("failTri.seg2d.gz",*sc);
            
            if(cutter) delete cutter;
            cutter=new CUTTING<TV>(sim_ta,sc);
            Run_Cutter();
            
            if(labels.Max()<2) {
                return;
            }
        }
        else {
            TV shift(5,5);
            if(sim_ta) delete sim_ta;
            sim_ta=TRIANGULATED_AREA<T>::Create();
            sim_ta->particles.Add_Elements(3);
            sim_ta->Update_Number_Nodes();
            sim_ta->particles.X(0) = TV(0.5,0.5)+shift;
            sim_ta->particles.X(1) = TV(0,0)+shift;
            sim_ta->particles.X(2) = TV(.5,-0.5)+shift;
            sim_ta->mesh.elements.Append(I3(0,1,2));
            
            if(sc) delete sc;
            sc=SEGMENTED_CURVE_2D<T>::Create();
            sc->particles.Add_Elements(3);
            sc->Update_Number_Nodes();
            sc->particles.X(0) = TV(-1,0)+shift;
            sc->particles.X(1) = TV(0,0)+shift;
            sc->particles.X(2) = TV(1,0)+shift;
            sc->mesh.elements.Append(I2(0,1));
            sc->mesh.elements.Append(I2(1,2));
            
            //perturb
            T scale = 1e-6;
            cout.precision(15);
            cout << "ta particles: " ;
            for (int i = 0; i < sim_ta->particles.X.m; ++i) {
                //cout << sim_ta->particles.X(i) << endl;
                //cout << r*scale << endl;
                sim_ta->particles.X(i) += TV(rn.Get_Uniform_Number(-1,1)*scale, rn.Get_Uniform_Number(-1,1)*scale);
                cout << sim_ta->particles.X(i) << " ";
            }
            cout << endl;
            TRIANGULATED_AREA<T> *tac = TRIANGULATED_AREA<T>::Create();
            tac->particles.Add_Elements(sim_ta->particles.X.m);
            tac->Update_Number_Nodes();
            for (int i = 0; i < sim_ta->particles.X.m; ++i) {
                tac->particles.X(i) = sim_ta->particles.X(i);
            }
            tac->mesh.elements = sim_ta->mesh.elements;
            cout << "sc particles: " ;
            for (int i = 0; i < sc->particles.X.m; ++i) {
                sc->particles.X(i) += TV(rn.Get_Uniform_Number(-1,1)*scale, rn.Get_Uniform_Number(-1,1)*scale);
                cout << sc->particles.X(i) << " ";
            }
            cout << endl;
            if(cutter) delete cutter;
            cutter=new CUTTING<TV>(sim_ta,sc);
            Run_Cutter();
            cout << t << " " << labels.Max() << endl;
//            if(labels.Max()<2) {
//                FILE_UTILITIES::Write_To_File<T>("failTri.tri2d.gz",*tac);
//                FILE_UTILITIES::Write_To_File<T>("failTri.seg2d.gz",*sc);
//                return;
//            }
        }
    }
}

void Step12()
{
    int t = 0;
    RANDOM_NUMBERS<T> rn;
    //deleteing nodes
    if (0) {
        if(sc) delete sc;
        sc=SEGMENTED_CURVE_2D<T>::Create();
        FILE_UTILITIES::Read_From_File<T>("fail.seg2d.gz",*sc);
        int np = sc->particles.X.m;
        HASHTABLE<int> s;
        for (int i = 0; i < np; ++i) {
            s.Insert(i);
        }
        for (int i = 1; i < np-1; ++i) {
            if(sim_ta) delete sim_ta;
            sim_ta=TRIANGULATED_AREA<T>::Create();
            FILE_UTILITIES::Read_From_File<T>("fail.tri2d.gz",*sim_ta);

            SEGMENTED_CURVE_2D<T>* sc2=SEGMENTED_CURVE_2D<T>::Create();
            for (int j = 0; j < np; ++j) {
                if (j != i && s.Contains(j)) {
                    sc2->particles.Add_Elements(1);
                    sc2->particles.X(sc2->particles.X.m-1)=sc->particles.X(j);
                }
            }
            sc2->Update_Number_Nodes();
            
            for (int j = 0; j < sc2->particles.X.m-1; ++j) {
                sc2->mesh.elements.Append(I2(j,j+1));
            }
            
            if(cutter) delete cutter;
            cutter=new CUTTING<TV>(sim_ta,sc2);
            Run_Cutter();
            cout << "*************************************CCs: " << labels.Max() << endl;
            if(labels.Max()==1) {
                s.Delete(i);
            }
        }
        cout << s << endl;
        return;
    }
    if (1) {
        if(sc) delete sc;
        sc=SEGMENTED_CURVE_2D<T>::Create();
        FILE_UTILITIES::Read_From_File<T>("fail.seg2d.gz",*sc);
        if(sim_ta) delete sim_ta;
        sim_ta=TRIANGULATED_AREA<T>::Create();
        FILE_UTILITIES::Read_From_File<T>("fail.tri2d.gz",*sim_ta);

        SEGMENTED_CURVE_2D<T>* sc2=SEGMENTED_CURVE_2D<T>::Create();
        sc2->particles.Add_Elements(4);
        sc2->particles.X(0)=sc->particles.X(0);
        sc2->particles.X(1)=sc->particles.X(4289);
        sc2->particles.X(2)=sc->particles.X(4290);
        sc2->particles.X(3)=sc->particles.X(8579);
        sc2->Update_Number_Nodes();
        for (int j = 0; j < sc2->particles.X.m-1; ++j) {
            sc2->mesh.elements.Append(I2(j,j+1));
        }
        
        if(cutter) delete cutter;
        cutter=new CUTTING<TV>(sim_ta,sc2);
        Run_Cutter();
        return;
    }
    while (t++ < 1e9) {
        if(sim_ta) delete sim_ta;
        sim_ta=TRIANGULATED_AREA<T>::Create();
        sim_ta->particles.Add_Elements(9);
        sim_ta->Update_Number_Nodes();
        sim_ta->particles.X(0)=TV();
        sim_ta->particles.X(1)=TV(0.25,0);
        sim_ta->particles.X(2)=TV(0.5,0);
        sim_ta->particles.X(3)=TV(0.125,0.25);
        sim_ta->particles.X(4)=TV(0.375,0.25);
        sim_ta->particles.X(5)=TV(0.625,0.25);
        sim_ta->particles.X(6)=TV(0.25,0.5);
        sim_ta->particles.X(7)=TV(0.5,0.5);
        sim_ta->particles.X(8)=TV(0.75,0.5);
        sim_ta->mesh.elements.Append(I3(0,1,3));
        sim_ta->mesh.elements.Append(I3(3,1,4));
        sim_ta->mesh.elements.Append(I3(4,1,2));
        sim_ta->mesh.elements.Append(I3(4,2,5));
        sim_ta->mesh.elements.Append(I3(6,3,4));
        sim_ta->mesh.elements.Append(I3(6,4,7));
        sim_ta->mesh.elements.Append(I3(7,4,5));
        sim_ta->mesh.elements.Append(I3(7,5,8));
        
        //perturb mesh
        cout << "ta particles: " ;
        for (int i = 0; i < sim_ta->particles.X.m; ++i) {
            sim_ta->particles.X(i) += TV(pow(10, rn.Get_Uniform_Number(-18,-12)), pow(10,rn.Get_Uniform_Number(-18,-12)));
        }
        cout << endl;
        
        //make a copy
        TRIANGULATED_AREA<T> *tac = TRIANGULATED_AREA<T>::Create();
        tac->particles.Add_Elements(sim_ta->particles.X.m);
        tac->Update_Number_Nodes();
        for (int i = 0; i < sim_ta->particles.X.m; ++i) {
            tac->particles.X(i) = sim_ta->particles.X(i);
        }
        tac->mesh.elements = sim_ta->mesh.elements;
        
        //curve
        if(sc) delete sc;
        sc=SEGMENTED_CURVE_2D<T>::Create();
        T x1=0.125,x2=0.625,y1=0.75,y2=-0.25;
        int n=rn.Get_Uniform_Integer(1,1e4);
        T dx=(x2-x1)/n, dy=(y2-y1)/n;
        sc->particles.Add_Elements(n+1);
        sc->Update_Number_Nodes();
        sc->particles.X(0)=TV(x1,y1);
        sc->particles.X(n)=TV(x2,y2);
        for (int i = 1; i < n; ++i) {
            sc->particles.X(i)=TV(x1+i*dx,y1+i*dy);
        }
        for (int i = 0; i < n; ++i) {
            sc->mesh.elements.Append(I2(i,i+1));
        }
        
        cout << sc->particles.X.m << " sc particles: " ;
        //shift curve
        T shift=pow(10,rn.Get_Uniform_Number(-20,-12));
        cout << "shift: " << shift << endl;
        for (int i = 0; i < sc->particles.X.m; ++i) {
            sc->particles.X(i) += TV(shift,0);
        }
        //jitter curve
        for (int i = 0; i < sc->particles.X.m; ++i) {
            sc->particles.X(i) += TV(pow(10, rn.Get_Uniform_Number(-18,-12)), pow(10,rn.Get_Uniform_Number(-18,-12)));
        }
        cout << endl;
        
        if(cutter) delete cutter;
        cutter=new CUTTING<TV>(sim_ta,sc);
        Run_Cutter();
        cout << t << " " << labels.Max() << endl;
        if(labels.Max()<2) {
            FILE_UTILITIES::Write_To_File<T>("fail.tri2d.gz",*tac);
            FILE_UTILITIES::Write_To_File<T>("fail.seg2d.gz",*sc);
            return;
        }
    }
}

static void Key(unsigned char key, int x, int y)
{
    switch( key ) {
        case 033: // Escape Key
            exit(EXIT_SUCCESS);
            break;
        case '1':
            Step1();
            break;
        case '2':
            Step2();
            break;
        case '3':
            Step3();
            break;
        case '4':
            Step4();
            break;
        case '5':
            Step5();
            break;
        case '6':
            Step6();
            break;
        case '7':
            Step7();
            break;
        case '8':
            Step8();
            break;
        case '9':
            Step9();
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


int main(int argc, char **argv)
{
    argc1 = argc;
    argv1 = argv;
    
    if (!recording) {
        Step12();
    }
    else {
        Step8();
        VS::start_timer();
    }
    
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
