//
//  mesh_cutting_test.cpp
//
//
//  Created by Yuting Wang on 5/24/12.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.
//

#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/ROTATION.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include "CONSISTENT_INTERSECTIONS.h"

#include <fstream>
#include <sstream>

#include "DEFORMABLE_OBJECTS.h"
#include "mesh_cutting_subd.h"

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#include <X11/Xlib.h>
#ifdef __APPLE__  // include Mac OS X verions of headers
//#include <GL/glew.h>
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>

#else // non-Mac OS X operating systems
#include <GL/freeglut.h>
#include <GL/freeglut_ext.h>
#endif  // __APPLE__
namespace VS{
void read_tsc(){__asm__("rdtsc");}
struct timeval starttime,stoptime;
void start_timer(){gettimeofday(&starttime,NULL);read_tsc();}
void stop_timer(){gettimeofday(&stoptime,NULL);read_tsc();}
double get_time(){return ((double)stoptime.tv_sec-(double)starttime.tv_sec)*1000+(double)1e-3*(double)stoptime.tv_usec-(double)1e-3*(double)starttime.tv_usec;}
}

using namespace PhysBAM;
using namespace std;


typedef double T;
typedef PhysBAM::VECTOR<int, 2> I2;
typedef PhysBAM::VECTOR<int, 3> I3;
typedef PhysBAM::VECTOR<int, 4> I4;

typedef PhysBAM::VECTOR<T, 3> TV;
typedef PhysBAM::MATRIX<T, 3> TM;

//global variables
int argc1;
char **argv1;

//mouse motion modes
bool cutting = 1;
bool dragging = 0;
bool rotating = 0;
bool translating = 0;

int dragging_id = -1;
T K = 50000;

PhysBAM::VECTOR<T,2> starting_rotation;
PhysBAM::VECTOR<T,2> starting_dragging;


bool cubes = 0;
int picked_cc_id = -1;
TV starting_position, end_position;
T timestep = 30;
T rotate_speed = 10/57.;
T intrude = 2;

int cutting_surface_id = 0;
fstream interaction_file;
bool writing_interaction_to_file = 0;
int current_frame = 0;
bool run_sim = 0;
int cutting_ratio = 10;

ARRAY<TV> cutting_curve;
MESH_CUTTING<T> *mcut = NULL;
TRIANGULATED_SURFACE<T> *cutting_tri_mesh = NULL;
TETRAHEDRALIZED_VOLUME<T> *sim_volume = NULL;

ARRAY<int> labels;
HASHTABLE<int> picked_nodes;

int window_height = 600;
int window_width = 600;

#define DE cout<<"file "<<__FILE__<<"   line "<<__LINE__<<"  "<<&mcut->volume->particles.X<<"   "<<mcut->volume->particles.X<<endl;

bool recording = false;
T timestamp = 0;
int fi = 0;
string writing_directory = "zoom_in_3d_new";
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
    fs << "TR " << 1 << " " << 1 << " " << 1 << endl;
    fs << "BL " << -1 << " " << -1 << " " << -1 << endl;
    fs << "TIME(millisec) " << timestamp << endl;
    fs << "CUT " << cut << endl;
    fs.close();
    
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/parent-"+fis+string(".tet.gz"), *sim_volume);
    FILE_UTILITIES::Write_To_File<T>(writing_directory+"/child-"+fis+string(".tet.gz"), *(mcut->volume));
    
    if (cut) {
        //write degeneracy
        
        //have to clear boundary and hierarchy information in the volume before computing intersections
        TETRAHEDRALIZED_VOLUME<T> *tv = TETRAHEDRALIZED_VOLUME<T>::Create();
        tv->mesh.elements = mcut->volume->mesh.elements;
        tv->particles.Add_Elements(mcut->volume->particles.X.m);
        tv->Update_Number_Nodes();
        for (int i = 0; i < mcut->volume->particles.X.m; ++i) {
            tv->particles.X(i) = mcut->volume->particles.X(i);
        }
        
        CONSISTENT_INTERSECTIONS<TV> intersections(*tv,*cutting_tri_mesh);
        intersections.Compute();
        FILE_UTILITIES::Write_To_File<T>(writing_directory+"/hash1-"+fis+string(".data"),intersections.hash_vv,intersections.hash_ve,intersections.hash_ev,intersections.hash_ee);
        FILE_UTILITIES::Write_To_File<T>(writing_directory+"/hash2-"+fis+string(".data"),intersections.hash_fv,intersections.hash_vf,intersections.hash_fe,intersections.hash_ef,intersections.hash_tv);

        //reinitialize tv
        TETRAHEDRALIZED_VOLUME<T>* ntv=TETRAHEDRALIZED_VOLUME<T>::Create();
        ntv->particles.Resize(mcut->volume->particles.X.m);
        for(int i=0;i<mcut->volume->particles.X.m;++i)
            ntv->particles.X(i)=mcut->volume->particles.X(i);
        ntv->mesh.elements=mcut->volume->mesh.elements;
        ntv->Update_Number_Nodes();
        delete mcut->volume;
        mcut->volume=ntv;
        
        //reinitialize ts
        TRIANGULATED_SURFACE<T>* nts=TRIANGULATED_SURFACE<T>::Create();
        nts->particles.Resize(cutting_tri_mesh->particles.X.m);
        for(int i=0;i<cutting_tri_mesh->particles.X.m;++i)
            nts->particles.X(i)=cutting_tri_mesh->particles.X(i);
        nts->mesh.elements=cutting_tri_mesh->mesh.elements;
        nts->Update_Number_Nodes();
        delete cutting_tri_mesh;
        cutting_tri_mesh = nts;
        
        FILE_UTILITIES::Write_To_File<T>(writing_directory+"/surface-"+fis+string(".tri3d.gz"), *cutting_tri_mesh);
    }
    
    ++fi;
}

void copy(TETRAHEDRALIZED_VOLUME<T>*& vt, TETRAHEDRALIZED_VOLUME<T> *vf)
{
    if (vt) {
        delete vt;
    }
    vt = TETRAHEDRALIZED_VOLUME<T>::Create();
    vt->particles.Resize(vf->particles.X.m);
    vt->Update_Number_Nodes();
    for (int i = 0; i < vf->particles.X.m; ++i) {
        vt->particles.X(i) = vf->particles.X(i);
    }
    vt->mesh.elements = vf->mesh.elements;
}

void Reshape(GLint newWidth,GLint newHeight) {
    glViewport(0,0,newWidth,newHeight);
    window_width=newWidth;
    window_height=newHeight;
}

bool draw_sim = 1, draw_material_edges = 1, drawing_cutting = 1;
void Render(){
    if (recording) {
        Write_Recording_Info();
    }
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnableClientState(GL_VERTEX_ARRAY);
    ARRAY<TV> vertices;
    
    //cutting mesh
    if (drawing_cutting && cutting_tri_mesh->mesh.elements.m) {
        vertices.Remove_All();
        for(int t=0;t<cutting_tri_mesh->mesh.elements.m;t++){
            I3 e=cutting_tri_mesh->mesh.elements(t);
            for(int i=0;i<3;++i) {
                vertices.Append(cutting_tri_mesh->particles.X(e(i)));
                vertices.Append(cutting_tri_mesh->particles.X(e((i+1)%3)));
            }
        }
        glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
        glColor4f(1, 0, 1, 1.0);
        glDrawArrays(GL_LINES,0,vertices.m);
        
        vertices.Remove_All();
        for(int t=0;t<cutting_tri_mesh->mesh.elements.m;t++){
            I3 e=cutting_tri_mesh->mesh.elements(t);
            for(int i=0;i<3;++i)
                vertices.Append(cutting_tri_mesh->particles.X(e(i)));
        }
        glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
        glColor4f(0.0, 1.0, 1.0, 1.0);
        glDrawArrays(GL_TRIANGLES,0,vertices.m);
    }
    
    //edges of sim mesh
    if(draw_sim){
        vertices.Remove_All();
        ARRAY<I2> segments = sim_volume->mesh.boundary_mesh->segment_mesh->elements;
        for (int i = 0; i < segments.m; ++i) {
            vertices.Append(sim_volume->particles.X(segments(i)(0)));
            vertices.Append(sim_volume->particles.X(segments(i)(1)));
        }
        glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
        glColor4f(0.0, 0.0, 1.0, 1.0);
        glDrawArrays(GL_LINES,0,vertices.m);
    }
    
    //edges material mesh
    if(draw_material_edges){
        vertices.Remove_All();
        ARRAY<I2> segments = mcut->volume->mesh.boundary_mesh->segment_mesh->elements;
        for (int i = 0; i < segments.m; ++i) {
            vertices.Append(mcut->volume->particles.X(segments(i)(0)));
            vertices.Append(mcut->volume->particles.X(segments(i)(1)));
        }
        glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
        glColor4f(0, 1, 0, 1.0);
        glDrawArrays(GL_LINES,0,vertices.m);
    }
    
    //material elements
    vertices.Remove_All();
    ARRAY<I3> boundary_tri = mcut->volume->mesh.boundary_mesh->elements;
    for(int t=0;t<boundary_tri.m;t++){
        I3 tri=boundary_tri(t);
        for(int i=0;i<3;++i)
            vertices.Append(mcut->volume->particles.X(tri(i)));
    }
    glVertexPointer(TV::m, GL_DOUBLE,0,vertices.base_pointer);
    glColor4f(1.0, 0.0, 0.0, 0.0);
    glDrawArrays(GL_TRIANGLES,0,vertices.m);
    
    glutSwapBuffers();
    glDisableClientState(GL_VERTEX_ARRAY);
}


static void SpecialKey( int key, int x, int y )
{
}

void Initialize(bool);

void Translate_Volume(const TV& translation)
{
    for (int i = 0; i < sim_volume->particles.X.m; i++) {
        sim_volume->particles.X(i) += translation;
    }
    mcut->Update_Cutting_Particles();
    
    for (int i = 0; i < mcut->dragging_targets.m; i++) {
        mcut->dragging_targets(i) += translation;
    }
    
    for (int i = 0; i < cutting_tri_mesh->particles.X.m; i++) {
        cutting_tri_mesh->particles.X(i) += translation;
    }
}

void Scale_Volume(T scale) {
    for (int i = 0; i < sim_volume->particles.X.m; i++) {
        sim_volume->particles.X(i) *= scale;
    }
    mcut->Update_Cutting_Particles();
    
    for (int i = 0; i < cutting_tri_mesh->particles.X.m; i++) {
        cutting_tri_mesh->particles.X(i) *= scale;
    }
}

static void Key( unsigned char key, int x, int y )
{
    TM r;
    switch( key ) {
        case 033: // Escape Key
            exit( EXIT_SUCCESS );
            break;
        case 'q': 
            Initialize(1);
            break;
        case 'm':
            cutting = 0;
            translating = 0;
            dragging = 1;
            rotating = 0;
            break;
        case 't':
            cutting = 0;
            translating = 1;
            dragging = 0;
            rotating = 0;
            break;
        case 'c': 
            cutting = 1;
            translating = 0;
            dragging = 0;
            rotating = 0;
            break;
        case 'r': 
            cutting = 0;
            translating = 0;
            dragging = 0;
            rotating = 1;
            break;
        case 'w': 
            Translate_Volume(TV(0,0.05,0));
            break;
        case 's': 
            Translate_Volume(TV(0,-0.05,0));
            break;
        case 'a': 
            Translate_Volume(TV(-0.05,0,0));
            break;
        case 'd': 
            Translate_Volume(TV(0.05,0,0));
            break;
        case 'f':
            Scale_Volume(1.1);
            break;
        case 'g':
            Scale_Volume(0.9);
            break;
        case 'e': 
            drawing_cutting = !drawing_cutting;
            break;
        case 'z':
            draw_sim = !draw_sim;
            break;
        case 'x':
            draw_material_edges = !draw_material_edges;
            break;
    }
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    T xcoor = 2 * x / T(window_width) - 1;
    T ycoor = 1 - 2 * y / T(window_height);
    if (state == GLUT_DOWN) {
        if (button == GLUT_LEFT_BUTTON){    
            if (dragging) {
                //cout << PhysBAM::VECTOR<T,2>(xcoor,ycoor) << endl;
                if (writing_interaction_to_file) {
                    interaction_file << "D " << xcoor << " " << ycoor << endl;
                }
                starting_dragging[0] = xcoor;
                starting_dragging[1] = ycoor;
                dragging_id = mcut->Compute_Intersection(xcoor, ycoor);
                cout << "dragging id: " << dragging_id << endl;
            }
            else if (rotating) {
                starting_rotation[0] = xcoor;
                starting_rotation[1] = ycoor;
                if (writing_interaction_to_file) {
                    interaction_file << "R " << xcoor << " " << ycoor << endl;
                }
            }
            else if (cutting){
                delete cutting_tri_mesh;
                cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
                cutting_curve.Remove_All();
                starting_position = TV( xcoor, ycoor, -intrude/2);
                cutting_curve.Append(starting_position);
            }
            else {
                starting_position = TV( xcoor, ycoor, -intrude/2);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                glBegin(GL_TRIANGLES);
                ARRAY<I4>& tets = mcut->volume->mesh.elements;
                for(int t=0;t<tets.m;t++){
                    I4 tet=tets(t);
                    glColor4d(labels(t)/255.,0,0,0);
                    TV p = mcut->volume->particles.X(tet(0));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(1));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(2));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(0));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(1));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(3));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(1));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(2));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(3));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(2));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(0));
                    glVertex3d(p(0), p(1), p(2));
                    p = mcut->volume->particles.X(tet(3));
                    glVertex3d(p(0), p(1), p(2));
                }
                glEnd();
                unsigned char val;
                glReadPixels(x, window_height-y, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &val);
                cout << "picked cc: " << (int)val << endl;
                picked_nodes.Clean_Memory();
                for(int t=0;t<tets.m;t++)
                    if(labels(t)==(int)val)
                        for(int i=0;i<4;++i){
                            int f=mcut->weights_in_sim(tets(t)(i)).id;
                            for(int j=0;j<4;++j)
                                picked_nodes.Set(sim_volume->mesh.elements(f)(j));
                        }
            }
        }
    }
    else if (state == GLUT_UP) {
        if (button == GLUT_LEFT_BUTTON) { 
            if (translating) {
                picked_cc_id = -1;
            }
            else if (dragging) {
                dragging_id = -1;
            }
            else if (cutting) {     
                int n = cutting_curve.m;
                if (n>1){
                    int n_in = 20;
                    T d_in = intrude / n_in;
                    int pid1 = cutting_tri_mesh->particles.X.m;
                    cutting_tri_mesh->particles.Add_Elements((n_in + 1) * n);
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n_in+1; ++j) {
                            cutting_tri_mesh->particles.X(pid1) = cutting_curve(i) + TV(0,0,j * d_in);
                            pid1++;
                        }
                    }
                    cutting_tri_mesh->Update_Number_Nodes();
                    
                    n--;
                    int np = n_in + 1;
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n_in; ++j) {
                            int p = i * np + j;
                            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(p,p+1,p+np));
                            cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(p+1,p+np,p+np+1));
                        }
                    }
                    cutting_tri_mesh->mesh.Initialize_Boundary_Mesh();
                    
                    //writing cutting information to file for high res simulation
                    if (writing_interaction_to_file){
                        stringstream ss; ss<<cutting_surface_id;
                        interaction_file<< "c " << cutting_surface_id << endl;
                        FILE_UTILITIES::Write_To_File<T>(string("cutting_surfaces/")+ss.str()+string(".tet.gz"), cutting_tri_mesh->particles, cutting_tri_mesh->mesh);
                        cutting_surface_id++;
                    }
                    if (recording) {
                        Write_Recording_Info(true);
                    }
                    mcut->Cut(*cutting_tri_mesh, true, false, 0);
                    mcut->Connected_Components(mcut->volume, labels);
                    cout << labels.m << " labels max: " << labels.Max() << ", " << mcut->volume->mesh.elements.m << endl;
                    if (recording) {
                        Write_Recording_Info(false);
                    }
                }
                glutPostRedisplay();
            } 
        }                  
    }
}

void rotate_meshes(T xx, T yy)
{
    T dx = xx - starting_rotation[0];
    T dy = yy - starting_rotation[1];

    TM r = TM::Rotation_Matrix(TV(-dy,dx,0),-sqrt(dx*dx+dy*dy));
    for (int i = 0; i < sim_volume->particles.X.m; i++) {
        sim_volume->particles.X(i) = r * sim_volume->particles.X(i);
        sim_volume->particles.V(i) = r * sim_volume->particles.V(i);
    }
    mcut->Update_Cutting_Particles();

    for (int i = 0; i < mcut->dragging_targets.m; i++) {
        mcut->dragging_targets(i) = r * mcut->dragging_targets(i);
    }

    for (int i = 0; i < cutting_tri_mesh->particles.X.m; i++) {
        cutting_tri_mesh->particles.X(i) = r * cutting_tri_mesh->particles.X(i);
    }
    
    starting_rotation[0] = xx;
    starting_rotation[1] = yy;
}

void motion(int x, int y)
{
    T xx = 2 * x / T(window_width) - 1;
    T yy = 1 - 2 * y / T(window_height);
    end_position = TV( xx, yy, -intrude/2);

    if (cutting) {
        if ((end_position - cutting_curve(cutting_curve.m-1)).Magnitude() > 1e-3)
            cutting_curve.Append(end_position);
        return;
    }
    else if(dragging) {
        if(dragging_id >= 0){
            T shiftx = xx-starting_dragging[0];
            T shifty = yy-starting_dragging[1];
            starting_dragging[0] = xx;
            starting_dragging[1] = yy;
            if (run_sim) {
                mcut->dragging_targets(dragging_id)(0) += shiftx;
                mcut->dragging_targets(dragging_id)(1) += shifty;
                if (writing_interaction_to_file) {
                    interaction_file << "d " << dragging_id << " " << xx << " " << yy << endl;
                }
            }
        }
    }
    else if(rotating) {
        if (writing_interaction_to_file) {
            interaction_file << "r " << xx << " " << yy << endl;
        }
        rotate_meshes(xx, yy);
    }
    else {
        T shiftx = xx - starting_position(0);
        T shifty = yy - starting_position(1);
        starting_position(0) = xx;
        starting_position(1) = yy;
        for(typename HASHTABLE<int>::ITERATOR it(picked_nodes);it.Valid();it.Next()) {
            mcut->sim_volume->particles.X(it.Key())(0) += shiftx;
            mcut->sim_volume->particles.X(it.Key())(1) += shifty;
        }
        mcut->Update_Cutting_Particles();
    }
    glutPostRedisplay();
}

template<class T>
void perturb_particles(TETRAHEDRALIZED_VOLUME<T>& v)
{
    RANDOM_NUMBERS<T> r;
    for(int i = 0; i < v.particles.X.m; i++){
        v.particles.X(i) += r.Get_Number() * 0.01;
    }
}

void initialize_cutting_mesh()
{
    //can't handle this case...
//    cutting_tri_mesh->particles.Add_Elements(4);
//    cutting_tri_mesh->particles.X(0) = TV(-0.2, -0.2, -1);    
//    cutting_tri_mesh->particles.X(1) = TV(0.5, 0.5, -1);    
//    cutting_tri_mesh->particles.X(2) = TV(0.5, 0.5, 1);    
//    cutting_tri_mesh->particles.X(3) = TV(-0.2, -0.2, 1);    
//    cutting_tri_mesh->mesh.elements.Append(VECTOR<int,3>(0,1,3));
//    cutting_tri_mesh->mesh.elements.Append(VECTOR<int,3>(1,2,3));

    cutting_tri_mesh->particles.Add_Elements(4);
    cutting_tri_mesh->particles.X(0) = TV(-0.1, -0.2, -1);    
    cutting_tri_mesh->particles.X(1) = TV(0.5, 0.5, -1);    
    cutting_tri_mesh->particles.X(2) = TV(0.5, 0.5, 1);    
    cutting_tri_mesh->particles.X(3) = TV(-0.1, -0.2, 1);    
    cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(0,1,3));
    cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(1,2,3));
}

void initialize_volume1()
{       
    sim_volume->particles.Add_Elements(4);
    sim_volume->Update_Number_Nodes();
    sim_volume->particles.X(0) = TV(0.5, 0.0, 0.0);
    sim_volume->particles.X(1) = TV(0.0, 0, -0.5);
    sim_volume->particles.X(2) = TV(-0.5, 0.0, 0.0);
    sim_volume->particles.X(3) = TV(0.0, 0.5, 0);
    sim_volume->Update_Number_Nodes();
    
    sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0,1,2,3));
}

void initialize_volume2()
{
    sim_volume->particles.Add_Elements(5);
    sim_volume->Update_Number_Nodes();
    sim_volume->particles.X(0) = TV(-0.5, 0.0, -0.5);    
    sim_volume->particles.X(1) = TV(0.5, 0.0, -0.5);  
    sim_volume->particles.X(2) = TV(0.5, 0.0, 0.5);  
    sim_volume->particles.X(3) = TV(0.0, 0.0, 0.0);    
    sim_volume->particles.X(4) = TV(0.0, 0.5, 0.0);    
    for (int i = 0; i < 2; i++) {
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(i+1,i,3,4));
    }
}

void initialize_volume3()
{
    cutting_tri_mesh->particles.Add_Elements(4);
    cutting_tri_mesh->Update_Number_Nodes();
    cutting_tri_mesh->mesh.elements.Append(I3(0, 1, 2));
    cutting_tri_mesh->mesh.elements.Append(I3(1, 2, 3));
    RANDOM_NUMBERS<T> r;
    for (int i = 0; i < 1e4; ++i) {
        if (sim_volume) {
            delete sim_volume;
        }
        sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
        sim_volume->particles.Store_Velocity();
        
        sim_volume->particles.Add_Elements(6);
        sim_volume->Update_Number_Nodes();
        sim_volume->particles.X(0) = TV(-0.5, -0.5, 0.0);
        sim_volume->particles.X(1) = TV(0.5, -0.5, 0.0);
        sim_volume->particles.X(2) = TV(0.5, 0.5, 0.0);
        sim_volume->particles.X(3) = TV(-0.5, 0.5, 0.0);
        sim_volume->particles.X(4) = TV(0.0, 0.0, 0.0);
        sim_volume->particles.X(5) = TV(0.0, 0.0, 0.5);
        for (int i = 0; i < 4; i++) {
            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>((i+1)%4, i, 4, 5));
        }
        
        if (mcut) {
            delete mcut;
        }
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
        
        T rr = (r.Get_Number() - 0.5) * 2;
        cout << "rr" << i << ":" << rr << endl;
        cutting_tri_mesh->particles.X(0) = TV(1, rr, -1);
        cutting_tri_mesh->particles.X(1) = TV(1, rr, 1);
        cutting_tri_mesh->particles.X(2) = TV(-1, -rr, -1);
        cutting_tri_mesh->particles.X(3) = TV(-1, -rr, 1);
        mcut->Cut(*cutting_tri_mesh);
        ARRAY<int> l;
        sim_volume->mesh.boundary_mesh->Identify_Connected_Components(l);
        if (l.Max() != 2) {
            break;
        }
    }
}

void initialize_volume4()
{
    sim_volume->particles.Add_Elements(10);
    sim_volume->Update_Number_Nodes();
    sim_volume->particles.X(0) = TV(-0.5, 0.0, -0.5);
    sim_volume->particles.X(1) = TV(0.0, 0.0, -0.5);
    sim_volume->particles.X(2) = TV(0.5, 0.0, -0.5);
    sim_volume->particles.X(3) = TV(0.5, 0.0, 0.0);
    sim_volume->particles.X(4) = TV(0.5, 0.0, 0.5);
    sim_volume->particles.X(5) = TV(0.0, 0.0, 0.5);
    sim_volume->particles.X(6) = TV(-0.5, 0.0, 0.5);
    sim_volume->particles.X(7) = TV(-0.5, 0.0, 0.0);
    sim_volume->particles.X(8) = TV(0.0, 0.0, 0.0);
    sim_volume->particles.X(9) = TV(0.0, 0.5, 0.0);
    
    for (int i = 0; i < 2; i++) {
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>((i+1)%8, i, 8, 9));
    }
}

void initialize_volume5()
{
    cutting_tri_mesh->particles.Add_Elements(4);
    cutting_tri_mesh->Update_Number_Nodes();
    cutting_tri_mesh->mesh.elements.Append(I3(0, 1, 2));
    cutting_tri_mesh->mesh.elements.Append(I3(1, 2, 3));
    RANDOM_NUMBERS<T> r;
    for (int i = 0; i < 1e4; ++i) {
        if (sim_volume) {
            delete sim_volume;
        }
        sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
        sim_volume->particles.Store_Velocity();
        
        sim_volume->particles.Add_Elements(7);
        sim_volume->Update_Number_Nodes();
        sim_volume->particles.X(0) = TV(0.0, 0.0, 0.0);
        sim_volume->particles.X(1) = TV(0.5, 0.0, 0.0);
        sim_volume->particles.X(2) = TV(0.0, 0.0, 0.5);
        sim_volume->particles.X(3) = TV(0.0, 0.5, 0.0);
        sim_volume->particles.X(4) = TV(0.0, -0.5, 0.0);
        sim_volume->particles.X(5) = TV(-0.5, 0.0, 0.0);
        sim_volume->particles.X(6) = TV(0.0, 0.0, -0.5);
        //sim_volume->particles.X(7) = TV(0.5, 0.5, 0.5);
        
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 2, 1, 3));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 2, 4));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 5, 2, 3));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 2, 5, 4));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 6, 3));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 6, 1, 4));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 6, 5, 3));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 5, 6, 4));
        //sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(1, 2, 3, 7));
        
        if (mcut) {
            delete mcut;
        }
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
        
        T rr = (r.Get_Number() - 0.5) * 2;
        cout << "rr" << i << ":" << rr << endl;
        TV pert(r.Get_Number(), r.Get_Number(), r.Get_Number());
        pert *= 0.1;
        cutting_tri_mesh->particles.X(0) = TV(1, rr, -1) + pert;
        cutting_tri_mesh->particles.X(1) = TV(1, rr, 1) + pert;
        cutting_tri_mesh->particles.X(2) = TV(-1, -rr, -1) + pert;
        cutting_tri_mesh->particles.X(3) = TV(-1, -rr, 1) + pert;
        mcut->Cut(*cutting_tri_mesh);
        ARRAY<int> l;
        sim_volume->mesh.boundary_mesh->Identify_Connected_Components(l);
        if (l.Max() != 2) {
            break;
        }
    }
}

void initialize_volume6()
{
    sim_volume->particles.Add_Elements(6);
    sim_volume->particles.X(0) = TV(0.0, 0.0, 0.0);
    sim_volume->particles.X(1) = TV(0.5, 0.0, 0.0);
    sim_volume->particles.X(2) = TV(0.0, 0.0, 0.5);
    sim_volume->particles.X(3) = TV(0.0, 0.5, 0.0);
    sim_volume->particles.X(4) = TV(-0.5, 0.0, 0.0);
    sim_volume->particles.X(5) = TV(0.0, 0.0, -0.5);

    sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 2, 1, 3));
    sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 4, 2, 3));
    sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 5, 3));
    sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 5, 4, 3));
}

void initialize_volume7()
{
    cutting_tri_mesh->particles.Add_Elements(4);
    cutting_tri_mesh->Update_Number_Nodes();
    cutting_tri_mesh->mesh.elements.Append(I3(0, 1, 2));
    cutting_tri_mesh->mesh.elements.Append(I3(1, 2, 3));
    RANDOM_NUMBERS<T> r;
    for (int i = 0; i < 1e4; ++i) {
        if (sim_volume) {
            delete sim_volume;
        }
        sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
        sim_volume->particles.Store_Velocity();
        
        sim_volume->particles.Add_Elements(6);
        sim_volume->Update_Number_Nodes();
        sim_volume->particles.X(0) = TV(0, 0.5, 0.0);
        sim_volume->particles.X(1) = TV(0.5, 0, 0.0);
        sim_volume->particles.X(2) = TV(0, -0.5, 0.0);
        sim_volume->particles.X(3) = TV(-0.5, 0, 0.0);
        sim_volume->particles.X(4) = TV(0, 0, 0.0);
        sim_volume->particles.X(5) = TV(0.0, 0.0, -0.5);
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(1, 2, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(2, 3, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(3, 0, 4, 5));
        
        if (mcut) {
            delete mcut;
        }
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
    
        T rr = (r.Get_Number() - 0.5) * 2;
        cout << "rr" << i << ":" << rr << endl;
        TV pert(r.Get_Number(), r.Get_Number(), r.Get_Number());
        pert *= 0.1;
        
        cutting_tri_mesh->particles.X(0) = TV(1, rr, -1) + pert;
        cutting_tri_mesh->particles.X(1) = TV(1, rr, 1) + pert;
        cutting_tri_mesh->particles.X(2) = TV(-1, -rr, -1) + pert;
        cutting_tri_mesh->particles.X(3) = TV(-1, -rr, 1) + pert;
        mcut->Cut(*cutting_tri_mesh);
        ARRAY<int> l;
        sim_volume->mesh.boundary_mesh->Identify_Connected_Components(l);
        if (l.Max() != 2) {
            break;
        }
    }
}

void initialize_volume8()
{
    cutting_tri_mesh->particles.Add_Elements(4);
    cutting_tri_mesh->Update_Number_Nodes();
    cutting_tri_mesh->mesh.elements.Append(I3(0, 1, 2));
    cutting_tri_mesh->mesh.elements.Append(I3(1, 2, 3));
    RANDOM_NUMBERS<T> r;
    for (int i = 0; i < 1e6; ++i) {
        if (sim_volume) {
            delete sim_volume;
        }
        sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
        sim_volume->particles.Store_Velocity();
        
        sim_volume->particles.Add_Elements(7);
        sim_volume->Update_Number_Nodes();
        sim_volume->particles.X(0) = TV(0, 0.5, 0.0);
        sim_volume->particles.X(1) = TV(0.5, 0, 0.0);
        sim_volume->particles.X(2) = TV(0, -0.5, 0.0);
        sim_volume->particles.X(3) = TV(-0.5, 0, 0.0);
        sim_volume->particles.X(4) = TV(0, 0, 0.0);
        sim_volume->particles.X(5) = TV(0.0, 0.0, -0.5);
        sim_volume->particles.X(6) = TV(0.0, 0.0, 0.5);
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(1, 2, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(2, 3, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(3, 0, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 4, 6));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(1, 2, 4, 6));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(2, 3, 4, 6));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(3, 0, 4, 6));
        
        if (mcut) {
            delete mcut;
        }
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
        
        T rr = (r.Get_Number() - 0.5) * 2;
        cout << "rr" << i << ":" << rr << endl;
        TV pert(r.Get_Number(), r.Get_Number(), r.Get_Number());
        pert *= 0;
        
        cutting_tri_mesh->particles.X(0) = TV(1, rr, -1) + pert;
        cutting_tri_mesh->particles.X(1) = TV(1, rr, 1) + pert;
        cutting_tri_mesh->particles.X(2) = TV(-1, -rr, -1) + pert;
        cutting_tri_mesh->particles.X(3) = TV(-1, -rr, 1) + pert;
        mcut->Cut(*cutting_tri_mesh);
        ARRAY<int> l;
        sim_volume->mesh.boundary_mesh->Identify_Connected_Components(l);
        if (l.Max() < 2) {
            break;
        }
    }
}

void initialize_volume9()//cubes
{
    cutting_tri_mesh->particles.Add_Elements(4);
    cutting_tri_mesh->Update_Number_Nodes();
    cutting_tri_mesh->mesh.elements.Append(I3(0, 1, 2));
    cutting_tri_mesh->mesh.elements.Append(I3(1, 2, 3));
    RANDOM_NUMBERS<T> r;
    for (int i = 0; i < 1e6; ++i) {
        if (sim_volume) {
            delete sim_volume;
        }
        sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
        sim_volume->particles.Store_Velocity();
        
        sim_volume->particles.Add_Elements(7);
        sim_volume->Update_Number_Nodes();
        sim_volume->particles.X(0) = TV(0, 0.5, 0.0);
        sim_volume->particles.X(1) = TV(0.5, 0, 0.0);
        sim_volume->particles.X(2) = TV(0, -0.5, 0.0);
        sim_volume->particles.X(3) = TV(-0.5, 0, 0.0);
        sim_volume->particles.X(4) = TV(0, 0, 0.0);
        sim_volume->particles.X(5) = TV(0.0, 0.0, -0.5);
        sim_volume->particles.X(6) = TV(0.0, 0.0, 0.5);
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(1, 2, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(2, 3, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(3, 0, 4, 5));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0, 1, 4, 6));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(1, 2, 4, 6));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(2, 3, 4, 6));
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(3, 0, 4, 6));
        
        if (mcut) {
            delete mcut;
        }
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
        
        T rr = (r.Get_Number() - 0.5) * 2;
        cout << "rr" << i << ":" << rr << endl;
        TV pert(r.Get_Number(), r.Get_Number(), r.Get_Number());
        pert *= 0;
        
        cutting_tri_mesh->particles.X(0) = TV(1, rr, -1) + pert;
        cutting_tri_mesh->particles.X(1) = TV(1, rr, 1) + pert;
        cutting_tri_mesh->particles.X(2) = TV(-1, -rr, -1) + pert;
        cutting_tri_mesh->particles.X(3) = TV(-1, -rr, 1) + pert;
        mcut->Cut(*cutting_tri_mesh);
        ARRAY<int> l;
        sim_volume->mesh.boundary_mesh->Identify_Connected_Components(l);
        if (l.Max() < 2) {
            break;
        }
    }
}

void initialize_cubes()
{
    cubes = 1;
    int width = 15;
    int height =15;
    int depth = 1;
    T low = -0.5;
    T high = 0.5;
    T left = -0.5;
    T right = 0.5;
    T far = 0.05;
    T near = -0.05;
    T lw = (right-left)/width;
    T lh = (high-low)/height;
    T ld = (far-near)/depth;
    int offset = (width+1)*(height+1);
    
    sim_volume->particles.Add_Elements(offset*(depth+1));
    sim_volume->Update_Number_Nodes();
    int node_index = 0;
    
    for (int i = 0; i < height+1; i++) {
        for (int j = 0; j < width+1; j++) {
            for (int k = 0; k < depth+1; k++) {
                sim_volume->particles.X(node_index+offset*k) = TV(left+j*lw, low+i*lh,near+k*ld);
            }
            node_index++;
        }
    }  

    int tet_index = 0;
    for (int k = 0; k < depth; k++) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
    //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset, tet_index+1, tet_index, tet_index+width+1));
    //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+2, tet_index+width+1, tet_index+1, tet_index+offset));
    //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+width+2, tet_index+offset+width+1, tet_index+offset));
    //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+1, tet_index+offset, tet_index+offset+1, tet_index+width+2));
    //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset+width+1, tet_index+width+2, tet_index+offset+width+2, tet_index+offset+1));
    //            sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset+1, tet_index+width+2, tet_index+offset, tet_index+offset+width+1));
                
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+offset+1, tet_index+1, tet_index, tet_index+width+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index, tet_index+offset, tet_index+offset+1, tet_index+width+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+offset, tet_index+offset+1, tet_index+offset+width+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+1, tet_index+width+1, tet_index+offset+1, tet_index+width+2));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+width+2, tet_index+offset+width+2, tet_index+offset+1));
                sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(tet_index+width+1, tet_index+width+2+offset, tet_index+offset+width+1, tet_index+offset+1));
                tet_index++;
            }
            tet_index++;
        }
        tet_index += (width+1);    
    }
}

void Initialize_Sphere_Levelset()
{
    mcut->levelsets.Resize(1);
    int width = 5;
    int height = 5;
    int depth = 5;
    T low = -0.5;
    T high = 0.5;
    T left = -0.5;
    T right = 0.5;
    T far = 0.5;
    T near = -0.5;
    T lw = (right-left)/width;
    T lh = (high-low)/height;
    T ld = (far-near)/depth;
    int offset = (width+1)*(height+1);
    
    mcut->levelsets(0).Resize(offset*(depth+1));
    int node_index = 0;
    T r_square = 0.421 * 0.421;
    
    for (int i = 0; i < height+1; i++) {
        for (int j = 0; j < width+1; j++) {
            for (int k = 0; k < depth+1; k++) {
                mcut->levelsets(0)(node_index+offset*k) = TV(left+j*lw, low+i*lh,near+k*ld).Magnitude_Squared()-r_square;
            }
            node_index++;
        }
    }  
}

void Initialize(bool reinitialize_cutting_mesh)
{
    drawing_cutting = 1;
    sim_volume = TETRAHEDRALIZED_VOLUME<T>::Create();
    sim_volume->particles.Store_Velocity();
    
    if(reinitialize_cutting_mesh){
        cutting_tri_mesh = TRIANGULATED_SURFACE<T>::Create();   
    }
    
    if(argc1 == 1) {
        initialize_volume1();
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
    }
    else {
        //csg
        if (0) {
            const std::string filename(argv1[1]);
            FILE_UTILITIES::Read_From_File<T>(filename, *sim_volume);
            //sim_volume->Initialize_Cube_Mesh_And_Particles(GRID<TV>(PhysBAM::VECTOR<int,3>(208, 128, 68),RANGE<TV>(TV(-0.52,-0.32,-0.17),TV(0.52,0.32,0.17))));
            mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
            
            const std::string surface_filename(argv1[2]);
            TRIANGULATED_SURFACE<float> *ts_float = TRIANGULATED_SURFACE<float>::Create();
            FILE_UTILITIES::Read_From_File<float>(surface_filename, *ts_float);
            
            cutting_tri_mesh->particles.Add_Elements(ts_float->particles.X.m);
            cutting_tri_mesh->Update_Number_Nodes();
            for (int i = 0; i < ts_float->particles.X.m; ++i) {
                for (int j = 0; j < 3; ++j) {
                    cutting_tri_mesh->particles.X(i)(j) = ts_float->particles.X(i)(j);
                }
            }
            cutting_tri_mesh->mesh.elements = ts_float->mesh.elements;
            T w = 0.14;
            Fit_In_Box<TV>(cutting_tri_mesh->particles.X, RANGE<TV>(TV(-w,-w,-w),TV(w,w,w)));
            for (int i = 0; i < ts_float->particles.X.m; ++i) {
                cutting_tri_mesh->particles.X(i)(0) = -cutting_tri_mesh->particles.X(i)(0)-0.1;
                cutting_tri_mesh->particles.X(i)(1) = cutting_tri_mesh->particles.X(i)(1)+0.05;
                cutting_tri_mesh->particles.X(i)(2) = cutting_tri_mesh->particles.X(i)(2)+0.1;
            }
            
            mcut->Cut(*cutting_tri_mesh);
            FILE_UTILITIES::Write_To_File<T>("cow_cut_by_bunny.tet.gz", *(mcut->volume));
            
            mcut->volume->Update_Number_Nodes();
            mcut->volume->mesh.Identify_Face_Connected_Components(labels);
            ARRAY<I4> old_mesh, new_mesh;
            old_mesh = mcut->volume->mesh.elements;
            int b = 5;
            for (int i = 0; i < old_mesh.m; ++i) {
                if (labels(i) != b) {
                    new_mesh.Append(old_mesh(i));
                }
            }
            mcut->volume->mesh.elements = new_mesh;
            FILE_UTILITIES::Write_To_File<T>("cow_without_bunny.tet.gz", *(mcut->volume));
            
            new_mesh.Remove_All();
            for (int i = 0; i < old_mesh.m; ++i) {
                if (labels(i) == b) {
                    new_mesh.Append(old_mesh(i));
                }
            }
            mcut->volume->mesh.elements = new_mesh;
            FILE_UTILITIES::Write_To_File<T>("bunny.tet.gz", *(mcut->volume));

            sim_volume->Update_Number_Nodes();
            sim_volume->mesh.Initialize_Boundary_Mesh(); //cout << "sim boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
            sim_volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
            
            mcut->volume->mesh.elements = old_mesh;
            mcut->volume->Update_Number_Nodes();
            mcut->volume->mesh.Initialize_Boundary_Mesh(); //cout << "cutting boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
            mcut->volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
            mcut->volume->mesh.Identify_Face_Connected_Components(labels);
            return;
        }
        
        
        const std::string filename(argv1[1]);
        TETRAHEDRALIZED_VOLUME<float> *sim_volume_float;
        sim_volume_float = TETRAHEDRALIZED_VOLUME<float>::Create();
        FILE_UTILITIES::Read_From_File<float>(filename, *sim_volume_float);
        sim_volume->particles.Add_Elements(sim_volume_float->particles.X.m);
        sim_volume->Update_Number_Nodes();
        for (int i = 0; i < sim_volume_float->particles.X.m; ++i) {
            for (int j = 0; j < 3; ++j) {
                sim_volume->particles.X(i)(j) = sim_volume_float->particles.X(i)(j);
            }
        }
        sim_volume->mesh.elements = sim_volume_float->mesh.elements;
        {
            ARRAY<int> ll;
            sim_volume->mesh.Identify_Face_Connected_Components(ll);
            cout << ll.Max() << " CCs\n";
        }
        Fit_In_Box<TV>(sim_volume->particles.X, RANGE<TV>(TV(-0.6,-0.6,-0.6),TV(0.6,0.6,0.6)));
        
        //for reading mesh in double
//        FILE_UTILITIES::Read_From_File<T>(filename, *sim_volume);
        
        
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
        
        if(argc1 == 2) {
            if (writing_interaction_to_file){
                string interaction_file_name = "interactions/interaction2.txt";
                interaction_file.open(interaction_file_name.c_str());
            }
            
            if (0) {
                if (1) {
                    int l = 11;//number of cuts
                    int n = 5;//refine curve
                    int f1 = 1, f2 = 181;
                    int f = f2 - f1;
                    
                    T r = 0.58;
                    T w = sqrt(sqr(0.6)-sqr(r))*1.1;
                    cout << "half width " << w / 1.1 << endl;
                    
                    T total_theta = 4 * pi;
                    T init_theta = -pi * 3 / 4;
                    T dtheta_cut = total_theta / f / n;
                    T theta = init_theta;
                    
                    T init_phi = asin(-0.5 / 0.6);
                    T final_phi = asin(0.5 / 0.6);
                    T dphi_cut = (final_phi - init_phi) / f / n;
                    T phi = init_phi;
                    
                    T pt = (final_phi - init_phi) / total_theta;
                    int i = 0;
                    
                    //get L_input from complete cutting mesh
                    T L=-1;
                    bool uniform_L = 1;
                    if (uniform_L) {
                        TRIANGULATED_SURFACE<T>* ts = TRIANGULATED_SURFACE<T>::Create();
                        ts->particles.Add_Elements((n+1)*3*l);
                        ts->Update_Number_Nodes();
                        while (i < l) {
                            if (i == 0) {
                                TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                                TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                                TV dp = p.Cross(t).Normalized() * w;
                                for (int j = 0; j < n+1; ++j) {
                                    ts->particles.X(3*j) = p - dp - t * w * (n - j) / n;
                                    ts->particles.X(3*j+1) = p - t * w * (n - j) / n;
                                    ts->particles.X(3*j+2) = p + dp - t * w * (n - j) / n;
                                }
                            }
                            else if (i < l - 1) {
                                for (int j = 0; j < n+1; ++j) {
                                    TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                                    TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt);
                                    TV dp = p.Cross(t).Normalized() * w;
                                    ts->particles.X(3*j+(n + 1) * 3 * i) = p - dp;
                                    ts->particles.X(3*j+1+(n + 1) * 3 * i) = p;
                                    ts->particles.X(3*j+2+(n + 1) * 3 * i) = p + dp;
                                    if (j != n) {
                                        theta += dtheta_cut;
                                        phi += dphi_cut;
                                    }
                                }
                            }
                            else if (i == l - 1) {
                                TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                                TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                                TV dp = p.Cross(t).Normalized() * w;
                                for (int j = 0; j < n+1; ++j) {
                                    ts->particles.X(3*j+(n + 1) * 3 * i) = p - dp + t * w * j / n;
                                    ts->particles.X(3*j+1+(n + 1) * 3 * i) = p + t * w * j / n;
                                    ts->particles.X(3*j+2+(n + 1) * 3 * i) = p + dp + t * w * j / n;
                                }
                            }
                            ++i;
                        }
                        for (int i = 0; i < l; ++i) {
                            for (int j = 0; j < n; ++j) {
                                for (int k = 0; k < 2; ++k) {
                                    int s = 3 * j + k + (n + 1) * 3 * i;
                                    ts->mesh.elements.Append(I3(s, s+3, s+1));
                                    ts->mesh.elements.Append(I3(s+1, s+3, s+4));
                                }
                            }
                        }
                        
                        T L_tv=0,L_ts=0;
                        for(int i=0;i<sim_volume->mesh.elements.m;i++)
                            L_tv=std::max(L_tv, RANGE<TV>::Bounding_Box(sim_volume->particles.X.Subset(sim_volume->mesh.elements(i))).Edge_Lengths().Max());
                        for(int i=0;i<ts->mesh.elements.m;i++)
                            L_ts=std::max(L_ts, RANGE<TV>::Bounding_Box(ts->particles.X.Subset(ts->mesh.elements(i))).Edge_Lengths().Max());
                        L=L_tv+L_ts;
                    }
                    
                    
                    
                    cutting_tri_mesh->particles.Add_Elements((n+1)*3);
                    cutting_tri_mesh->Update_Number_Nodes();
                    for (int j = 0; j < n; ++j) {
                        for (int k = 0; k < 2; ++k) {
                            int s = 3 * j + k;
                            cutting_tri_mesh->mesh.elements.Append(I3(s, s+3, s+1));
                            cutting_tri_mesh->mesh.elements.Append(I3(s+1, s+3, s+4));
                        }
                    }
                    
                    i = 0;
                    theta = init_theta;
                    phi = init_phi;
                    while (i < l) {
                        if (i == 0) {
                            TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                            TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                            TV dp = p.Cross(t).Normalized() * w;
                            for (int j = 0; j < n+1; ++j) {
                                cutting_tri_mesh->particles.X(3*j) = p - dp - t * w * (n - j) / n;
                                cutting_tri_mesh->particles.X(3*j+1) = p - t * w * (n - j) / n;
                                cutting_tri_mesh->particles.X(3*j+2) = p + dp - t * w * (n - j) / n;
                            }
                            mcut->Cut(*cutting_tri_mesh, false, true, L);
                        }
                        else if (i < l - 1) {
                            for (int j = 0; j < n+1; ++j) {
                                TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                                TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt);
                                TV dp = p.Cross(t).Normalized() * w;
                                cutting_tri_mesh->particles.X(3*j) = p - dp;
                                cutting_tri_mesh->particles.X(3*j+1) = p;
                                cutting_tri_mesh->particles.X(3*j+2) = p + dp;
                                if (j != n) {
                                    theta += dtheta_cut;
                                    phi += dphi_cut;
                                }
                            }
                            mcut->Cut(*cutting_tri_mesh, false, true, L);
                        }
                        else if (i == l - 1) {
                            TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                            TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                            TV dp = p.Cross(t).Normalized() * w;
                            for (int j = 0; j < n+1; ++j) {
                                cutting_tri_mesh->particles.X(3*j) = p - dp + t * w * j / n;
                                cutting_tri_mesh->particles.X(3*j+1) = p + t * w * j / n;
                                cutting_tri_mesh->particles.X(3*j+2) = p + dp + t * w * j / n;
                            }
                            mcut->Cut(*cutting_tri_mesh, true, true, L);
                        }
                        ++i;
                    }
                }
                else {
                    int f1 = 1, f2 = 181;
                    int f = f2 - f1;
                    int n = 5;//refine curve
                    
                    T r = 0.58;
                    T w = sqrt(sqr(0.6)-sqr(r))*1.1;
                    cout << "half width " << w / 1.1 << endl;
                    
                    T total_theta = 4 * pi;
                    T init_theta = -pi * 3 / 4;
                    T dtheta_cut = total_theta / f / n;
                    T theta = init_theta;
                    
                    T init_phi = asin(-0.5 / 0.6);
                    T final_phi = asin(0.5 / 0.6);
                    T dphi_cut = (final_phi - init_phi) / f / n;
                    T phi = init_phi;
                    
                    T pt = (final_phi - init_phi) / total_theta;
                    
                    int l = 11;
                    cutting_tri_mesh->particles.Add_Elements((n+1)*3*l);
                    cutting_tri_mesh->Update_Number_Nodes();
                    for (int i = 0; i < l; ++i) {
                        for (int j = 0; j < n; ++j) {
                            for (int k = 0; k < 2; ++k) {
                                int s = 3 * j + k + (n + 1) * 3 * i;
                                cutting_tri_mesh->mesh.elements.Append(I3(s, s+3, s+1));
                                cutting_tri_mesh->mesh.elements.Append(I3(s+1, s+3, s+4));
                            }
                        }
                    }
                    
                    int i = 0;
                    while (i < l) {
                        if (i == 0) {
                            TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                            TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                            TV dp = p.Cross(t).Normalized() * w;
                            for (int j = 0; j < n+1; ++j) {
                                cutting_tri_mesh->particles.X(3*j) = p - dp - t * w * (n - j) / n;
                                cutting_tri_mesh->particles.X(3*j+1) = p - t * w * (n - j) / n;
                                cutting_tri_mesh->particles.X(3*j+2) = p + dp - t * w * (n - j) / n;
                            }
                        }
                        else if (i < l - 1) {
                            for (int j = 0; j < n+1; ++j) {
                                TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                                TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt);
                                TV dp = p.Cross(t).Normalized() * w;
                                cutting_tri_mesh->particles.X(3*j+(n + 1) * 3 * i) = p - dp;
                                cutting_tri_mesh->particles.X(3*j+1+(n + 1) * 3 * i) = p;
                                cutting_tri_mesh->particles.X(3*j+2+(n + 1) * 3 * i) = p + dp;
                                if (j != n) {
                                    theta += dtheta_cut;
                                    phi += dphi_cut;
                                }
                            }
                        }
                        else if (i == l - 1) {
                            TV p = TV(sin(phi), cos(phi) * cos(theta), cos(phi) * sin(theta)) * r;
                            TV t = TV(cos(phi), -sin(phi) * cos(theta) - cos(phi) * sin(theta) / pt, -sin(phi) * sin(theta) + cos(phi) * cos(theta) / pt).Normalized();
                            TV dp = p.Cross(t).Normalized() * w;
                            for (int j = 0; j < n+1; ++j) {
                                cutting_tri_mesh->particles.X(3*j+(n + 1) * 3 * i) = p - dp + t * w * j / n;
                                cutting_tri_mesh->particles.X(3*j+1+(n + 1) * 3 * i) = p + t * w * j / n;
                                cutting_tri_mesh->particles.X(3*j+2+(n + 1) * 3 * i) = p + dp + t * w * j / n;
                            }
                            mcut->Cut(*cutting_tri_mesh, true, true);
                        }
                        ++i;
                    }
                }
            }
            if (0) {//test whether cutting result is correct
                T vo = sim_volume->Total_Volume();
                int ii;
                T mvo = sim_volume->Minimum_Volume(&ii);
                TETRAHEDRALIZED_VOLUME<T>* v = NULL;
                copy(v, sim_volume);
                cutting_tri_mesh->particles.Add_Elements(4);
                cutting_tri_mesh->Update_Number_Nodes();
                cutting_tri_mesh->mesh.elements.Append(I3(0, 1, 2));
                cutting_tri_mesh->mesh.elements.Append(I3(1, 2, 3));
                RANDOM_NUMBERS<T> r;
                for (int i = 0; i < 1e4; ++i) {
                    copy(sim_volume, v);
                    sim_volume->particles.Store_Velocity();
                    
                    if (mcut) {
                        delete mcut;
                    }
                    mcut = new MESH_CUTTING<T>(sim_volume, timestep, cutting_ratio, true);
                    
                    T rr = (r.Get_Number() - 0.5) * 20;
                    cout << "rr" << i << ":" << rr << endl;
                    TV pert(r.Get_Number(), r.Get_Number(), r.Get_Number());
                    pert *= 0.1;
                    
                    cutting_tri_mesh->particles.X(0) = TV(1, rr, -1) + pert;
                    cutting_tri_mesh->particles.X(1) = TV(1, rr, 1) + pert;
                    cutting_tri_mesh->particles.X(2) = TV(-1, -rr, -1) + pert;
                    cutting_tri_mesh->particles.X(3) = TV(-1, -rr, 1) + pert;
                    mcut->Cut(*cutting_tri_mesh);
                    ARRAY<int> l;
                    mcut->Connected_Components(mcut->volume, l);
                    int cc = l.Max();
                    cout << cc << " CCs after cut\n";
                    ARRAY<T> vcc(cc);
                    for (int j = 0; j < l.m; ++j) {
                        I4 e = mcut->volume->mesh.elements(j);
                        vcc(l(j)-1) += TETRAHEDRON<T>(mcut->volume->particles.X(e(0)),
                                                    mcut->volume->particles.X(e(1)),
                                                    mcut->volume->particles.X(e(2)),
                                                    mcut->volume->particles.X(e(3))).Signed_Volume();
                    }
                    vcc.Sort();
                    int ecc = 0;
                    for (int j = 0; j < vcc.m; ++j) {
                        if (vcc(j) > mvo) {
                            ++ecc;
                        }
                    }
                    T s = vcc.Sum();
                    cout << ecc << " effective CCs after cut" << endl;
                    cout << vcc << " sums to " << s << ", while original volume is " << vo << endl;
                    if (ecc != 2 || fabs(s - vo) > 1e-6) {
                        break;
                    }
                }
            }
        }
        if(argc1 == 3) {
            const std::string surface_filename(argv1[2]);
//            TRIANGULATED_SURFACE<float> *ts_float = TRIANGULATED_SURFACE<float>::Create();
//            FILE_UTILITIES::Read_From_File<float>(surface_filename, *ts_float);
            
            
//            cutting_tri_mesh->particles.Add_Elements(ts_float->particles.X.m);
//            cutting_tri_mesh->Update_Number_Nodes();
//            for (int i = 0; i < ts_float->particles.X.m; ++i) {
//                for (int j = 0; j < 3; ++j) {
//                    cutting_tri_mesh->particles.X(i)(j) = ts_float->particles.X(i)(j);
//                }
//            }
//            cutting_tri_mesh->mesh.elements = ts_float->mesh.elements;
//            Fit_In_Box<TV>(cutting_tri_mesh->particles.X, RANGE<TV>(TV(-0.4,-0.4,-0.4),TV(0.4,0.4,0.4)));
            
            FILE_UTILITIES::Read_From_File<T>(surface_filename, *cutting_tri_mesh);
            mcut->Cut(*cutting_tri_mesh);
        }
    }
    sim_volume->Update_Number_Nodes();
    sim_volume->mesh.Initialize_Boundary_Mesh(); //cout << "sim boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
    sim_volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
    
    mcut->volume->Update_Number_Nodes();
    mcut->volume->mesh.Initialize_Boundary_Mesh(); //cout << "cutting boundary elements:" << sim_volume->mesh.boundary_mesh->elements.m << endl;
    mcut->volume->mesh.boundary_mesh->Initialize_Segment_Mesh();
    mcut->volume->mesh.Identify_Face_Connected_Components(labels);
    
    if (recording)
        VS::start_timer();
}

void display(){}
int main(int argc, char **argv)
{ 
    argc1 = argc;
    argv1 = argv;
    glutInit( &argc, argv );
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
    glutInitWindowSize( window_width, window_height );
    string window_name = "cutting";
    glutCreateWindow( window_name.c_str() );
    glClearColor(1, 1, 1, 1);
    glEnable( GL_DEPTH_TEST );
    
    Initialize(1);
    
    //glutTimerFunc(timestep, time_func, 0);
    glutSpecialFunc( SpecialKey );
    glutKeyboardFunc( Key );
    glutMouseFunc( mouse );
    glutMotionFunc( motion );
    glutDisplayFunc( Render );
    glutReshapeFunc(Reshape);
    
    glutMainLoop();
    
    return 0;
}
