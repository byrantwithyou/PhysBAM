//
//  mesh_cutting_test.cpp
//
//
//  Created by Yuting Wang on 5/24/12.
//  Copyright (c) 2012 __Yuting Wang__. All rights reserved.
//

#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/ROTATION.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Geometry/Topology/SEGMENT_MESH.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

#include <fstream>
#include <sstream>

#include "DEFORMABLE_OBJECTS.h"
#include "mesh_cutting_subd_old.h"

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
#include <GL/glew.h>
#endif  // __APPLE__
namespace VS{
    void read_tsc(){__asm__("rdtsc");}
    struct timeval starttime,stoptime;
    void start_timer(){gettimeofday(&starttime,NULL);read_tsc();}
    void stop_timer(){gettimeofday(&stoptime,NULL);read_tsc();}
    double get_time(){return (double)stoptime.tv_sec-(double)starttime.tv_sec+(double)1e-6*(double)stoptime.tv_usec-(double)1e-6*(double)starttime.tv_usec;}
}

using namespace PhysBAM;
using namespace std;

typedef double T;
typedef PhysBAM::VECTOR<T, 3> TV;
typedef PhysBAM::MATRIX<T, 3> TM;
typedef PhysBAM::VECTOR<int, 2> I2;
typedef PhysBAM::VECTOR<int, 3> I3;
typedef PhysBAM::VECTOR<int, 4> I4;

//global variables
int argc1;
char **argv1;

//mouse motion modes
bool cutting = 1;
bool dragging = 0;
bool rotating = 0;
bool translating = 0;
bool lighting_on = 0;

bool drawing_cutting = 1;
int dragging_id = -1;
T K = 50000;

PhysBAM::VECTOR<T,2> starting_rotation;
PhysBAM::VECTOR<T,2> starting_dragging;


bool cubes = 0;
int picked_cc_id = -1;
TV starting_position, end_position;
T timestep = 30;
T rotate_speed = 10/57.;
T intrude = -2;

int cutting_surface_id = 0;
fstream interaction_file;
bool writing_interaction_to_file = 0;
int current_frame = 0;

ARRAY<TV> cutting_curve;
MESH_CUTTING<T> *mcut;
TRIANGULATED_SURFACE<T> *cutting_tri_mesh;
TETRAHEDRALIZED_VOLUME<T> *sim_volume;

#define HEIGHT 512
#define WIDTH 512

#define DE cout<<"file "<<__FILE__<<"   line "<<__LINE__<<"  "<<&mcut->volume->particles.X<<"   "<<mcut->volume->particles.X<<endl;

int current = 0;
int ratio = 10;
bool draw_sim = 1, draw_material_edges = 1;

void Render(){
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

void time_func(int value)
{
    Render();
    glutTimerFunc(timestep, time_func, 0);
    
    if (writing_interaction_to_file) {
        interaction_file << "f " << current_frame << endl;
        current_frame++;
    }
    //cout << "*********************frame ends*********************************" << endl << endl;
}

static void SpecialKey( int key, int x, int y )
{
    VS::start_timer();
    TM r;
    switch( key ) {
        case GLUT_KEY_DOWN:
            r =TM::Rotation_Matrix_X_Axis(rotate_speed);
            break;
        case GLUT_KEY_UP:
            r = TM::Rotation_Matrix_X_Axis(-rotate_speed);
            break;
        case GLUT_KEY_RIGHT:
            r = TM::Rotation_Matrix_Y_Axis(rotate_speed);
            break;
        case GLUT_KEY_LEFT:
            r = TM::Rotation_Matrix_Y_Axis(-rotate_speed);
            break;
    }
    for (int i = 0; i < sim_volume->particles.X.m; i++) {
        sim_volume->particles.X(i) = r * sim_volume->particles.X(i);
        sim_volume->particles.V(i) = r * sim_volume->particles.V(i);
        
        for (int k = 0; k < 3; k++){
            mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
            mcut->deformable_object->Velocities()(i*3+k) = mcut->sim_volume->particles.V(i)(k);
        }
    }
    mcut->Update_Cutting_Particles();
    
    for (int i = 0; i < mcut->dragging_targets.m; i++) {
        mcut->dragging_targets(i) = r * mcut->dragging_targets(i);
    }
    
    for (int i = 0; i < mcut->my_constrained->n/3; i++){
        int fixed_node = mcut->my_constrained->operator()(3*i)/3;
        for (int k = 0; k < 3; k++){
            mcut->my_constrained_locations->operator()(3*i+k) = sim_volume->particles.X(fixed_node)(k);
        }
    }
    mcut->be->Set_Boundary_Conditions(*(mcut->my_constrained), *(mcut->my_constrained_locations));
    
    for (int i = 0; i < mcut->cutting_vertices.m; i++) {
        mcut->cutting_vertices(i) = r * mcut->cutting_vertices(i);
    }
    VS::stop_timer();
    printf("rotation time:    %f\n",VS::get_time());
}

void Initialize(bool);
void Recover_Volume();

void Translate_Volume(const TV& translation)
{
    for (int i = 0; i < sim_volume->particles.X.m; i++) {
        sim_volume->particles.X(i) += translation;
        for (int k = 0; k < 3; k++){
            mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
        }
    }
    mcut->Update_Cutting_Particles();
    
    for (int i = 0; i < mcut->dragging_targets.m; i++) {
        mcut->dragging_targets(i) += translation;
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
            return;
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
        case 'e':
            drawing_cutting = !drawing_cutting;
            break;
        case 'p':
            if(cubes){
                for (int i = 0; i < sim_volume->particles.X.m; i++) {
                    sim_volume->particles.X(i) *= 1.1;
                    for (int k = 0; k < 3; k++){
                        mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
            }
            break;
        case 'z':
            Recover_Volume();
            break;
    }
    
}

ARRAY<TV> cutting_copy;
void mouse(int button, int state, int x, int y)
{
    T xcoor = 2 * x / T(WIDTH) - 1;
    T ycoor = 1 - 2 * y / T(HEIGHT);
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
                starting_position = TV( xcoor, ycoor, intrude/2);
                cutting_curve.Append(starting_position);
            }
            else {
                starting_position = TV( xcoor, ycoor, intrude/2);
                mcut->Draw_For_Picking();
                unsigned char val;
                glReadPixels(x, HEIGHT-y, 1, 1, GL_RED, GL_UNSIGNED_BYTE, &val);
                picked_cc_id = (int)val;
                cout << "picked cc: " << picked_cc_id << endl;
            }
        }
        else if (button == 3) {
            //cout << "processing mouse rolling\n";
            if(cubes){
                for (int i = 0; i < sim_volume->particles.X.m; i++) {
                    sim_volume->particles.X(i) *= 1.1;
                    for (int k = 0; k < 3; k++){
                        mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
            }
            //            else {
            //                int i = 0;
            //                for (set<MESH_CUTTING<T>::DIRI,MESH_CUTTING<T>::DIRI_LESS>::iterator it = mcut->diri_nodes.begin(); it != mcut->diri_nodes.end(); it++) {
            //                    if((*it)(1)==1){
            //                        (*(mcut->my_constrained_locations))(3*i) += 0.05;}
            //                    else {
            //                        (*(mcut->my_constrained_locations))(3*i) -= 0.05;}
            //                    i++;
            //                }
            //                mcut->be->Set_Boundary_Conditions(*(mcut->my_constrained), *(mcut->my_constrained_locations));
            //            }
            for (int i = 0; i < mcut->cutting_vertices.m; i++) {
                mcut->cutting_vertices(i) *= 1.1;
            }
        }
        else if (button == 4) {
            //cout << "processing mouse rolling\n";
            if (cubes){
                for (int i = 0; i < sim_volume->particles.X.m; i++) {
                    sim_volume->particles.X(i) *= 0.9;
                    for (int k = 0; k < 3; k++){
                        mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
                    }
                }
                mcut->Update_Cutting_Particles();
            }
            //            else {
            //                int i = 0;
            //                for (set<MESH_CUTTING<T>::DIRI,MESH_CUTTING<T>::DIRI_LESS>::iterator it = mcut->diri_nodes.begin(); it != mcut->diri_nodes.end(); it++) {
            //                    if((*it)(1)==1){
            //                        (*(mcut->my_constrained_locations))(3*i) -= 0.05;}
            //                    else {
            //                        (*(mcut->my_constrained_locations))(3*i) += 0.05;}
            //                    i++;
            //                }
            //                mcut->be->Set_Boundary_Conditions(*(mcut->my_constrained), *(mcut->my_constrained_locations));
            //            }
            for (int i = 0; i < mcut->cutting_vertices.m; i++) {
                mcut->cutting_vertices(i) *= 0.9;
            }
        }
    }
    else if (state == GLUT_UP) {
        if (button == GLUT_LEFT_BUTTON) {
            if (translating) {
                if (picked_cc_id >= 0) {
                    mcut->Translate_CC(picked_cc_id, end_position - starting_position);
                    for (int i = 0; i < sim_volume->particles.X.m; i++) {
                        for (int k = 0; k < 3; k++){
                            mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
                        }
                    }
                }
                picked_cc_id = -1;
            }
            else if (dragging) {
                dragging_id = -1;
            }
            else if (cutting) {
                int n = cutting_curve.m;
                if (n>1){
                    int pid1 = cutting_tri_mesh->particles.X.m;
                    cutting_tri_mesh->particles.Add_Elements(2*n);
                    for (int i = 0; i < n; i++) {
                        cutting_tri_mesh->particles.X(pid1) = cutting_curve(i);
                        pid1++;
                        cutting_tri_mesh->particles.X(pid1) = cutting_curve(i) - TV(0,0,intrude);
                        pid1++;
                    }
                    pid1-=2*n;
                    cutting_tri_mesh->Update_Number_Nodes();
                    
                    n--;
                    for (int i = 0; i < n; i++) {
                        cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(pid1,pid1+1,pid1+2));
                        cutting_tri_mesh->mesh.elements.Append(PhysBAM::VECTOR<int,3>(pid1+1,pid1+3,pid1+2));
                        pid1+=2;
                    }
                    
                    //writing cutting information to file for high res simulation
                    if (writing_interaction_to_file){
                        stringstream ss; ss<<cutting_surface_id;
                        interaction_file<< "c " << cutting_surface_id << endl;
                        FILE_UTILITIES::Write_To_File<T>(string("cutting_surfaces/")+ss.str()+string(".tet.gz"), cutting_tri_mesh->particles, cutting_tri_mesh->mesh);
                        cutting_surface_id++;
                    }
                    
                    for (int i = 0; i < cutting_curve.m-1; i++) {
                        mcut->cutting_vertices.Append(cutting_curve(i));
                        mcut->cutting_vertices.Append(cutting_curve(i) - TV(0,0,intrude));
                        mcut->cutting_vertices.Append(cutting_curve(i));
                        mcut->cutting_vertices.Append(cutting_curve(i+1));
                        mcut->cutting_vertices.Append(cutting_curve(i) - TV(0,0,intrude));
                        mcut->cutting_vertices.Append(cutting_curve(i+1) - TV(0,0,intrude));
                        mcut->cutting_vertices.Append(cutting_curve(i));
                        mcut->cutting_vertices.Append(cutting_curve(i+1) - TV(0,0,intrude));
                    }
                    mcut->cutting_vertices.Append(cutting_curve(cutting_curve.m-1));
                    mcut->cutting_vertices.Append(cutting_curve(cutting_curve.m-1) - TV(0,0,intrude));
                    
                    cutting_copy = mcut->cutting_vertices;
                    mcut->Cut(*cutting_tri_mesh);
                }
                delete cutting_tri_mesh;
                cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
                cutting_curve.Remove_All();
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
        
        for (int k = 0; k < 3; k++){
            mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
            mcut->deformable_object->Velocities()(i*3+k) = mcut->sim_volume->particles.V(i)(k);
        }
    }
    mcut->Update_Cutting_Particles();
    
    for (int i = 0; i < mcut->dragging_targets.m; i++) {
        mcut->dragging_targets(i) = r * mcut->dragging_targets(i);
    }
    
    for (int i = 0; i < mcut->my_constrained->n/3; i++){
        int fixed_node = mcut->my_constrained->operator()(3*i)/3;
        for (int k = 0; k < 3; k++){
            mcut->my_constrained_locations->operator()(3*i+k) = sim_volume->particles.X(fixed_node)(k);
        }
    }
    mcut->be->Set_Boundary_Conditions(*(mcut->my_constrained), *(mcut->my_constrained_locations));
    
    for (int i = 0; i < mcut->cutting_vertices.m; i++) {
        mcut->cutting_vertices(i) = r * mcut->cutting_vertices(i);
    }
    
    starting_rotation[0] = xx;
    starting_rotation[1] = yy;
}

void motion(int x, int y)
{
    T xx = 2 * x / T(WIDTH) - 1;
    T yy = 1 - 2 * y / T(HEIGHT);
    end_position = TV( xx, yy, intrude/2);
    
    if (cutting) {
        if ((end_position - cutting_curve(cutting_curve.m-1)).Magnitude_Squared() > 0.0005)
            cutting_curve.Append(end_position);
    }
    else if(dragging) {
        if(dragging_id >= 0){
            mcut->dragging_targets(dragging_id)(0) += (xx-starting_dragging[0]);
            mcut->dragging_targets(dragging_id)(1) += (yy-starting_dragging[1]);
            starting_dragging[0] = xx;
            starting_dragging[1] = yy;
            if (writing_interaction_to_file) {
                interaction_file << "d " << dragging_id << " " << xx << " " << yy << endl;
            }
        }
    }
    else if(rotating) {
        //cout << "processing rotation" << endl;
        if (writing_interaction_to_file) {
            interaction_file << "r " << xx << " " << yy << endl;
        }
        rotate_meshes(xx, yy);
    }
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
    sim_volume->particles.X(0) = TV(-0.5, 0.0, -0.5);
    sim_volume->particles.X(1) = TV(0.5, 0.0, -0.5);
    sim_volume->particles.X(2) = TV(0.0, 0.0, 0.0);
    sim_volume->particles.X(3) = TV(0.0, 0.5, 0.0);
    sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>(0,2,1,3));
    perturb_particles(*sim_volume);
}

void initialize_volume2()
{
    sim_volume->particles.Add_Elements(5);
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
    sim_volume->particles.Add_Elements(6);
    sim_volume->particles.X(0) = TV(-0.5, 0.0, -0.5);
    sim_volume->particles.X(1) = TV(0.5, 0.0, -0.5);
    sim_volume->particles.X(2) = TV(0.5, 0.0, 0.5);
    sim_volume->particles.X(3) = TV(-0.5, 0.0, 0.5);
    sim_volume->particles.X(4) = TV(0.0, 0.0, 0.0);
    sim_volume->particles.X(5) = TV(0.0, 0.5, 0.0);
    for (int i = 0; i < 4; i++) {
        sim_volume->mesh.elements.Append(PhysBAM::VECTOR<int,4>((i+1)%4, i, 4, 5));
    }
}

void initialize_volume4()
{
    sim_volume->particles.Add_Elements(10);
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
    sim_volume->particles.Add_Elements(7);
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

void Recover_Volume()
{
    cout << cutting_tri_mesh->mesh.elements.m << endl;
    Initialize(0);
    cout << cutting_tri_mesh->mesh.elements.m << endl;
    mcut->cutting_vertices = cutting_copy;
    mcut->Cut(*cutting_tri_mesh);
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
        initialize_cubes();
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, ratio);
        mcut->Initialize_Elasticity();
    }
    else {
        const std::string filename(argv1[1]);
        TETRAHEDRALIZED_VOLUME<float> *sim_volume_float;
        sim_volume_float = TETRAHEDRALIZED_VOLUME<float>::Create();
        FILE_UTILITIES::Read_From_File<float>(filename, *sim_volume_float);
        sim_volume->particles.Add_Elements(sim_volume_float->particles.X.m);
        for (int i = 0; i < sim_volume_float->particles.X.m; ++i) {
            for (int j = 0; j < 3; ++j) {
                sim_volume->particles.X(i)(j) = sim_volume_float->particles.X(i)(j);
            }
        }
        sim_volume->mesh.elements = sim_volume_float->mesh.elements;
        Fit_In_Box<TV>(sim_volume->particles.X, RANGE<TV>(TV(-0.6,-0.6,-0.6),TV(0.6,0.6,0.6)));
        mcut = new MESH_CUTTING<T>(sim_volume, timestep, ratio);
        mcut->Initialize_Elasticity();
        if(argc1 == 2) {
            if (writing_interaction_to_file){
                string interaction_file_name = "interactions/interaction2.txt";
                interaction_file.open(interaction_file_name.c_str());
            }
        }
        if(argc1 == 3) {
            const std::string surface_filename(argv1[2]);
            FILE_UTILITIES::Read_From_File<T>(surface_filename,*cutting_tri_mesh);
            Fit_In_Box<TV>(cutting_tri_mesh->particles.X, RANGE<TV>(TV(-0.8,-0.8,-0.8),TV(0.8,0.8,0.8)));
            mcut->Cut(*cutting_tri_mesh);
        }
        else if(argc1 == 4) {
            //./virtual_surgery tet_mesh_input interaction_events_input output_directory
            const string interaction_file_name(argv1[2]);
            ifstream ifs;
            ifs.open(interaction_file_name.c_str());
            const string writing_directory(argv1[3]);
            string event;
            int frame = 0;
            K = 10000;
            while(getline(ifs, event)){
                //                if(frame == 1){
                //                    for (int i = 0; i < sim_volume->particles.X.m; i++) {
                //                        sim_volume->particles.X(i) *= 0.7;
                //                        for (int k = 0; k < 3; k++){
                //                            mcut->deformable_object->Positions()(i*3+k) = mcut->sim_volume->particles.X(i)(k);
                //                        }
                //                    }
                //                    mcut->Update_Cutting_Particles();
                //                }
                
                if (event[0] == 'f'){
                    VS::start_timer();
                    for (int i = 0; i < ratio; i++) {
                        mcut->be->Advance_One_Time_Step(*mcut, K, 1, 1);
                    }
                    VS::stop_timer();
                    printf("frame %d time:    %f\n", frame, VS::get_time());
                    for (int i = 0; i<mcut->sim_volume->particles.X.m; i++){
                        for (int k = 0; k<3; k++){
                            mcut->sim_volume->particles.X(i)(k) = mcut->deformable_object->X(i)(k);
                            mcut->sim_volume->particles.V(i)(k) = mcut->deformable_object->V(i)(k);
                        }
                    }
                    mcut->Update_Cutting_Particles();
                    mcut->Write_To_File(writing_directory, frame);
                    frame++;
                }
                else if (event[0] == 'c'){
                    string cutting_surface_file_name = string("cutting_surfaces/")+event[2]+string(".tet.gz");
                    cout << "Cutting by " + cutting_surface_file_name << endl;
                    FILE_UTILITIES::Read_From_File<T>(cutting_surface_file_name, cutting_tri_mesh->particles, cutting_tri_mesh->mesh);
                    mcut->Cut(*cutting_tri_mesh);
                    delete cutting_tri_mesh;
                    cutting_tri_mesh = new TRIANGULATED_SURFACE<T>();
                }
                else if (event[0] == 'D'){
                    stringstream ss(event);
                    string d;
                    T xx, yy;
                    ss >> d;
                    ss >> xx >> yy;
                    starting_dragging[0] = xx;
                    starting_dragging[1] = yy;
                    if (mcut->Compute_Intersection(xx, yy) == -1) {
                        cout << "not cut on high res!\n";
                        exit(1);
                    }
                }
                else if (event[0] == 'd'){
                    stringstream ss(event);
                    string d;
                    int dragging_ID;
                    T xx, yy;
                    ss >> d;
                    ss >> dragging_ID;
                    ss >> xx >> yy;
                    mcut->dragging_targets(dragging_ID)(0) += (xx-starting_dragging[0]);
                    mcut->dragging_targets(dragging_ID)(1) += (yy-starting_dragging[1]);
                    starting_dragging[0] = xx;
                    starting_dragging[1] = yy;
                }
                else if (event[0] == 'R'){
                    stringstream ss(event);
                    string r;
                    T xx, yy;
                    ss >> r;
                    ss >> xx >> yy;
                    starting_rotation[0] = xx;
                    starting_rotation[1] = yy;
                }
                else if (event[0] == 'r'){
                    stringstream ss(event);
                    string r;
                    T xx, yy;
                    ss >> r;
                    ss >> xx >> yy;
                    rotate_meshes(xx, yy);
                }
            }
            ifs.close();
            exit(0);
        }
    }
}

void display(){}
int main(int argc, char **argv)
{
    argc1 = argc;
    argv1 = argv;
    glutInit( &argc, argv );
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA);
    glutInitWindowSize( WIDTH, HEIGHT );
    string window_name = "cutting";
    glutCreateWindow( window_name.c_str() );
    glClearColor( 1.0, 1.0, 1.0, 1.0 );
    glEnable( GL_DEPTH_TEST );
    
    if (lighting_on) {
        glEnable(GL_LIGHTING);
        GLfloat dif[] = {1.f, 1.f, 1.f, 1.f};
        glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
        GLfloat lightpos[] = {0.f, 0.f, -1.f, 0.f};
        glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
        glEnable(GL_LIGHT0);
        glEnable(GL_COLOR_MATERIAL);
    }
    
    Initialize(1);
    
    glutTimerFunc(timestep, time_func, 0);
    glutSpecialFunc( SpecialKey );
    glutKeyboardFunc( Key );
    glutMouseFunc( mouse );
    glutMotionFunc( motion );
    glutDisplayFunc( display );
    
    glutMainLoop();
    
    return 0;
}
